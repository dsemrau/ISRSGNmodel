# -*- coding: utf-8 -*-
""" ISRS GN model implementation 

This module implements the function that returns the nonlinear interference 
power and coefficient for each WDM channel. This function implements the 
ISRS GN model in closed-form published in:

D. Semrau, R. I. Killey, P. Bayvel, "A Closed-Form Approximation of the
Gaussian Noise Model in the Presence of Inter-Channel Stimulated Raman
Scattering, " J. Lighw. Technol., Early Access, Jan. 2019

Author: Daniel Semrau, Eric Sillekens, R. I. Killey, P. Bayvel, Jan 2019.
"""

from numpy import abs,arange,arcsinh,arctan,isfinite,log,mean,newaxis,pi,sum,zeros

def ISRSGNmodel(
        Att, Att_bar, Cr, Pch, fi, Bch, Length, D, S, gamma, RefLambda, coherent=True,
        **P_unused):
    """
    Returns nonlinear interference power and coefficient for each WDM
    channel. This function implements the ISRS GN model in closed-form
    published in:
    
    D. Semrau, R. I. Killey, P. Bayvel, "A Closed-Form Approximation of the
    Gaussian Noise Model in the Presence of Inter-Channel Stimulated Raman
    Scattering, " J. Lighw. Technol., vol. xx, no. xx, pp.xxxx-xxxx, Jan. 2019
    
    Format:

    - channel dependent quantities have the format of a N_ch x n matrix,
      where N_ch is the number of channels slots and n is the number of spans.
    - channel independent quantities have the format of a 1 x n matrix
    - channel and span independent quantities are scalars
    
    INPUTS:
        Att: attenuation coefficient [Np/m] of channel i of span j,
            format: N_ch x n matrix
        Att_bar: attenuation coefficient (bar) [Np/m] of channel i of span j,
            format: N_ch x n matrix
        Cr[i,j]: the slope of the linear regression of the normalized Raman gain spectrum [1/W/m/Hz] of channel i of span j, 
            format: N_ch x n matrix
        
        Pch[i,j]:  the launch power [W] of channel i of span j, 
            format: N_ch x n matrix
        fi[i,j]: center frequency relative to the reference frequency (3e8/RefLambda) [Hz]
            of channel i of span j, format: N_ch x n matrix
        Bch[i,j]: the bandwidth [Hz] of channel i of span j, 
            format: N_ch x n matrix
        
        Length[j]: the span length [m] of span j,
            format: 1 x n vector
        D[j]: the dispersion coefficient [s/m^2] of span j,
            format: 1 x n vector
        S[j]: the span length [s/m^3] of span j,
            format: 1 x n vector
        gamma[j]: the span length [1/W/m] of span j,
            format: 1 x n vector
        RefLambda: is the reference wavelength (where beta2, beta3 are defined) [m],
            format: 1 x 1 vector

        coherent: boolean for coherent or incoherent NLI accumulation across multiple fiber spans 
        
    RETURNS:
        NLI: Nonlinear Interference Power[W],
            format: N_ch x 1 vector 
        eta_n: Nonlinear Interference coeffcient [1/W^2],
            format: N_ch x 1 matrix 
    """
    channels,n = fi.shape
    c = 3e8;
    
    a     = Att
    a_bar = Att_bar
    L     = Length
    P_ij  = Pch
    Ptot  = sum(P_ij,axis=0)
    
    beta2 = -D*RefLambda**2/(2*pi*c)
    beta3 = RefLambda**2/(2*pi*c)**2*(RefLambda**2*S+2*RefLambda*D)

    # Average Coherence Factor
    mean_att_i = mean(a, axis=1) #average attenuation coefficent for channel i
    mean_L     = mean(L) # average fiber length 

    if coherent:
        # closed-for formula for average coherence factor extended by dispersion slope, cf. Ref. [2, Eq. (22)]
        eps = lambda B_i, f_i, a_i: (3/10)*log(1+(6/a_i)/(mean_L*arcsinh(pi**2/2*abs( mean(beta2) + 2*pi*mean(beta3)*f_i )/a_i*B_i**2)))     
    else:
        eps = lambda B_i, f_i, a_i: 0
    
    # SPM and XPM Closed-form Formula Definition
    SPM = lambda phi_i, T_i, B_i, a, a_bar, gamma:\
            4/9*gamma**2/B_i**2*pi/(phi_i*a_bar*(2*a+a_bar)) \
              *( (T_i-a**2)/a*arcsinh(phi_i*B_i**2/a/pi) + ((a+a_bar)**2-T_i)/(a+a_bar)*arcsinh(phi_i*B_i**2/(a+a_bar)/pi) ) 
        # closed-form formula for SPM contribution, see Ref. [1, Eq. (9-10)]

    XPM = lambda Pi, Pk, phi_ik, T_k, B_i, B_k, a, a_bar, gamma:\
            32/27*sum(removenan( (Pk/Pi)**2*gamma**2 / ( B_k*phi_ik*a_bar*(2*a+a_bar) )
                                 *( (T_k-a**2)/a*arctan(phi_ik*B_i/a) 
                                    +((a+a_bar)**2-T_k)/(a+a_bar)*arctan(phi_ik*B_i/(a+a_bar)) )
                               )) 
        # closed-form formula for XPM contribution, see Ref. [1, Eq. (11)] 
    
    NLI     = []
    eta_n   = []
    eta_SPM = zeros([channels,n])
    eta_XPM = zeros([channels,n])
    for j in arange(n):
        """ Calculation of nonlinear interference (NLI) power in fiber span j """
        for i in arange(channels):
            """ Compute the NLI of each COI """
            not_i   = arange(channels)!=i
            a_i     = a[i,j]                                                                             # \alpha of COI in fiber span j
            a_k     = a[not_i,j]                                                                         # \alpha of INT in fiber span j
            a_bar_i = a_bar[i,j]                                                                         # \bar{\alpha} of COI in fiber span j
            a_bar_k = a_bar[not_i,j]                                                                     # \bar{\alpha} of INT in fiber span j
            f_i     = fi[i,j]                                                                            # f_i of COI in fiber span j
            f_k     = fi[not_i,j]                                                                        # f_k of INT in fiber span j
            B_i     = Bch[i,j]                                                                           # B_i of COI in fiber span j
            B_k     = Bch[not_i,j]                                                                       # B_k of INT in fiber span j
            Cr_i    = Cr[i,j]                                                                            # Cr  of COI in fiber span j
            Cr_k    = Cr[not_i,j]                                                                        # Cr  of INT in fiber span j
            P_i     = P_ij[i,j]                                                                          # P_i of COI in fiber span j
            P_k     = P_ij[not_i,j]                                                                      # P_k of INT in fiber span j
            
            phi_i  = 3/2*pi**2 *( beta2[j] + pi*beta3[j]*(f_i + f_i) )                                   # \phi_i  of COI in fiber span j
            phi_ik = 2*pi**2 *( f_k - f_i )*( beta2[j] + pi*beta3[j]*(f_i + f_k) )                       # \phi_ik of COI-INT pair in fiber span j

            T_i    = (a_i + a_bar_i - f_i*Ptot[j]*Cr_i)**2                                               # T_i of COI in fiber span j
            T_k    = (a_k + a_bar_k - f_k*Ptot[j]*Cr_k)**2                                               # T_k of INT in fiber span j        
            
            eta_SPM[i,j] = SPM(phi_i, T_i, B_i, a_i, a_bar_i, gamma[j]) *n**eps(B_i, f_i, mean_att_i[i]) # computation of SPM contribution in fiber span j
            eta_XPM[i,j] = XPM(P_i, P_k, phi_ik, T_k, B_i, B_k, a_k, a_bar_k, gamma[j])                  # computation of XPM contribution in fiber span j

    eta_n = sum( ( P_ij/P_ij[:,0,newaxis] )**2 * ( eta_SPM + eta_XPM ) , axis=1)                         # computation of NLI normalized to transmitter power, see Ref. [1, Eq. (5)]
    NLI   = P_ij[:,0]**3 *eta_n                                                                          # Ref. [1, Eq. (1)]
    
    return NLI, eta_n

def removenan(x):
    x[~isfinite(x)] = 0
    return x

def todB(x):
    return 10*log(x)/log(10)