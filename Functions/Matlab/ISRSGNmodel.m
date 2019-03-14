function [ NLI, eta_n] = ISRSGNmodel( P )
% Returns nonlinear interference power and coefficient for each WDM
% channel. This function implements the ISRS GN model in closed-form
% published in:
%
% D. Semrau, R. I. Killey, P. Bayvel, "A Closed-Form Approximation of the
% Gaussian Noise Model in the Presence of Inter-Channel Stimulated Raman
% Scattering, " J. Lighw. Technol., Early Access, Jan. 2019
% 
%
% INPUTS:
% P - parameter structure in SI units with attributes:
% 
% P.Att       - P.Att(i,j)     is the attenuation coefficient                                                   [Np/m]     of channel i of span j, format: N_ch x n matrix
% P.Att_bar   - P.Att_bar(i,j) is the attenuation coefficient (bar)                                             [Np/m]     of channel i of span j, format: N_ch x n matrix
% P.Cr        - P.Cr(i,j)      is the slope of the linear regression of the normalized Raman gain spectrum      [1/W/m/Hz] of channel i of span j, format: N_ch x n matrix
% 
% P.Pch       - P.Pch(i,j)     is the launch power                                                              [W]        of channel i of span j, format: N_ch x n matrix
% P.fi        - P.fi(i,j)      is center frequency relative to the reference frequency (3e8/P.RefLambda)        [Hz]       of channel i of span j, format: N_ch x n matrix
% P.Bch       - P.Bch(i,j)     is the bandwidth                                                                 [Hz]       of channel i of span j, format: N_ch x n matrix
%
% P.Length    - P.Length(1,j)  is the span length                                                               [m]        of span j             , format: 1 x n matrix   
% P.D         - P.D(1,j)       is the dispersion coefficient                                                    [s/m^2]    of span j             , format: 1 x n matrix
% P.S         - P.S(1,j)       is the span length                                                               [s/m^3]    of span j             , format: 1 x n matrix
% P.Gamma     - P.gamma(1,j)  is the span length                                                                [1/W/m]    of span j             , format: 1 x n matrix
% P.RefLambda - P.RefLambda    is the reference wavelength (where beta2, beta3 are defined)                     [m]                              , format: 1 x 1 matrix
%
% P.coherent  - P.coherent is 1 or 0 for coherent or incoherent NLI accumulation across multiple fiber spans 
%
% RETURNS:
% NLI   - Nonlinear Interference Power      [W]     , format: Nch x 1 matrix 
% eta_n - Nonlinear Interference coeffcient [1/W^2] , format: Nch x 1 matrix 
%
% Author: Daniel Semrau, Eric Sillekens, Jan 2019.


%% Parameter Definition
                                                                                               
a      = P.Att;                                                                                                    % attenuation coefficent       of channel i in fiber span j, format: N_ch x n matrix
a_bar  = P.Att_bar;                                                                                                % attenuation coefficent (bar) of channel i in fiber span j, format: N_ch x n matrix
Cr     = P.Cr;                                                                                                     % Raman gain spectrum slope    of channel i in fiber span j, format: N_ch x n matrix

P_ij   = P.Pch;                                                                                                    % Launch power                 of channel i in fiber span j, format: N_ch x n matrix 

fi     = P.fi;                                                                                                     % Center frequency             of channel i in fiber span j, format: N_ch x n matrix 
Bch    = P.Bch;                                                                                                    % Channel bandwidth            of channel i in fiber span j, format: N_ch x n matrix

n      = P.n;                                                                                                      % number of spans                                          , format: 1 x 1 matrix
lambda = P.RefLambda;                                                                                              % Reference wavelength                                     , format: 1 x 1 matrix
beta2  = -P.D*lambda^2/(2*pi*3e8);                                                                                 % beat2                                     in fiber span j, format: 1 x n matrix
beta3  = lambda.^2/(2*pi*3e8).^2*(lambda.^2*P.S+2*lambda*P.D);                                                     % beat3                                     in fiber span j, format: 1 x n matrix
gamma  = P.Gamma;                                                                                                  % Nonlinaerity coefficient                  in fiber span j, format: 1 x n matrix 
L      = P.Length;                                                                                                 % Fiber length                              in fiber span j, format: 1 x n matrix
Ptot   = sum(P_ij,1);                                                                                              % Total power launched                    into fiber span j, format: 1 x n matrix 

%% NLI computation of each span and normalization to transmitter power level 

for j=1:n                                                                                                          % compute the NLI for each span, see [1, Eq. (5)]

    % SPM and XPM Closed-form Formula Definition

    SPM = @(phi_i, T_i, B_i, a, a_bar, gamma)            ...
          4/9*gamma^2/B_i^2*pi/(phi_i*a_bar*(2*a+a_bar)) ...  
          *( (T_i-a^2)/a*asinh(phi_i*B_i^2/a/pi) + ((a+a_bar)^2-T_i)/(a+a_bar)*asinh(phi_i*B_i^2./(a+a_bar)/pi) ); % closed-form formula for SPM contribution, see Ref. [1, Eq. (9-10)]

    XPM = @(Pi, Pk, phi_ik, T_k, B_i, B_k, a, a_bar, gamma)                                                       ...
          32/27*sum(removenan(  (Pk./Pi).^2.*gamma^2 ./ ( B_k.*phi_ik.*a_bar.*(2.*a+a_bar) )                      ...
          .*( (T_k-a.^2)./a.*atan(phi_ik.*B_i./a) + ((a+a_bar).^2-T_k)./(a+a_bar).*atan(phi_ik.*B_i./(a+a_bar)) ) ...
                   )         );                                                                                    % closed-form formula for XPM contribution, see Ref. [1, Eq. (11)] 
    % Average Coherence Factor

    mean_att_i = mean(a,2);                                                                                    % average attenuation coefficent for channel i
    mean_L     = mean(L,2);                                                                                    % average fiber length 
    
    if P.coherent == 1
        
        eps = @(B_i, f_i, a_i) ...
              (3/10)*log(1+(6/a_i)/(mean_L*asinh(pi^2/2*abs( mean(beta2) + 2*pi*mean(beta3)*f_i )/a_i*B_i^2)));   % closed-for formula for average coherence factor extended by dispersion slope, cf. Ref. [2, Eq. (22)]
            
    else
        
        eps = @(B_i, f_i, a_i) 0;
        
    end
    % Calculation of nonlinear interference (NLI) power in fiber span j

    for i=1:length(fi)                                                                                             % compute the NLI of each COI

        % define variables for [1, Eq. (10-11)]
        
        a_i     = a(i,j);                                                                                          % \alpha of COI in fiber span j
        a_k     = a(:,j);                                                                                          % \alpha of INT in fiber span j
        a_bar_i = a_bar(i,j);                                                                                      % \bar{\alpha} of COI in fiber span j
        a_bar_k = a_bar(:,j);                                                                                      % \bar{\alpha} of INT in fiber span j
        f_i     = fi(i,j);                                                                                         % f_i of COI in fiber span j
        f_k     = fi(:,j);                                                                                         % f_k of INT in fiber span j
        B_i     = Bch(i,j);                                                                                        % B_i of COI in fiber span j
        B_k     = Bch(:,j);                                                                                        % B_k of INT in fiber span j
        Cr_i    = Cr(i,j);                                                                                         % Cr  of COI in fiber span j
        Cr_k    = Cr(:,j);                                                                                         % Cr  of INT in fiber span j
        P_i     = P_ij(i,j);                                                                                       % P_i of COI in fiber span j
        P_k     = P_ij(:,j);                                                                                       % P_k of INT in fiber span j


        phi_i  = 3/2*pi^2              .*( beta2(j) + pi*beta3(j)*(f_i + f_i) );                                   % \phi_i  of COI in fiber span j
        phi_ik = 2  *pi^2*( f_k - f_i ).*( beta2(j) + pi*beta3(j)*(f_i + f_k) );                                   % \phi_ik of COI-INT pair in fiber span j

        T_i    = (a_i + a_bar_i - f_i.*Ptot(j).*Cr_i).^2;                                                          % T_i of COI in fiber span j
        T_k    = (a_k + a_bar_k - f_k.*Ptot(j).*Cr_k).^2;                                                          % T_k of INT in fiber span j

        eta_SPM(i,j) = SPM(phi_i, T_i, B_i, a_i, a_bar_i, gamma(j))                *n^eps(B_i, f_i, mean_att_i(i)); % computation of SPM contribution in fiber span j
        eta_XPM(i,j) = XPM(P_i, P_k, phi_ik, T_k, B_i, B_k, a_k, a_bar_k, gamma(j));                               % computation of XPM contribution in fiber span j

    end
    
end

%% Normalizatio of NLI power generated in each span to transmitter power 

eta_n = sum( ( P_ij./P_ij(:,1) ).^2 .* ( eta_SPM + eta_XPM ) ,2);                                                   % computation of NLI normalized to transmitter power, see Ref. [1, Eq. (5)]
NLI   = P_ij(:,1).^3.*eta_n;                                                                                        % Ref. [1, Eq. (1)]

       
end
       
function x = removenan(x)                                                                                           % removes invalid k=i term in the XPM contribution

    x(isnan(x)) = 0;

end

%% References
%
% [1] D. Semrau, R. I. Killey, P. Bayvel, "A Closed-Form Approximation of the Gaussian Noise Model in the Presence of Inter-Channel Stimulated Raman Scattering, "  J. Lighw. Technol., vol. x, no. x, pp.xxxx-xxxx, Jan. 2019
% [2] P. Poggiolini, "The GN model of non-linear propagation in uncompensated coherent optical systems, " J. Lighw. Technol., vol. 30, no. 24, pp.3857-3879, Dec. 2012
%


