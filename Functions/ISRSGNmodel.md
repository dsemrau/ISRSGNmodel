## The ISRS GN model function (syntax)

The signal-to-noise ratio (SNR) after coherent detection and digital dispersion compensation is given by
SNR = <sup>P</sup>&frasl;<sub>(P<sub>ASE</sub>+P<sub>NLI</sub>)</sub>, where P is the launch power, P<sub>ASE</sub> is the linear noise power and P<sub>NLI</sub> is the nonlinear interference power. All quantities are channel (frequency) dependent. 

The ISRS GN model function returns the nonlinear interference (NLI) power P<sub>NLI</sub> = ηP<sup>3</sup> [W] and the NLI coefficient η [<sup>1</sup>&frasl;<sub>W<sup>2</sup></sub>] for every channel slot (N<sub>ch</sub> × 1). The function needs the following variables as input parameters:

| Input Parameter                                          | Variable | Dimension | Format   |
|----------------------------------------------------------|----------|-----------|----------|
| Attenuation coefficient (α)                              | Att      | <sup>Np</sup>&frasl;<sub>m</sub>                    | N<sub>ch</sub> × n |
| Attenuation coefficient (bar) (ᾱ)                        | Att_bar  | <sup>Np</sup>&frasl;<sub>m</sub>                    | N<sub>ch</sub> × n |
| Linear regression of Raman gain function (C<sub>r</sub>) | Cr       | <sup>1</sup>&frasl;<sub>(W·m·Hz)</sub>| N<sub>ch</sub> × n |         
| Channel launch power (P<sub>i</sub>)                     | Pch      | W                                                   | N<sub>ch</sub> × n |
| Channel center frequency (f<sub>i</sub>)                 | fi       | Hz                                                  | N<sub>ch</sub> × n |
| Channel bandwidth (B<sub>i</sub>)                        | Bch      | Hz                                                  | N<sub>ch</sub> × n |
| Number of spans (n)                                      | n        |                                                     | scalar      |
| Fiber length (L)                                         | L        | m                                                   | 1 × n |
| Dispersion coefficient (D)                               | D        | <sup>s</sup>&frasl;<sub>m<sup>2</sup></sub>         | 1 × n |
| Dispersion slope (S)                                     | S        | <sup>s</sup>&frasl;<sub>m<sup>3</sup></sub>         | 1 × n |
| Nonlinearity coefficient (γ)                             | γ        | <sup>1</sup>&frasl;<sub>(W·m)</sub>                 | 1 × n |
| Reference wavelength (λ)                                 | λ        | m                                                   | scalar |
| coherent NLI accumulation                                | coherent | boolean                                             | scalar |

The format N<sub>ch</sub> × n is represented by a matrix, where N<sub>ch</sub> is the number of channel slots and n is the number of fiber spans. Essentially, this expressed the channel and span dependency on certain parameters. Generally, it holds that ᾱ=α. For more information and explanation, please see [1]. The input parameters are passed to the function as an object which handled slightly differently in Python and in Matlab. 

#### Python

In Python, the input parameters are passed as a dictionary with the input parameters as keys/values. An example is listed below

```Python
n = 1                                                                 # number of spans
Bch = 40.004e9                                                        # WDM channel bandwidth
channels = 251                                                        # number of channels
spacing = 40.005e9                                                    # WDM channel spacing

P = {
       
    'fi'       : np.repeat(np.reshape( 
                 (np.arange(channels) - (channels-1)/2)*spacing
                 ,[-1,1]),n,axis=1),                                  # center frequencies 
    'n'        : n,                                                   # number of spans
    'Bch'      : np.tile(40.004e9,[channels,n]),                      # channel bandwith
    'RefLambda': 1550e-9,                                             # reference wavelength 
    'D'        :  17    *1e-12/1e-9/1e3      * np.ones(n),            # dispersion coefficient   
    'S'        : 0.067 *1e-12/1e-9/1e3/1e-9  * np.ones(n),            # dispersion slope  
    'Att'      : 0.2   /4.343/1e3 * np.ones([channels, n]),           # attenuation coefficient     
    'Cr'       : 0.028 /1e3/1e12  * np.ones([channels, n]),           # Raman gain spectrum slope
    'gamma'    : 1.2 /1e3  * np.ones(n),                              # nonlinearity coefficient
    'Length'   : 100 *1e3  * np.ones(n),                              # fiber length
    'coherent' : 1                                                    # coherent NLI accumulation
    'Pch'      : 10**(0/10)*0.001 * np.ones([channels, n])            # channel launch power
    }

P['Att_bar'] = P['Att']                                               # attenuation (bar)

NLI_power, NLI_coefficient = ISRSGNmodel(**P)                         # ISRS GN model function call
```

#### Matlab

In Matlab, the input parameters are passed as a struct with the input parameters as attributes. An example is listed below

```Matlab
P.n = 1;                                                             % number of spans

Bch      = 40.004e9;                                                 % WDM channel bandwidth
channels = 251;                                                      % number of channels
spacing  = 40.005e9;                                                 % WDM channel spacing

P.fi  = -channels*spacing/2+0.5*spacing+(0:(channels-1))'*spacing;   % center frequencies
P.fi  = repmat(P.fi, 1, P.n);                                        % center frequencies 

P.Bch       = Bch * ones(channels, P.n);                             % channel bandwidth
P.RefLambda = 1550e-9;                                               % reference wavelength 

P.D = 17    *1e-12/1e-9/1e3      * ones(1, P.n);                     % dispersion coefficient   
P.S = 0.067 *1e-12/1e-9/1e3/1e-9 * ones(1, P.n);                     % dispersion slope       

P.Att     = 0.2   /4.343/1e3 * ones(channels, P.n);                  % attenuation coefficient   
P.Att_bar = P.Att;                                                   % attenuation coefficient (bar)
P.Cr      = 0.028 /1e3/1e12  * ones(channels, P.n);                  % Raman gain spectrum slope 

P.Gamma  = 1.2 /1e3 * ones(1, P.n);                                  % nonlinearity coefficient

P.Length = 100 *1e3 * ones(1, P.n);                                  % fiber length

P.coherent = 1;                                                      % coherent NLI accumulation

P.Pch      = 10.^(0/10)*0.001 * ones(channels, P.n);                 % launch power per channel

[NLI_power, NLI_coeffcient] = ISRSGNmodel( P );               % ISRS GN model function call
```

>[1] D. Semrau, R. I. Killey, P. Bayvel, "[A Closed-Form Approximation of the Gaussian Noise Model in the Presence of Inter-Channel Stimulated Raman Scattering]((https://ieeexplore.ieee.org/document/8625492)), " J. Lighw. Technol., Early Access, DOI: 10.1109/JLT.2019.2895237, Jan. 2019. 
