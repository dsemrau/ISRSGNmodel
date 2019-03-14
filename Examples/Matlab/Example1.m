clear
close all

addpath('../../Functions/Matlab')
addpath('../../Data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute and plot results of Fig. 5a) and Fig. 5b) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define parameters

P.n = 1;                                                                   % number of spans

Bch      = 40.004e9;                                                       % WDM channel bandwidth
channels = 251;                                                            % number of channels
spacing  = 40.005e9;                                                       % WDM channel spacing

P.fi  = -channels*spacing/2+0.5*spacing+(0:(channels-1))'*spacing;         % center frequencies of WDM channels (relative to reference frequency)
P.fi  = repmat(P.fi, 1, P.n);                                              % same center frequencies for each span

P.Bch = Bch              * ones(channels, P.n);                            % same channel bandwidth  for each span 

P.RefLambda = 1550e-9;                                                     % reference wavelength 

P.D = 17    *1e-12/1e-9/1e3      * ones(1, P.n);                           % dispersion coefficient (same) for each span  
P.S = 0.067 *1e-12/1e-9/1e3/1e-9 * ones(1, P.n);                           % dispersion slope       (same) for each span

P.Att     = 0.2   /4.343/1e3 * ones(channels, P.n);                        % attenuation coefficient     (same) for each channel and span 
P.Att_bar = P.Att;                                                         % attenuation coefficient bar (same) for each channel and span
P.Cr      = 0.028 /1e3/1e12 * ones(channels, P.n);                         % Raman gain spectrum slope   (same) for each channel and span

P.Gamma = 1.2 /1e3  * ones(1, P.n);                                        % nonlinearity coefficient (same) for each span

P.Length = 100 *1e3 * ones(1, P.n);                                        % fiber length (same) for each span

P.coherent = 1;                                                            % NLI is added coherently across multiple spans

%% Computation of NLI coefficients

% compute results of Fig. 5a)

P.Pch          = 10.^(0/10)*0.001 * ones(channels, P.n);                   % launch power per channel (same) for each channel and span 
eta_1span_0dBm = ISRSGNmodel( P )./P.Pch.^3;                        % compute NLI

P.Pch          = 10.^(2/10)*0.001 * ones(channels, P.n);                   % launch power per channel (same) for each channel and span 
eta_1span_2dBm = ISRSGNmodel( P )./P.Pch.^3;                        % compute NLI

P_noISRS          = P;                                                     % copy parameter struct
P_noISRS.Cr       = zeros(channels, P.n);                                  % set Raman gain spectrum to zero (switch off ISRS)
eta_1span_noISRS  = ISRSGNmodel( P_noISRS )./P_noISRS.Pch.^3;       % compute NLI

% compute results of Fig. 5b)
 
Leff = (1-exp(-P.Att(1)*P.Length))/P.Att(1);                               % effective length

rho_array = linspace(1e-2,13.2,1e2);                                       % x-axis in Fig. 5b)
Pch_array = rho_array/4.3/channels/Leff/P.Cr(1)/channels/P.Bch(1);         % convert to launch powers, see. [1, Eq. (3)]

for i = 1:length(Pch_array)                                                % compute NLI for every launch power (x-axis in Fig. 5b)
   
    P.Pch    = Pch_array(i)*ones(channels,1);                              % launch power per channel (same) for each channel and span 
    NLI(:,i) = ISRSGNmodel( P )./P.Pch.^3;                          % compute NLI
    
end

%% Plot results

figure;
subplot(1,2,1)
plot(P.fi(:,1)*1e-12,10*log10(eta_1span_noISRS),'r--','LineWidth',1.5,'DisplayName','no ISRS')
hold all
plot(P.fi(:,1)*1e-12,10*log10(eta_1span_0dBm),'b--','LineWidth',1.5,'DisplayName','0 dBm/ch.')
plot(P.fi(:,1)*1e-12,10*log10(eta_1span_2dBm),'g--','LineWidth',1.5,'DisplayName','2 dBm/ch.')
grid on
xlabel('Channel frequency f_i [THz]')
ylabel('NLI coefficient \eta_1 [dB(1/W^2)]')
legend('show','Location','best')
title('Fig. 5 a)')
axis([-5.3 5.3 26 32.5])

subplot(1,2,2)
plot(rho_array,10*log10(NLI(13,:)/NLI(13,1)),'--','LineWidth',1.5,'DisplayName','-4500 Ghz')
hold all
plot(rho_array,10*log10(NLI(76,:)/NLI(76,1)),'--','LineWidth',1.5,'DisplayName','-2000 Ghz')
plot(rho_array,10*log10(NLI(126,:)/NLI(126,1)),'--','LineWidth',1.5,'DisplayName','0 Ghz')
plot(rho_array,10*log10(NLI(176,:)/NLI(176,1)),'--','LineWidth',1.5,'DisplayName','2000 Ghz')
plot(rho_array,10*log10(NLI(238,:)/NLI(238,1)),'--','LineWidth',1.5,'DisplayName','4500 Ghz')
grid on
xlabel('ISRS power transfer \Delta \rho(L) [dB]')
ylabel('NLI coefficient \Delta\eta_1 [dB(1/W^2)]')
legend('show','Location','best')
title('Fig. 5 b)')
axis([0 13.2 -4 4])


clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute and plot results of Fig. 6a) and Fig. 6b) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define parameters

P.n = 6;                                                                   % number of spans

Bch      = 40.004e9;                                                       % WDM channel bandwidth
channels = 251;                                                            % number of channels
spacing  = 40.005e9;                                                       % WDM channel spacing

P.fi  = -channels*spacing/2+0.5*spacing+(0:(channels-1))'*spacing;         % center frequencies of WDM channels (relative to reference frequency)
P.fi  = repmat(P.fi, 1, P.n);                                              % same center frequencies for each span

P.Bch = Bch              * ones(channels, P.n);                            % same channel bandwidth  for each span 

P.Pch                = 10.^(0/10)*0.001 * ones(channels, P.n);             % launch power per channel (same) for each channel and span

P.RefLambda = 1550e-9;                                                     % reference wavelength 

P.D = 17    *1e-12/1e-9/1e3      * ones(1, P.n);                           % dispersion coefficient (same) for each span  
P.S = 0.067 *1e-12/1e-9/1e3/1e-9 * ones(1, P.n);                           % dispersion slope       (same) for each span

P.Att     = 0.2   /4.343/1e3 * ones(channels, P.n);                        % attenuation coefficient     (same) for each channel and span 
P.Att_bar = P.Att;                                                         % attenuation coefficient bar (same) for each channel and span
P.Cr      = 0.028 /1e3/1e12 * ones(channels, P.n);                         % Raman gain spectrum slope   (same) for each channel and span

P.Gamma = 1.2 /1e3  * ones(1, P.n);                                        % nonlinearity coefficient (same) for each span

P.Length = 100 *1e3 * ones(1, P.n);                                        % fiber length (same) for each span


%% Computation of NLI coefficients

% compute results with coherent accumultion

P.coherent = 1;                                                            % NLI is added coherently across multiple spans

eta_6span_coh = ISRSGNmodel( P )./P.Pch(:,1).^3;                         % compute NLI

P_noISRS             = P;                                                  % copy parameter struct
P_noISRS.Cr          = zeros(channels, P.n);                               % set Raman gain spectrum to zero (switch off ISRS)
eta_6span_noISRS_coh = ISRSGNmodel( P_noISRS )./P_noISRS.Pch(:,1).^3;% compute NLI

% compute results with coherent accumultion

P.coherent = 0;                                                            % NLI is added incoherently across multiple spans

eta_6span_inc        = ISRSGNmodel( P )./P.Pch(:,1).^3;             % compute NLI

P_noISRS             = P;                                                  % copy parameter struct
P_noISRS.Cr          = zeros(channels, P.n);                               % set Raman gain spectrum to zero (switch off ISRS)
eta_6span_noISRS_inc = ISRSGNmodel( P_noISRS )./P_noISRS.Pch(:,1).^3;% compute NLI

%% Plot results

figure;
subplot(1,2,1)
plot(P.fi(:,1)*1e-12,10*log10(eta_6span_noISRS_coh),'b--','LineWidth',1.5,'DisplayName','coherent, \epsilon \neq 0')
hold all
plot(P.fi(:,1)*1e-12,10*log10(eta_6span_noISRS_inc),'b:','LineWidth',1.5,'DisplayName','incoherent, \epsilon = 0')
grid on
xlabel('Channel frequency f_i [THz]')
ylabel('NLI coefficient \eta_6 [dB(1/W^2)]')
legend('show','Location','best')
title('Fig. 6 a), without ISRS')
axis([-5.3 5.3 34.5 39.5])

subplot(1,2,2)
plot(P.fi(:,1)*1e-12,10*log10(eta_6span_coh),'b--','LineWidth',1.5,'DisplayName','coherent, \epsilon \neq 0')
hold all
plot(P.fi(:,1)*1e-12,10*log10(eta_6span_inc),'b:','LineWidth',1.5,'DisplayName','incoherent, \epsilon = 0')
grid on
xlabel('Channel frequency f_i [THz]')
ylabel('NLI coefficient \eta_6 [dB(1/W^2)]')
legend('show','Location','best')
title('Fig. 6 b), with ISRS')
axis([-5.3 5.3 34.5 39.5])


clear 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute and plot results of Fig. 8a) and Fig. 8b) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define parameters

P.n = 6;                                                                   % number of spans

Bch      = 40.004e9;                                                       % WDM channel bandwidth
channels = 251;                                                            % number of channels
spacing  = 40.005e9;                                                       % WDM channel spacing

P.fi  = -channels*spacing/2+0.5*spacing+(0:(channels-1))'*spacing;         % center frequencies of WDM channels (relative to reference frequency)
P.fi  = repmat(P.fi, 1, P.n);                                              % same center frequencies for each span

P.Bch = Bch              * ones(channels, P.n);                            % same channel bandwidth  for each span 

P.RefLambda = 1550e-9;                                                     % reference wavelength 

P.D = 17    *1e-12/1e-9/1e3      * ones(1, P.n);                           % dispersion coefficient (same) for each span  
P.S = 0.067 *1e-12/1e-9/1e3/1e-9 * ones(1, P.n);                           % dispersion slope       (same) for each span

P.Att     = 0.2   /4.343/1e3 * ones(channels, P.n);                        % attenuation coefficient     (same) for each channel and span 
P.Att_bar = P.Att;                                                         % attenuation coefficient bar (same) for each channel and span
P.Cr      = 0.028 /1e3/1e12 * ones(channels, P.n);                         % Raman gain spectrum slope   (same) for each channel and span

P.Gamma = 1.2 /1e3  * ones(1, P.n);                                        % nonlinearity coefficient (same) for each span

P.Length = 100 *1e3 * ones(1, P.n);                                        % fiber length (same) for each span

P.coherent = 1; 

%% Computation of NLI coefficients

COI_idx = 1:5:channels;

% compute 80% network utilization 

P.Pch        = getfield(load('data_Fig8.mat'),'Pch_80');                   % launch power per channel for each channel and span
eta_6span_80 = ISRSGNmodel( P )./P.Pch(:,1).^3;                     % compute NLI
eta_6span_80 = eta_6span_80(COI_idx);                                      % disable non-COI channels (add/drop channels)

% compute 90% network utilization 

P.Pch        = getfield(load('data_Fig8.mat'),'Pch_90');                   % launch power per channel for each channel and span
eta_6span_90 = ISRSGNmodel( P )./P.Pch(:,1).^3;                     % compute NLI
eta_6span_90 = eta_6span_90(COI_idx);                                      % disable non-COI channels (add/drop channels)

%% Plot results

figure;
subplot(1,2,1)
plot(P.fi(COI_idx,1)*1e-12,10*log10(eta_6span_80),'b--','LineWidth',1.5,'DisplayName','coherent, \epsilon \neq 0')
grid on
xlabel('Channel frequency f_i [THz]')
ylabel('NLI coefficient \eta_6 [dB(1/W^2)]')
legend('show','Location','best')
title('Fig. 8 a), 80% network utilization')
axis([-5.3 5.3 33.5 39])

subplot(1,2,2)
plot(P.fi(COI_idx,1)*1e-12,10*log10(eta_6span_90),'b--','LineWidth',1.5,'DisplayName','coherent, \epsilon \neq 0')
grid on
xlabel('Channel frequency f_i [THz]')
ylabel('NLI coefficient \eta_6 [dB(1/W^2)]')
legend('show','Location','best')
title('Fig. 8 b), 90% network utilization')
axis([-5.3 5.3 33.5 39])


clear 

%% References
%
% [1] D. Semrau, R. I. Killey, P. Bayvel, "A Closed-Form Approximation of the Gaussian Noise Model in the Presence of Inter-Channel Stimulated Raman Scattering, "  J. Lighw. Technol., Early Access, Jan. 2019
%






