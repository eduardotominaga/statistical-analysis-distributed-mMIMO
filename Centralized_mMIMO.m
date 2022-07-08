%% Reset
clearvars;
close all;
clc;

%% Parameters
d=100;                      % Length of one side of the square area [m]
dmin=5;                     % Minimum distance between the BS and the MTD
h_BS=6;                     % Heigth of the BS [m]
h_MTD=1.5;                  % Height of the MTD [m]
M=64;                       % Number of antenna elements at the BS
fc=3.5;                     % Carrier frequency [GHz]
N0=10^-14.4;                % PSD of the thermal noise [W/Hz]
Pt_dBm=23;                  % Transmit power of the MTD [dBm]
Pt=10^((Pt_dBm-30)/10);     % Transmit power of the MTD [W]
Pt_dB=Pt_dBm-30;            % Transmit power [dB]
B=10e6;                     % Bandwidth [Hz]
NF_dB=7;                    % Noise Figure [dB]
NF=10^(NF_dB/10);           % Noise Figure [linear scale]
MC=1e4;                     % Number of Monte Carlo runs
c=3e8;                      % Speed of the light [m/s]
Position_MTD="Typical";     % Typical case
% Position="Worst";           % Worst case
rng(1)                      % Set the random number generator

%% Define the color and marker maps for the plots:
colorMarkerMap;

%% Determining the distances of the BS and of the MTD

% Determining the position of the BS:
dx_BS=d/2;      % x-position [m]
dy_BS=d/2;      % x-position [m]

% Determining the position of the MTD:
if Position_MTD=="Typical"
    dx_MTD=55;              % x-position of the MTD [m]
    dy_MTD=75;              % y-position of the MTD [m]
elseif Position_MTD=="Worst"
    dx_MTD=d;               % x-position of the MTD [m]
    dy_MTD=d;               % y-position of the MTD [m]
else
    error("Error! Please check the parameters of the script.")
end

% 3D distance between the BS and the MTD:
d_3D=sqrt((dx_BS-dx_MTD).^2+(dy_BS-dy_MTD).^2+(h_BS-h_MTD).^2);

%% Computation of the channel gains

% Small-scale fading samples:
H=(1/sqrt(2))*(randn(M,MC)+1i*randn(M,MC));

% Path-loss model:
n_NLOS=3.19;                                        % Path-loss exponent
sigma_NLOS_dB=7.56;                                 % Variance of the shadowing [dB]
X_dB=sigma_NLOS_dB*randn(1,MC);                     % Shadowing samples [dB]
X=10.^(X_dB/10);                                    % Shadowing samples (linear scale)
A_dB=32.45+20*log10(fc)+10*n_NLOS*log10(d_3D);      % Attenuation due to the distance [dB]
A=10^(A_dB/10);                                     % Attenuation due to the distance (linear scale)
PL_dB=A_dB+X_dB;                                    % Path-loss [dB]
PL=10.^(PL_dB/10);                                  % Path-loss [linear scale]
beta=1./PL;                                         % Large-scale fading

X_bar=10^((1/400)*2*(sigma_NLOS_dB^2)*log(10));     % Average shadowing (linear scale)

G=sqrt(beta).*H;                            % Channel vectors
channelGains=vecnorm(G,2,1).^2;             % Channel gains (linear scale)
channelGains_dB=10*log10(channelGains);     % Channel gains [dB]

% Average channel gain [dB]:
avgChannelGain_simulation=10*log10(mean(channelGains));     % Simulation
avgChannelGain_theory=10*log10((M/A)*X_bar);                % Theory

% Standard deviation of the channel gains (simulation):
stdDevChannelGains_simulation=sqrt(var(channelGains));
stdDevChannelGains_dB_simulation=10*log10(stdDevChannelGains_simulation);

% Coefficient of Variation:
CV=sqrt(var(channelGains))/mean(channelGains);

% rho=Pt/(N0*B*NF);                   % Normalized Transmit SNR
% rho_dB=10*log10(rho);               % Normalized Transmit SNR [dB]
% SNR=rho*ChannelGains;         % Received SNR (linear scale)
% SNR_dB=10*log10(SNR);               % Received SNR [dB]
% disp(mean(SNR_dB))                    % Average SNR

%% Write the results on the console:
disp("Centralized mMIMO")
disp(['M = ',num2str(M)])
disp(['Case = ',num2str(Position_MTD)])
disp(['Average Channel Gain (theory) = ',num2str(avgChannelGain_theory),' dB'])
disp(['Average Channel Gain (simulation) = ',num2str(avgChannelGain_simulation),' dB'])
disp(['Standard Deviation (simulation) = ',num2str(stdDevChannelGains_dB_simulation),' dB'])
disp(['Coefficient of Variation (CV) = ',num2str(CV)])

%% Plotting the results

% Plot the PDF of the channel gains:
fig1=figure(1);
    fig1.Position=[200 300 550 400];
    [curve1,xi]=ksdensity(channelGains_dB);
    plot(xi,curve1,'LineWidth',1.0,'Color',black);
    grid on
    hold on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Channel gain [dB]','Interpreter','latex','FontSize',12)
    ylabel('PDF','Interpreter','latex','FontSize',12)

% Plot the complementary CDF of the channel gains:
fig2=figure(2);
    fig2.Position=[800 300 550 400];
    curve1=ccdfplot(channelGains_dB,1000);
    set(curve1,'LineWidth',1.0,'Color',black)
    grid on
    hold on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Channel gain [dB]','Interpreter','latex','FontSize',12)
    ylabel('CCDF','Interpreter','latex','FontSize',12)

