%% Reset
clearvars;
% close all;
clc;

%% Fixed parameters
d=100;                      % Length of one side of the square area [m]
h_AP=6;                     % Heigth of the AP [m]
h_MTD=1.5;                  % Height of the MTD [m]
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
rng(1)                      % Set the random number generator

%% Define the color and marker maps for the plots:
colorMarkerMap;

%% Configurable parameters
Q=64;                           % Number of APs
S=1;                            % Number of antennas at each AP
M=Q*S;                          % Total number of antennas
Position_MTD="Typical";         % Typical position for the MTD
% Position_MTD="Worst";         % Worst case position for the MTD
% Deployment_APs="Grid";          % Define the distribution of APs
Deployment_APs="Radio Stripe";  % Define the distribution of APs
color=blue;                     % Color of the plot

%% Determining the distances of the BS and of the MTD (worst case)

% Determining the position of the MTD:
if Position_MTD=="Typical"
    dx_MTD=55;                  % x-position of the MTD [m]
    dy_MTD=75;                  % y-position of the MTD [m]
elseif Position_MTD=="Worst"
    if Deployment_APs=="Grid"
        dx_MTD=d;               % x-position of the MTD [m]
        dy_MTD=d;               % y-position of the MTD [m]
    elseif Deployment_APs=="Radio Stripe"
        dx_MTD=d/2;             % x-position of the MTD [m]
        dy_MTD=d/2;             % y-position of the MTD [m]
    else
        error("Error! Please check the configurable parameters.")
    end        
else
    error("Error! Please check the configurable parameters.")
end

% Determining the positions of the APs:
if Deployment_APs=="Grid"
    range_AP=d/(2*sqrt(Q)):d/sqrt(Q):d-d/(2*sqrt(Q));     
    dx_AP=zeros(1,Q);                           % x-position of the APs [m]
    dy_AP=zeros(1,Q);                           % y-position of the APs [m]
    for idx=1:sqrt(Q)
       dx_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP(idx); 
    end
    for idx=1:sqrt(Q)
       dy_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP; 
    end
elseif Deployment_APs=="Radio Stripe"
    range_AP=d/(2*(Q/4)):d/(Q/4):d-d/(2*(Q/4));   
    dx_AP=zeros(1,Q);                   % x-position of the APs [m]
    dy_AP=zeros(1,Q);                   % y-position of the APs [m]
    for idx=1:4         % A square area has four sides.
        switch idx
            case 1
                dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=range_AP;
            case 2
                dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=range_AP;
                dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=d*ones(1,Q/4);
            case 3
                dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=d*ones(1,Q/4);
                dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=flip(range_AP);
            case 4
                dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=flip(range_AP);
        end
    end
else
    error("Error! Please check the configurable parameters.")
end

% 3D distances between the APs and the MTD:
d_3D=zeros(1,Q);
for idx=1:Q
    d_3D(idx)=sqrt((dx_AP(idx)-dx_MTD(1))^2+(dy_AP(idx)-dy_MTD(1))^2+(h_AP-h_MTD)^2);
end

%% Computation of the channel gains

% Small-scale fading samples:
H=(1/sqrt(2))*(randn(M,MC)+1i*randn(M,MC));

% Path-loss model:
n_NLOS=3.19;                        % Path-loss exponent
sigma_NLOS_dB=7.56;                 % Variance of the shadowing [dB]
X_dB=sigma_NLOS_dB*randn(Q,MC);     % Shadowing samples

% Attenuation due to the distance [dB]:
A_dB=zeros(Q,1);
for idx=1:Q
    A_dB(idx)=32.45+20*log10(fc)+10*n_NLOS*log10(d_3D(idx));
end

A=10.^(A_dB/10);        % Attenuation due to distance (linear scale)

PL_dB=A_dB+X_dB;        % Path-loss [dB]
PL=10.^(PL_dB/10);      % Path-loss [linear scale]
beta=1./PL;             % Large-scale fading

% Generate the channel vectors:
G=zeros(M,MC);
for idx1=1:MC
    for idx2=1:Q
        G(1+S*(idx2-1):S*idx2,idx1)=sqrt(beta(idx2,idx1))*H(1+S*(idx2-1):S*idx2,idx1);
    end
end

% Computing the channel gains for all the Monte Carlo runs:
channelGains=vecnorm(G,2,1).^2;             % Linear Scale
channelGains_dB=10*log10(channelGains);     % [dB]

% Average shadowing (linear scale):
X_bar=10^((1/400)*2*(sigma_NLOS_dB^2)*log(10));

% Average channel gain [dB]:
avgChannelGain_simulation=10*log10(mean(channelGains));     % Simulation
avgChannelGain_theory=10*log10(S*sum(1./A)*X_bar);          % Theory
    
% Standard deviation of the channel gains (simulation):
stdDevChannelGains_simulation=sqrt(var(channelGains));                      % Linear Scale
stdDevChannelGains_dB_simulation=10*log10(stdDevChannelGains_simulation);   % [dB]

% Coefficient of Variation:
CV=sqrt(var(channelGains))/mean(channelGains);

%% Write the results on the console:
disp("Distributed mMIMO")
disp(['Configuration = ',num2str(Deployment_APs)])
disp(['Q = ',num2str(Q)])
disp(['S = ',num2str(S)])
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
    plot(xi,curve1,'LineWidth',1.0,'Color',color);
    grid on
    hold on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Channel gain [dB]','Interpreter','latex','FontSize',12)
    ylabel('PDF','Interpreter','latex','FontSize',12)

% Plot the complementary CDF of the channel gains:
fig2=figure(2);
    fig2.Position=[800 300 550 400];
    curve1=ccdfplot(channelGains_dB,1000);
    set(curve1,'LineWidth',1.0,'Color',color)
    grid on    
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Channel gain [dB]','Interpreter','latex','FontSize',12)
    ylabel('CCDF','Interpreter','latex','FontSize',12)

