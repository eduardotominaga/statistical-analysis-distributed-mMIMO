%% Reset
clearvars;
close all;
clc;

%% Parameters
d=100;                      % Length of one side of the square area [m]
h_AP=6;                     % Heigth of the AP [m]
h_MTD=1.5;                  % Height of the MTD [m]
Q=64;                       % Number of APs
S=1;                        % Number of antennas at each AP
M=Q*S;                      % Total number of antennas
fc=3.5;                     % Carrier frequency [GHz]
MC=1e4;                     % Number of Monte Carlo runs
Position="Typical";         % Typical case
% Position="Worst";           % Worst case
N_Q_AP=[1 4 8 16];          % Number of APs in the subset Q_AP
rng(1)                      % Set the random number generator

%% Define the color and marker maps for the plots:
colorMarkerMap;

%% Writing the parameters on the console:
disp("Distributed mMIMO - Grid Configuration")
disp("Parameters:")
disp(['Q = ',num2str(Q)])
disp(['S = ',num2str(S)])

%% Determining the positions of the APs and of the MTD

% Determining the position of the MTD:
if Position=="Typical"
    dx_MTD=55;              % x-position of the MTD [m]
    dy_MTD=75;              % y-position of the MTD [m]
elseif Position=="Worst"
    dx_MTD=d;               % x-position of the MTD [m]
    dy_MTD=d;               % y-position of the MTD [m]
else
    error("Error! Please check the parameters of the script.")
end

% Determining the position of the APs:
range_AP=d/(2*sqrt(Q)):d/sqrt(Q):d-d/(2*sqrt(Q));       % Auxiliary range   
dx_AP=zeros(1,Q);                                       % x-position of the APs [m]
dy_AP=zeros(1,Q);                                       % y-position of the APs [m]
% x-positions of the APs:
for idx=1:sqrt(Q)
   dx_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP(idx); 
end
% y-positions of the APs:
for idx=1:sqrt(Q)
   dy_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP; 
end

% 3D distances between the APs and the MTD:
d_3D=zeros(1,Q);
for idx=1:Q
    d_3D(idx)=sqrt((dx_AP(idx)-dx_MTD)^2+(dy_AP(idx)-dy_MTD)^2+(h_AP-h_MTD)^2);
end

%% Computation of the channel gains

% Small-scale fading samples:
H=(1/sqrt(2))*(randn(M,MC)+1i*randn(M,MC));     

% Path-loss model:
n_NLOS=3.19;                                % Path-loss exponent
sigma_NLOS_dB=7.56;                         % Variance of the shadowing [dB]
X_dB=sigma_NLOS_dB*randn(Q,MC);             % Shadowing samples

% Attenuation due to the distance for all the APs:
A_dB=zeros(Q,1);
for idx=1:Q
    A_dB(idx)=32.45+20*log10(fc)+10*n_NLOS*log10(d_3D(idx));
end

A=10.^(A_dB/10);        % Attenuation due to distance (linear scale)

PL_dB=A_dB+X_dB;        % Path-loss [dB]
PL=10.^(PL_dB/10);      % Path-loss [linear scale]
beta=1./PL;             % Large-scale fading coefficients

% Generate the wireless channel vectors:
G=zeros(M,MC);
for idx1=1:MC
    for idx2=1:Q
        G(1+S*(idx2-1):S*idx2,idx1)=sqrt(beta(idx2,idx1))*H(1+S*(idx2-1):S*idx2,idx1);
    end
end

X_bar=10^((1/400)*2*(sigma_NLOS_dB^2)*log(10));     % Average shadowing (linear scale)

% Computing the channel gains for all the Monte Carlo runs:
channelGains=vecnorm(G,2,1).^2;             % Linear scale
channelGains_dB=10*log10(channelGains);     % [dB]

% Average channel gain [dB]:
avgChannelGain_simulation=10*log10(mean(channelGains));     % Simulation
avgChannelGain_theory=10*log10(S*sum(1./A)*X_bar);                 % Theory

% Standard deviation  of the channel gains (simulation):
stdDevChannelGains_simulation=sqrt(var(channelGains));
stdDevChannelGains_dB_simulation=10*log10(stdDevChannelGains_simulation);

CV=sqrt(var(channelGains))/mean(channelGains);        % Coefficient of Variation

%% Subset Q_AP os APs closest to the MTD:

[~,idx_ord_AP]=sort(d_3D);      % Indices of the APs in increasing order of distances
for x=1:length(N_Q_AP)    
    Q_AP=idx_ord_AP(1:N_Q_AP(x));   % Indices of the APs belonging to Q_AP

    % Select the channel gains associated to the APs in Q_AP:
    G_Q_AP=zeros(N_Q_AP(x)*S,MC);
    for idx=1:N_Q_AP(x)
        G_Q_AP(1+S*(idx-1):S*idx,:)=G(1+S*(Q_AP(idx)-1):S*Q_AP(idx),:);
    end

    % Computing the channel gains for all the Monte Carlo runs:
    channelGains_Q_AP=vecnorm(G_Q_AP,2,1).^2;               % Linear scale
    channelGains_Q_AP_dB=10*log10(channelGains_Q_AP);       %[dB]
    
    % Average channel gain [dB]:
    avgChannelGain_Q_AP_simulation=10*log10(mean(channelGains_Q_AP));   % Simulation
    avgChannelGain_Q_AP_theory=10*log10(S*sum(1./A(Q_AP))*X_bar);       % Theory 

    % Standard deviation  of the channel gains (simulation):
    stdDevNormChannelVectors_Q_AP_simulation=sqrt(var(channelGains_Q_AP));
    stdDevNormChannelVectors_dB_Q_AP_simulation=10*log10(stdDevNormChannelVectors_Q_AP_simulation);

    % Coefficient of Variation:
    CV_Q_AP=sqrt(var(channelGains_Q_AP))/mean(channelGains_Q_AP);
        
    % Ratio of the channel gains:
    ratioChannelGains=mean(channelGains_Q_AP./channelGains);
    
    % Writing the results on the console:
    fprintf("\n")
    disp(['|Q_AP| = ',num2str(N_Q_AP(x))])
    disp(['Average Channel Gain (theory) = ',num2str(avgChannelGain_Q_AP_theory),' dB'])
    disp(['Average Channel Gain (simulation) = ',num2str(avgChannelGain_Q_AP_simulation),' dB'])
    disp(['Standard Deviation (simulation) = ',num2str(stdDevChannelGains_dB_simulation),' dB'])
    disp(['Coefficient of Variation (CV) = ',num2str(CV)])
    disp(['Ratio of the channel gains = ',num2str(ratioChannelGains)])

    % Plot the PDF of the channel gains:
    figure(1);
        [curve1,xi]=ksdensity(channelGains_Q_AP_dB);
        plot(xi,curve1,'LineWidth',1.5,'color',myColorMap(x,:),'LineWidth',1.5,...
                 'DisplayName',strcat('$|Q_{AP}|=',num2str(N_Q_AP(x)),'$'))
        hold on
    
    % Plot the complementary CDF of the channel gains:
    figure(2);
        curve1=ccdfplot(channelGains_Q_AP_dB,1000);
        set(curve1,'LineWidth',1.5,'color',myColorMap(x,:),'LineWidth',1.5,...
                 'DisplayName',strcat('$|Q_{AP}|=',num2str(N_Q_AP(x)),'$'))
        hold on
end

%% Writing the results on the console ("All APs" case):
fprintf("\n")
disp("All APs:")
disp(['Average Channel Gain (theory) = ',num2str(avgChannelGain_theory),' dB'])
disp(['Average Channel Gain (simulation) = ',num2str(avgChannelGain_simulation),' dB'])
disp(['Standard Deviation (simulation) = ',num2str(stdDevChannelGains_dB_simulation),' dB'])
disp(['Coefficient of Variation (CV) = ',num2str(CV)])

%% Plotting the results

% Plot the PDF for the case of all APs:
fig1=figure(1);
    fig1.Position=[200 300 550 400];
    [curve1,xi]=ksdensity(channelGains_dB);
    plot(xi,curve1,'LineWidth',1.5,'color',yellow,'LineWidth',1.5,...
             'DisplayName',strcat('All APs'))
    grid on
    xlim([-110 -50])
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    xlabel('Channel gain [dB]','Interpreter','latex','FontSize',14)
    ylabel('PDF','Interpreter','latex','FontSize',14)
    leg=legend('show');
    set(leg,'Interpreter','latex','FontSize',12,'Location','Northwest')
    str=strcat('PDF_Subset_APs_',Position,'_Case.png');
    saveas(fig1,str)
    str=strcat('PDF_Subset_APs_',Position,'_Case.eps');
    saveas(fig1,str,'epsc')
    
% Plot the CCDF for the case of all APs:
fig2=figure(2);
    fig2.Position=[800 300 550 400];
    curve1=ccdfplot(channelGains_dB,1000);
    set(curve1,'LineWidth',1.5,'color',yellow,'LineWidth',1.5,...
             'DisplayName',strcat('All APs'))
    grid on
    xlim([-110 -50])
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    xlabel('Channel gain [dB]','Interpreter','latex','FontSize',14)
    ylabel('CCDF','Interpreter','latex','FontSize',14)
    leg=legend('show');
    set(leg,'Interpreter','latex','FontSize',12,'Location','Southwest')
    str=strcat('CCDF_Subset_APs_',Position,'_Case.png');
    saveas(fig2,str)
    str=strcat('CCDF_Subset_APs_',Position,'_Case.eps');
    saveas(fig2,str,'epsc')
