%% Reset
clearvars;
close all;
clc;

%% Parameters
d=100;              % Length of one side of the square area [m]
dmin=5;             % Minimum distance between the BS and the MTD
K=20;               % Number of MTDs
rng(3)              % Set the random number generator

%% Colors for the plots
blue = [57 106 177]./255;

%% Define the color and marker maps for the plots:
colorMarkerMap;

%% Illustration: Centralized mMIMO

% Determining the position of the BS:
dx_BS=d/2;      % x-position of the BS [m]
dy_BS=d/2;      % x-position of the BS [m]

% Determining the position of the MTDs:
dx_MTD=dmin+(d-dmin)*rand(1,K);     % x-position of the MTD [m]
dy_MTD=dmin+(d-dmin)*rand(1,K);     % y-position of the MTD [m]

% Plot the spatial distribution of BS and MTD:
fig1=figure(1);
    fig1.Position=[200 200 480 440];
    scatter(dx_BS,dy_BS,400,'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
    hold on    
    scatter(dx_MTD,dy_MTD,50,'MarkerEdgeColor',black,'LineWidth',1.5)
    grid on
    xlim([0 d])
    ylim([0 d])
    set(gca,'XTick',[], 'YTick', [])
    box on
    grid on
    % xlabel('$d_X$ [m]','Interpreter','latex','FontSize',16)
    % ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',16)
    leg=legend('BS','MTD');
    set(leg,'Interpreter','latex','FontSize',16,'Location','Southeast')
    saveas(fig1,'mMIMO.png')
    saveas(fig1,'mMIMO.eps','epsc')
    
%% Illustration: Partially Distributed mMIMO

Q=16;                                       % Number of APs
% Determining the position of the APs:
range_AP=d/(2*sqrt(Q)):d/sqrt(Q):d-d/(2*sqrt(Q));     
dx_AP=zeros(1,Q);                           % x-position of the APs [m]
dy_AP=zeros(1,Q);                           % y-position of the APs [m]
for idx=1:sqrt(Q)
   dx_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP(idx); 
end
for idx=1:sqrt(Q)
   dy_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP; 
end

% Plot the spatial distribution of APs and MTD:
fig2=figure(2);
    fig2.Position=[200 200 480 440];
    scatter(dx_AP,dy_AP,100,'MarkerEdgeColor',green,'MarkerFaceColor',green)
    hold on    
    scatter(dx_MTD,dy_MTD,50,'MarkerEdgeColor',black,'LineWidth',1.5)
    box on
    grid on
    set(gca,'XTick',[], 'YTick', [])
    ylim([0 100])
    % xlabel('$d_X$ [m]','Interpreter','latex','FontSize',16)
    % ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',16)
    leg=legend('AP','MTD');
    set(leg,'Interpreter','latex','FontSize',16,'Location','Southeast')
    saveas(fig2,'PD_mMIMO.png')
    saveas(fig2,'PD_mMIMO.eps','epsc')
    
%% Illustration: Totally Distributed mMIMO

Q=64;                                       % Number of APs
% Determining the position of the APs:
range_AP=d/(2*sqrt(Q)):d/sqrt(Q):d-d/(2*sqrt(Q));     
dx_AP=zeros(1,Q);                           % x-position of the APs [m]
dy_AP=zeros(1,Q);                           % y-position of the APs [m]
for idx=1:sqrt(Q)
   dx_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP(idx); 
end
for idx=1:sqrt(Q)
   dy_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP; 
end

% Plot the spatial distribution of APs and MTD:
fig3=figure(3);
    fig3.Position=[200 200 480 440];
    scatter(dx_AP,dy_AP,25,'MarkerEdgeColor',red,'MarkerFaceColor',red)
    hold on    
    scatter(dx_MTD,dy_MTD,50,'MarkerEdgeColor',black,'LineWidth',1.5)
    box on
    grid on
    set(gca,'XTick',[], 'YTick', [])
    % xlabel('$d_X$ [m]','Interpreter','latex','FontSize',16)
    % ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',16)
    leg=legend('AP','MTD');
    set(leg,'Interpreter','latex','FontSize',16,'Location','Southeast')
    saveas(fig3,'TD_mMIMO.png')
    saveas(fig3,'TD_mMIMO.eps','epsc')
    
%% Illustration: Partially Distributed Radio Stripes

Q=16;                                       % Number of APs
% Define the positions of the APs:
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

% Plot the geographic distribution of APs and MTD:
fig4=figure(4);
    fig4.Position=[200 200 480 440];
    scatter(dx_AP,dy_AP,100,'MarkerEdgeColor',green,'MarkerFaceColor',green)
    hold on    
    scatter(dx_MTD,dy_MTD,50,'MarkerEdgeColor',black,'LineWidth',1.5)
    box on
    grid on
    set(gca,'XTick',[], 'YTick', [])
    % xlabel('$d_X$ [m]','Interpreter','latex','FontSize',16)
    % ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',16)
    leg=legend('AP','MTD');
    set(leg,'Interpreter','latex','FontSize',16,'Location','Southeast')
    saveas(fig4,'PD_Radio_Stripes.png')
    saveas(fig4,'PD_Radio_Stripes.eps','epsc')

%% Illustration: Totally Distributed Radio Stripes

Q=64;                                       % Number of APs
% Define the positions of the APs:
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

% Plot the geographic distribution of APs and MTD:
fig5=figure(5);
    fig5.Position=[200 200 480 440];
    scatter(dx_AP,dy_AP,25,'MarkerEdgeColor',red,'MarkerFaceColor',red)
    hold on    
    scatter(dx_MTD,dy_MTD,50,'MarkerEdgeColor',black,'LineWidth',1.5)
    box on
    grid on
    set(gca,'XTick',[], 'YTick', [])
    % xlabel('$d_X$ [m]','Interpreter','latex','FontSize',16)
    % ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',16)
    leg=legend('AP','MTD');
    set(leg,'Interpreter','latex','FontSize',16,'Location','Southeast')
    saveas(fig5,'TD_Radio_Stripes.png')
    saveas(fig5,'TD_Radio_Stripes.eps','epsc')