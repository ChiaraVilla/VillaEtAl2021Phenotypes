%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%    "Modelling phenotypic heterogeneity in vascularised tumours"     %%%
%%%                                                                     %%%
%%%              C. Villa (*), M.A.J. Chaplain, T. Lorenzi              %%%
%%%                                                                     %%%
%%%             SIAM Journal on Applied Mathematics (2021)              %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%% (*) Email: cv23[at]st-andrews.ac.uk                                 %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   Code for 2D simulations with random blood vessel distribution     %%%
%%%   (see Figure 1 in villa2021modelling), solving systems (2.3)       %%%
%%%   and (2.9) together with definitions (2.4), (2.5), (2.10), (4.1),  %%%
%%%   (4.4) and initial conditions (4.5).                               %%%
%%%                                                                     %%%
%%%   Analogous code can be employed to obtain the 2D results presented %%%   
%%%   in 'Villa, Chaplain, Lorenzi (2021), Evolutionary dynamics of     %%%
%%%   vascularised tumours under chemotherapy: mathematical modelling,  %%%
%%%   analysis and numerical simulations, Viet. J. Math.' by:           %%%
%%%   1) setting beta=0; 2) compute the chemotherapeutic agent          %%%
%%%   concentration in the way as the oxygen concentration below.       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                       %
%%% Simulations_2D.m: simulates a 2D mathematical model of cancer cell    %
%%% evolutionary dynamics in vascularised tumours (villa2021modelling)    %
%%% Copyright (C) 2021 C. Villa                                           %
%%%                                                                       %
%%% This program is free software: you can redistribute it and/or modify  %
%%% it under the terms of the GNU General Public License as published by  %
%%% the Free Software Foundation, either version 3 of the License, or     %
%%% (at your option) any later version.                                   %
%%%                                                                       %
%%% This program is distributed in the hope that it will be useful,       %
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%%% GNU General Public License for more details.                          %
%%%                                                                       %
%%% You should have received a copy of the GNU General Public License     %
%%% along with this program.  If not, see <https://www.gnu.org/licenses/>.%
%%%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%    Disclaimer: the code is not fully optimised (e.g. changes in     %%%
%%%    domain will require changes in the plot function code, ...)      %%% 
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 16)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')

%% Time discretisation
tf = 400000;                  % Final time
dt = 0.001;                   % Timestep size
time = 0:dt:tf;               % Time span
Nt = tf/dt;                   % Number of timesteps

%% Spacial discretisation
l1 = 0;                       % 2D domain [l1,l2]x[l1,l2]
l2 = 0.5;                     % 2D domain [l1,l2]x[l1,l2]
Nr = 100;                     % Number or spatial grid points
r = linspace(l1,l2,Nr);       % Discretised spatial domain (1D)
dr = abs(r(2)-r(1));          % Spatial grid step size

%% Phenotype discretisation
xm = 0;                       % Phenotypic domain [xm,xM]
xM = 1;                       % Phenotypic domain [xm,xM]  
Nx = 100;                     % Number of phenotypic grid points
x = linspace(xm,xM,Nx);       % Discretised phenotypic domain
dx = abs(x(2)-x(1));          % Phenotypic grid step size

%% Parameter values 
phi = 10^(-5);                % Maximum rate of cell division via anaerobic pathways
gamma = 10^(-4);              % Maximum rate of cell division via aerobic pathways
alphas = 1.5*10^(-7);         % Michaelis-Menten constant of oxygen
S0 = 6.3996*10^(-7);          % Constant rate of inflow of oxygen through blood vessels
d = 2*10^(-14);               % Rate of cell death due to competition for space
betas = 2*10^(-5);            % Oxygen diffusivity
etas = 2*10^(-11);            % Conversion factor for cell consumption of oxygen
lambdas = 2.78*10^(-6);       % Rate of natural decay of oxygen
theta = 10^(-13);             % Rate of spontaneous phenotypic changes
beta = 10^(-13);              % Cell motility
sig0 = 0.1;                   % Initial variance sigma^2 
xbar0 = 0.5;                  % Initial mean trait
rho0 = 10^8;                  % Scaling factor for the initial population density 
bv_num = 6;                   % Number of blood vessels

%% Initial conditions
% Cell phenotypic distribution with Gaussian profile   
n = initial_population(Nr,Nx,x,sig0,xbar0,rho0);
I = sum(n,3)*dx;              % Local cell density
chi = bv_rand(Nr,bv_num);     % Randomly placed blood vessels
s = S0*chi;                   % Initial oxygen distribution

%% Time interation
figure('Units','normalized','Position',[0 0 1 1]) 
for i=1:Nt
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dynamics of tumour cells - system (2.3) in villa2021modelling
    for j=2:Nr-1  
        for k=2:Nr-1
            % Fitness function R - eq. (2.4)-(2.5) and (4.1)
            gjk = gamma*(s(j,k)/(alphas+s(j,k)))*(1-x.^2); 
            R(j,k,:) = phi*(1-(1-x).^2) + gjk - d*I(j,k);
            % Second order derivatives in space and in the phenotypic dimension
            der2r = (n(j-1,k,:)+n(j+1,k,:)+n(j,k-1,:)+n(j,k+1,:)-4*n(j,k,:))/dr^2;
            der2x = (n(j,k,1:Nx-2)-2*n(j,k,2:Nx-1)+n(j,k,3:Nx))/dx^2;
            % Phenotypic distribution in Omega - eq. (2.3)_1
            n(j,k,2:Nx-1) = n(j,k,2:Nx-1) + dt*(R(j,k,2:Nx-1).*n(j,k,2:Nx-1) + theta*der2x + beta*der2r(2:Nx-1));
        end
    end
    % Update local cell density - eq. (2.3)_2
    I=sum(n,3)*dx;
    % Zero flux boundary conditions - eq. (2.3)_3 and (2.3)_4
    n = BC_zeroflux_n(n,Nx,Nr); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dynamics of oxygen - system (2.9) in villa2021modelling
    % Nonlocal term
    T = sum(gjk.*n,3)*dx;         
    % Second order derivative in space
    sder2r = (s(1:Nr-2,2:Nr-1)+s(3:Nr,2:Nr-1)+s(2:Nr-1,1:Nr-2)+s(2:Nr-1,3:Nr)-4*s(2:Nr-1,2:Nr-1))/dr^2;
    % Oxygen concentration in Omega - eq. (2.9)_1
    s(2:Nr-1,2:Nr-1) = s(2:Nr-1,2:Nr-1)+(-lambdas*s(2:Nr-1,2:Nr-1)+betas*sder2r-etas*T(2:Nr-1,2:Nr-1)+S0*chi(2:Nr-1,2:Nr-1))*dt;
    % Zero flux boundary conditions - eq. (2.9)_2
    s = BC_zeroflux_s(s,Nr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save and plot solution
    if mod(i,500)==0 
        disp(['Time = ' num2str(i*dt) 's']); % Or i*dt/(3600*24) days
        xbar = mean_trait(n,Nr,Nx,x,I,dx);   % Local mean phenotypic state
        save('saved_var.mat','chi','I','s','xbar','n','i'); 
        plot_solution(x,dx,xbar,Nr,Nx,s,I,gamma,alphas,phi,d,r,n,chi)
    end
    
end


%% Annexed functions

%%% Initial cell phenotypic distribution with Gaussian profile   
function n = initial_population(Nr,Nx,x,sig0,xbar0,rho0);
    for i=1:Nr
        for j=1:Nr
            n(i,j,1:Nx) = rho0*(exp(-(x-xbar0).^2/sig0));
        end
    end 
end

%%% Mean phenotype at each point in space  
function xbar = mean_trait(n,Nr,Nx,x,I,dx);
    for i=1:Nr
        for j=1:Nr
            nsup(1:Nx) = n(i,j,:);
            xbar(i,j) = sum(x.*nsup)*dx/I(i,j);
        end
    end
end

%%% Argmax of n at each point in space
function xmax = argmax(n,Nx,Nr,dx);
    xmax = zeros(Nr);
    for i=1:Nr
        for j=1:Nr
            nij(1:Nx) = n(i,j,1:Nx);
            temp = find(nij==max(nij));
            xmax(i,j) = dx*temp(1);
        end
    end
end

%%% Randomly selected positions of blood vessels
function chi = bv_rand(Nr,bv_num);
    chi = zeros(Nr); 
    pos = randperm(Nr^2,bv_num);
    chi(pos) = 1;
end

%%% Zero-flux boundary conditions on n
function n = BC_zeroflux_n(n,Nx,Nr);
    n(:,:,1) = n(:,:,2);
    n(:,:,Nx) = n(:,:,Nx-1);
    n(1,:,:) = n(2,:,:);
    n(Nr,:,:) = n(Nr-1,:,:);    
    n(:,1,:) = n(:,2,:);
    n(:,Nr,:) = n(:,Nr-1,:);
end 

%%% Zero-flux boundary conditions on s
function s = BC_zeroflux_s(s,Nr);
    s(Nr,:)=s(Nr-1,:);
    s(1,:)=s(2,:);
    s(:,1)=s(:,2);
    s(:,Nr)=s(:,Nr-1);
end 

%%% Plot solution
function plot_solution(x,dx,xbar,Nr,Nx,s,I,gamma,alphas,phi,d,r,n,chi)
    
    % 2D plots
    clf
    [chi1,chi2]=find(chi==1);
    chi1p=chi1/(2*Nr);
    chi2p=chi2/(2*Nr);
    subplot(2,4,1)
    sz=8;
    scatter(chi2p,chi1p,sz,'filled','r')
    xlim([r(1) r(end)])
    ylim([r(1) r(end)])
    xticks([0 0.5])
    xlabel('$r_1$')
    yticks([0 0.5])
    ylabel('$r_2$') %'rot', 0
    title(['$\omega$']) %Blood vessels
    shading flat
    view(0,90)
    axis square
    hold on
    re = rectangle('Position',[0 0 0.5 0.5],'LineWidth',1);
    subplot(2,4,2)
    pcolor(r,r,s)
    colorbar
    caxis([0 max(max(s))])
    title(['$s($T$,r_{1}$,$r_{2})$'])
    shading flat
    view(0,90)
    xticks([0 0.5])
    xlabel('$r_1$')
    yticks([0 0.5])
    ylabel('$r_2$') %'rot', 0
    axis square
    subplot(2,4,3)
    pcolor(r,r,I)
    colorbar
    caxis([0 max(max(I))])
    title(['$I($T$,r_{1}$,$r_{2})$'])
    shading flat
    view(0,90)
    xticks([0 0.5])
    xlabel('$r_1$')
    yticks([0 0.5])
    ylabel('$r_2$') %'rot', 0
    axis square
    subplot(2,4,4)
    pcolor(r,r,xbar)
    colorbar('XTick',0:0.5:1)
    caxis([0 1])
    shading flat
    xticks([0 0.5])
    xlabel('$r_1$')
    yticks([0 0.5])
    ylabel('$r_2$') %'rot', 0
    view(0,90)
    title(['$X($T$,r_{1},r_{2})$'])
    axis square
    
    % 1D cross section to compare with analytic solution
    ss(1:Nr)=s(chi1(1),1:Nr);
    II(1:Nr)=I(chi1(1),1:Nr);
    xxbar(1:Nr)=xbar(chi1(1),1:Nr);
    % Compute asymptotic solution of formal analysis - eq. (3.16)
    % (NOTE: solutions match at equilibrium, not earlier)
    As = gamma*s./(alphas + s);
    Ibar = (As+phi.^2./(phi+As))/d;  
    xbarend = phi./(phi+As);
    IIb(1:Nr)=Ibar(chi1(1),1:Nr);
    xxbare(1:Nr)=xbarend(chi1(1),1:Nr);
    % Plot
    subplot(2,4,5)
    pcolor(r,r,s)
    title(['$s($T$,r_{1}$,$r_{2})$'])
    shading flat
    view(0,90)
    xticks([0 0.5])
    xlabel('$r_1$')
    yticks([0 0.5])
    ylabel('$r_2$') %'rot', 0
    line(r,chi1p(1)*ones(1,length(r)),'Color','white','LineWidth',2,'LineStyle','--')
    axis square
    subplot(2,4,6)
    plot(r,ss)
    xlim([0 0.5])
    title(['$s($T$,r_{1}$,0.4)'])
    xticks([0 0.5])
    xlabel('$r_1$')
    axis square
    subplot(2,4,7)
    plot(r,II,'b')
    xlim([0 0.5])
    hold on
    plot(r,IIb,'r--')
    title(['$I($T$,r_{1},$0.4)'])
    xticks([0 0.5])
    xlabel('$r_1$')
    axis square
    subplot(2,4,8)
    plot(r,xxbar,'b')
    xlim([0 0.5])
    hold on
    plot(r,xxbare,'r--')
    title(['$X($T$,r_{1},0.4)$'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Code to compare xbar^infinity with argmax(n) instead 
    %%% (cf. insets in Figure 1 in villa2021modelling)
    % max = argmax(n,Nx,Nr,dx);   
    % xxmax(1:Nr)=xmax(chi1(1),1:Nr);
    % hold on
    % plot(r,xxmax,'g--')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xticks([0 0.5])
    xlabel('$r_1$')
    axis square
    drawnow
end