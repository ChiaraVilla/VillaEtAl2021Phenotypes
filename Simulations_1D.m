%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% "Evolutionary Dynamics in Vascularised Tumours under Chemotherapy:  %%%
%%%           Mathematical Modelling, Asymptotic Analysis               %%%
%%%                   and Numerical Simulations"                        %%%
%%%                                                                     %%%
%%%             C. Villa (*), M.A.J. Chaplain, T. Lorenzi               %%%
%%%                                                                     %%%
%%%              Vietnamese Journal of Mathematics (2021)               %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%% (*) Email: chiara.villa.1[at]sorbonne-universite.fr                 %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   Code for 1D simulations under given distributions of oxygen and   %%%
%%%   chemotherapeutic agent concentrations (see Figures 25 in          %%%
%%%   villa2021evolutionary), solving system (1) together with          %%%
%%%   definitions (2), (3), (9) and under initial conditions (20).      %%%
%%%   The solution is compared with that of system (23), together       %%%
%%%   with definitions (11)-(13).                                       %%%
%%%   For the implementation of the numerical scheme described in the   %%%
%%%   paper for the 2D problem under dynamical concentrations of oxygen %%%
%%%   and chemotherapeutic agent concentrations, see code available for %%%
%%%   villa2021modelling (Simulations_2D.m in the same repository).     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                       %
%%% Simulations_1D.m: simulates a 1D mathematical model of cancer cell    %
%%% evolutionary dynamics in vascularised tumours under chemotherapy      %
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
%%% Disclaimer: the code is not optimised (e.g. changes in choice of    %%%
%%% points r1-3 will require changes in various parts of the code, ...) %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')

%% Time discretisation
tf = 1000000;                 % Final time
dt = 0.1;                     % Timestep size
t = 0:dt:tf;                  % Time span
Nt = tf/dt;                   % Number of timesteps

%% Spacial discretisation
l1 = 0;                       % 2D domain [l1,l2]x[l1,l2]
l2 = 0.05;                     % 2D domain [l1,l2]x[l1,l2]
Nr = 101;                     % Number or spatial grid points
r = linspace(l1,l2,Nr);       % Discretised spatial domain (1D)

%% Phenotypic discretisation 
ym = -7;                      % Phenotypic domain [ym,yM]
yM = 7;                       % Phenotypic domain [ym,yM]  
ny = 1000;                    % Number of phenotypic grid points
y = linspace(ym,yM,ny);       % Discretised phenotypic domain
dy = abs(y(2)-y(1));          % Phenotypic grid step size

%% Parameters
beta = 10^(-6);               % Phenotypic diffusion
gammas = 10^(-4);             % Max aerobic prolif rate
phi = 10^(-5);                % Max anaerobic prolif rate
d = 2*10^(-13);               % Cell death rate
mms = 1.5*10^(-7);            % Michaelis Menten constant for O2
gammac = 1.8*10^(-4);         % Max death rate by drug
mmc = 10^(-6);                % Michaelis Menten constant for drug
v0 = 1;                       % Initial inverse variance
mu0 = 0.5;                    % Initial mean
rho0 = 10^(8);                % Initial cell density


%% Initial conditions
n0 = exp(-(v0/2)*(y-mu0).^2); 
n0 = rho0*n0/(sum(n0,2)*dy);  % Initial phenotypic distribution
% Store initial condition at three points in the domain for comparison
v = [v0,v0,v0];
mu = [mu0,mu0,mu0]; 
rho = [rho0,rho0,rho0];
n = [n0;n0;n0];
rho2 = sum(n,2)*dy;           % Cell density computed directly from n

%% Given oxygen and chemotherapeutic agent concentrations
% Chosen to fit experimentally obtained profiles in Helmlinger (1997)
% (cf. Figure 2 in ville2021evolutionary)
% Oxygen concentration profile
sr = (4.5*1e-8)*22.6*exp(-4885*(r+0.01).^2);
% Drug concentration profile (multiply by appropriate factor for Figure 5)
cr = 10^(-5)*exp(-6885*(r+0.01).^2);
%Select three points for simulations
r3 = [r(15),r(31),r(71)]; 
s = [sr(15),sr(31),sr(71)];
c = [cr(15),cr(31),cr(71)];
 
%% Initialise arrays to store evolutions of rho, mu, sigma over time
rhon1 = [rho0];
rhon2 = rhon1;
rhon3 = rhon1;
rhoc1 = rhon1;
rhoc2 = rhon1;
rhoc3 = rhon1;
mun1 = [mu0];
mun2 = [mu0];
mun3 = [mu0];
muc1 = [mu0];
muc2 = [mu0];
muc3 = [mu0];
sign1 = [1/v0];
sign2 = [1/v0];
sign3 = [1/v0];
sigc1 = [1/v0];
sigc2 = [1/v0];
sigc3 = [1/v0];
time = [t(1)];

%% Time iteration
% Given a, b, h - defintions (11)-(13) 
ac = gammac*c./(mmc+c);
as = gammas*s./(mms+s);
a2 = as -ac + ((phi+ac).^2)./(phi+as+ac);
h2 = (phi+ac)./(phi+as+ac);
b2 = as + ac + phi;
figure('Units','normalized','Position',[0 0 1 1])
for i=1:length(t) 
    
    %%% Numerical solutions
    for j=1:3
        % Solve original PDE - eq. (1)_1 
        R(j,1:ny) = a2(j) - b2(j).*(y(1:ny)-h2(j)).^2 - d*rho2(j);  
        der2y = (n(j,1:ny-2)+n(j,3:ny)-2*n(j,2:ny-1))/(dy^2);
        n(j,2:ny-1) = n(j,2:ny-1) + dt*(R(j,2:ny-1).*n(j,2:ny-1) +beta*der2y);
        % Solve ODE system - system (23) 
        v(j) = v(j) + dt*2*(b2(j)-beta*(v(j))^2);
        mu(j) = mu(j) + dt*(2/v(j))*(h2(j)*b2(j)-mu(j)*b2(j));
        rho(j) = rho(j) + dt*rho(j)*(-b2(j)/v(j)+a2(j)-b2(j)*(mu(j)-h2(j))^2 -d*rho(j));
    end
    % Zero flux boundary conditions 
    % (on numerically imposed boundary of the phenotypic domain)
    n(:,1) = n(:,2);
    n(:,ny) = n(:,ny-1); 
    % Update local cell density - eq. (1)_2
    rho2 = sum(n,2)*dy; 

    %%% Save and plot
    if mod(t(i),5000)==0
        % Construct n from solution of ODE system
        nc = exp(-(v(1)/2)*(y-mu(1)).^2); 
        nc1 = rho(1)*nc/(sum(nc,2)*dy); 
        nc = exp(-(v(2)/2)*(y-mu(2)).^2); 
        nc2 = rho(2)*nc/(sum(nc,2)*dy);
        nc = exp(-(v(3)/2)*(y-mu(3)).^2); 
        nc3 = rho(3)*nc/(sum(nc,2)*dy);
        % Save variables to plot over time rho, mu, sigma
        rhon1 = [rhon1,rho2(1)];
        rhon2 = [rhon2,rho2(2)];
        rhon3 = [rhon3,rho2(3)];
        rhoc1 = [rhoc1,rho(1)];
        rhoc2 = [rhoc2,rho(2)];
        rhoc3 = [rhoc3,rho(3)];
        munn1 = sum(y.*n(1,:))*dy/rho2(1);
        munn2 = sum(y.*n(2,:))*dy/rho2(2);
        munn3 = sum(y.*n(3,:))*dy/rho2(3);
        mun1 = [mun1,munn1];
        mun2 = [mun2,munn2];
        mun3 = [mun3,munn3];
        muc1 = [muc1,mu(1)];
        muc2 = [muc2,mu(2)];
        muc3 = [muc3,mu(3)];
        sign1 = [sign1,sum((y.^2).*n(1,:))*dy/rho2(1)-munn1^2];
        sign2 = [sign2,sum((y.^2).*n(2,:))*dy/rho2(2)-munn2^2];
        sign3 = [sign3,sum((y.^2).*n(3,:))*dy/rho2(3)-munn3^2];
        sigc1 = [sigc1,1/v(1)];
        sigc2 = [sigc2,1/v(2)];
        sigc3 = [sigc3,1/v(3)];
        time = [time,t(i)];
        % Plot
        clf
        subplot(2,3,1)
        plot(time*10^(-4),rhoc1,'r')
        hold on
        plot(time*10^(-4),rhoc2,'b')
        hold on
        plot(time*10^(-4),rhoc3,'Color',[0.04 0.75 0.01])
        legend({'$x=0.007$','$x=0.015$','$x=0.035$'},'Interpreter','latex')
        hold on
        plot(time*10^(-4),rhon1,'--k')
        hold on
        plot(time*10^(-4),rhon2,'--k')
        hold on
        plot(time*10^(-4),rhon3,'--k')
        title('$\rho(t)$')
        xlim([0,100])
        axis square
        subplot(2,3,2)
        plot(time*10^(-4),muc1,'r')
        hold on
        plot(time*10^(-4),muc2,'b')
        hold on
        plot(time*10^(-4),muc3,'Color',[0.04 0.75 0.01])
        legend({'$x=0.007$','$x=0.015$','$x=0.035$'},'Interpreter','latex')
        hold on
        plot(time*10^(-4),mun1,'--k')
        hold on
        plot(time*10^(-4),mun2,'--k')
        hold on
        plot(time*10^(-4),mun3,'--k')
        title('$\mu(t)$')
        axis square
        xlim([0,100])
        ylim([0,1])
        subplot(2,3,3)
        plot(time*10^(-4),sigc1,'r')
        hold on
        plot(time*10^(-4),sigc2,'b')
        hold on
        plot(time*10^(-4),sigc3,'Color',[0.04 0.75 0.01])
        legend({'$x=0.007$','$x=0.015$','$x=0.035$'},'Interpreter','latex')
        hold on
        plot(time*10^(-4),sign1,'--k')
        hold on
        plot(time*10^(-4),sign2,'--k')
        hold on
        plot(time*10^(-4),sign3,'--k')
        title('$\sigma^2(t)$')
        xlim([0,100])
        axis square
        subplot(2,3,4)
        plot(y,n(1,:),'r')
        xlabel('$y$')
        axis square
        hold on
        plot(y,nc1,'--k')
        xlim([-4,4])
        title(['$n($',num2str(round(t(i)*10^(-4),1)),'$,0.007,y)$'])
        subplot(2,3,5)
        plot(y,n(2,:),'b')
        xlabel('$y$')
        axis square
        hold on
        plot(y,nc2,'--k')
        xlim([-4,4])
        title(['$n($',num2str(round(t(i)*10^(-4),1)),'$,0.015,y)$'])
        subplot(2,3,6)
        plot(y,n(3,:),'Color',[0.04 0.75 0.01])
        hold on
        plot(y,nc3,'--k')
        axis square
        xlabel('$y$')
        xlim([-4,4])
        title(['$n($',num2str(round(t(i)*10^(-4),1)),'$,0.035,y)$'])
        drawnow
    end
end
