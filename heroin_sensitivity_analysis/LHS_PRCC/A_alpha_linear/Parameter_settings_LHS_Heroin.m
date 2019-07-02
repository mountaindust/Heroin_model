%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.0156;
beta_A=0.00235;
beta_P=0.000141;
theta_1=0.000507;
epsilon=2.54;
gamma=0.00115;
sigma=0.0284;
mu=0.00868;
mu_A=0.00870;
mu_H=0.0507;
theta_2=0.0370;
zeta=0.265;
theta_3=3.51;
nu=0.00657;
omega=0.0000000001;
b=0.303;
% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_A,mu_H,theta_2,zeta,theta_3,nu,omega,b];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

S0=0.858016;
P0=0.0835;
A0=0.00671;
H0=0.000874;
R0=0.0509;

total=S0+P0+A0+H0+R0;

%y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
%y_var_label={'S0', 'P0', 'A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_A','mu_H','theta_2','zeta','theta_3','nu','omega','b','P0','A0','H0','R0'};

