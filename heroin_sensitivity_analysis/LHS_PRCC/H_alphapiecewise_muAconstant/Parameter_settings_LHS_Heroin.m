%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.00541;
beta_A=0.000193;
beta_P=0.0000138;
theta_1=0.0984;
epsilon=2.50;
gamma=0.00101;
sigma=1.14;
mu=0.00868;
mu_A=0.0109;
mu_H=0.0507;
theta_2=1.98;
zeta=0.199;
theta_3=3.77;
nu=0.000309;
omega=0.0000000001;
b=0.267;
c=-0.0266;
% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_A,mu_H,theta_2,zeta,theta_3,nu,omega,b,c];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

P0=0.0918;
A0=0.00557;
H0=0.000366;
R0=0.00240;
S0=1-P0-A0-H0-R0;

total=S0+P0+A0+H0+R0;

%y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
% y_var_label={'S0','P0','A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_A','mu_H','theta_2','zeta','theta_3','nu','omega','b','c','P0','A0','H0','R0'};


