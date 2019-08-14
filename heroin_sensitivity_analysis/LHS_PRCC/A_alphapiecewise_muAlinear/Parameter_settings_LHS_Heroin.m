%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.003132770417329;
beta_A=0.004579901231381;
beta_P=3.725469731616883e-04;
theta_1=5.053609307592831e-04;
epsilon=2.519035790399939;
mu=0.00868;      
mu_H=0.0507;
gamma=0.001389474357184;
theta_2=1.178529767918704;
sigma=0.035409435842082;
zeta=0.419579274859810;
theta_3=2.371839014674933;
nu=1.726536501674095e-04;
omega=0.0000000001;
b=0.281054576347342;
c=-0.030243908757160;
d=0.002942643755155;
e=0.009160945129219;
% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_H,theta_2,zeta,theta_3,nu,omega,b,c,d,e];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

P0=0.100203815696228;
A0=0.005627109256495;
H0=4.043433686669488e-04;
R0=0.066239648842513;
S0=1-P0-A0-H0-R0;

total=S0+P0+A0+H0+R0;

%y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
% y_var_label={'S0','P0','A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_H','theta_2','zeta','theta_3','nu','omega','b','c','d','e','P0','A0','H0','R0'};


