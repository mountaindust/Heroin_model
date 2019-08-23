%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.00597;%-0.005967450875825;
beta_A=0.00479;%0.004792878305594;
beta_P=0.00139;%0.001385562582746;
theta_1=0.0686;%0.068622166795905;
epsilon=2.520;%2.520349321095830;
mu=0.00868;
mu_H=0.0507;
gamma=0.00238;%0.002382471532991;
theta_2=0.356;%0.356456354919115;
sigma=0.0278;%0.027757623491680;
zeta=0.474;%0.473812996408801;
theta_3=1.87;%1.867973157342380;
nu=0.000462;%4.621392068173277e-04;
omega=0.0000000001;
b=0.295;%0.294801569102472;
c=-0.0298;%-0.029767190605917;
d=0.00305;%0.003050773511563;
e=0.00952;%0.009515476679820;
% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_H,theta_2,zeta,theta_3,nu,omega,b,c,d,e];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

P0=0.0937;%0.093687513067027;
A0=0.00543;%0.005426282173042;
H0=0.000408;%4.080401954644302e-04;
R0=0.0861;%0.086079957903333;
S0=1-P0-A0-H0-R0;

total=S0+P0+A0+H0+R0;

%y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
% y_var_label={'S0','P0','A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_H','theta_2','zeta','theta_3','nu','omega','b','c','d','e','P0','A0','H0','R0'};


