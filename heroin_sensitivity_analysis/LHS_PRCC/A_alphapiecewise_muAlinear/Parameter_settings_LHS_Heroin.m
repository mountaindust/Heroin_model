%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.005595650279299;
beta_A=8.784316423507084e-04;
beta_P=6.543137174001155e-05;
theta_1=0.222457489109919;
epsilon=2.527869965595365;
mu=0.00710; 
mu_H=0.0466; 
gamma=0.005050794777625;
theta_2=0.236479520411597;
sigma=0.101518004918260;
zeta=0.198182427387906;
theta_3=19.726408301325770;
nu=5.312631489285296e-04;
omega=0.0000000001;
b=0.270110337915851;
c=-0.026969098706352;
d=9.774825266577510e-04;
e=0.008831387924813;



% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_H,theta_2,zeta,theta_3,nu,omega,b,c,d,e];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

P0=0.094972745098928;
A0=0.007097422873023;
H0=4.648950554349272e-04;
R0=0.005072290169507;
S0=1-P0-A0-H0-R0;


total=S0+P0+A0+H0+R0;

%y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
% y_var_label={'S0','P0','A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_H','theta_2','zeta','theta_3','nu','omega','b','c','d','e','P0','A0','H0','R0'};


