%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.00559565027929907;
beta_A=0.000878431642350708;
beta_P=6.54313717400116e-05;
theta_1=0.222457489109919;
epsilon=2.52786996559537;
mu=0.00710; 
mu_H=0.0466; 
gamma=0.00505079477762453;
theta_2=0.236479520411597;
sigma=0.101518004918260;
zeta=0.198182427387906;
theta_3=19.7264083013258;
nu=0.000531263148928530;
omega=0.0000000001;
b=0.270110337915851;
c=-0.0269690987063522;
d=0.000977482526657751;
e=0.00883138792481281;


% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_H,theta_2,zeta,theta_3,nu,omega,b,c,d,e];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

P0=0.0949727450989279;
A0=0.00709742287302280;
H0=0.000464895055434927;
R0=0.00507229016950725;
S0=1-P0-A0-H0-R0;


total=S0+P0+A0+H0+R0;

y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
y_var_label={'S0','P0','A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_H','theta_2','zeta','theta_3','nu','omega','b','c','d','e','P0','A0','H0','R0'};


