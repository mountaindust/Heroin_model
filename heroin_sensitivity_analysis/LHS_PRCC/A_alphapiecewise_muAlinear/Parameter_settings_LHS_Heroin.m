%Need this file to set parameter baselines and timespan

% PARAMETER BASELINE VALUES

m=-0.004995197428455;%-0.004995197428455;
beta_A=0.003075241841496;%0.003075241841496;
beta_P=3.088422784921360e-04;%3.088422784921360e-04;
theta_1=0.214850659242197;%0.214850659242197;
epsilon=2.543540776849738;%2.543540776849738;
mu=0.00710; 
mu_H=0.0466; 
gamma=0.005281682020894;%0.005281682020894;
theta_2=0.665253867081153;%0.665253867081153;
sigma=0.109745817098382;%0.109745817098382;
zeta=0.188896020171160;%0.188896020171160;
theta_3=17.900486383845042;%17.900486383845042;
nu=0.002278030004504;%0.002278030004504;
omega=0.0000000001;
b=0.268584385064102;%0.268584385064102;
c=-0.027158918855258;%-0.027158918855258;
d=9.984897519858505e-04;%9.984897519858505e-04;
e=0.008708050841388;%0.008708050841388;


% %

params=[m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_H,theta_2,zeta,theta_3,nu,omega,b,c,d,e];


%% TIME SPAN OF THE SIMULATION
t_end= 6; % length of the simulations- 30 days
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 6; % time points of interest for the distemper analysis- amount of infected on last day = 30

% INITIAL CONDITIONS FOR THE ODE MODEL

P0=0.095416167677016;%0.095416167677016;
A0=0.007233998280094;%0.007233998280094;
H0=4.687664094772614e-04;%4.687664094772614e-04;
R0=0.002880151501971;%0.002880151501971;
S0=1-P0-A0-H0-R0;


total=S0+P0+A0+H0+R0;

%y0 = [S0,P0,A0,H0,R0]; 
% Variables Labels
% y_var_label={'S0','P0','A0','H0','R0'};

PRCC_var={'m','beta_A','beta_P','theta_1','epsilon','gamma','sigma','mu','mu_H','theta_2','zeta','theta_3','nu','omega','b','c','d','e','P0','A0','H0','R0'};


