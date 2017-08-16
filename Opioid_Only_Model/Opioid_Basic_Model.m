%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code solves a model of a basic opioid addiction epidemic
%
%  
%
%
%
%
% Author: Nick Battista
% Date: August 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Opioid_Basic_Model()

%
% STOCHASTIC? TIME DELAY?
%
global stochastic_flag;
stochastic_flag = 0;
time_delay_flag = 0;

%
% Temporal information
%
Tstart = 0;
Tstop = 100;

%
% initial conditions
%
S_0 = 0.9; 
G_0 = 0.10; 
R_0 = 0.00; 


%
Initial_Values = [S_0 G_0 R_0];


%
% ode45 is matlab's ode solver
%
if time_delay_flag == 0
    options=odeset('RelTol',1e-3);
	[t,sol] = ode45(@f,[Tstart Tstop],Initial_Values,options);
else
    options=odeset('RelTol',1e-3);
    [t,sol] = ode45(@f,[Tstart Tstop],Initial_Values,options); 
end


%
% storing solutions for each variable, theta_k.
% 
S  = sol(:,1);    %gives us S(t)
G =  sol(:,2);    %gives us G(t)
R  = sol(:,3);    %gives us R(t)
H = 1 - S - G - R; %gives us H(t)

%
% Plotting solutions
%
plot_Phase_Planes(S,G,H,R);
plot_Time_Evolutions(t,S,G,H,R)


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: RHS vector of the problem: this function evaluates the 
%           RHS of the ODEs and passes it back to the ode45 solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dvdt = f(t,sol)

global stochastic_flag;

%(d/dt)^2 x_n + x_n - eps (1-x_n^2) (d/dt) x_n = eps alpha (d/dt) x_i (t-tau); i\neq n

%
% ODE Coupling Parameters
%
%
% alpha = 0.1;   % S->G : people who are prescribed prescription opioids
% eps = 0.8;     % G->S : people who use their prescriptions and then go back to susceptible
% delta = 0.1;   % R->S : people who finish their treatment and then go back to susceptible class
% beta = 0.05;    % S->H : people who get opioids from their relatives/friends/etc to abuse them
% d = 0.01;       %      : natural death rate
% dSTAR = 0.02 ;  %      : enhanced death rate for opioid abusers
% gamma =(1-eps);  % G->H : percent of prescribed opioid class who get addicted to opioids
% mu = 1-delta/2;       % R->H : rate at which users in treatment fall back into drug use
% zeta = 0.9;      % H->R : rate at which Opioid abusers start treatment
% sigma = 1-delta-mu; % R->H : rate at which people in treatment fall back into use themselves.


alpha = 0.2;      % S->G : people who are prescribed prescription opioids
eps = 0.76;        % G->S : people who use their prescriptions and then go back to susceptible
beta = 0.006;     % S->H : people who get opioids from their relatives/friends/etc to abuse them
d = 0.008237;     %      : natural death rate
dSTAR = 0.00834 ; %      : enhanced death rate for opioid abusers
gamma =(1-eps);   % G->H : percent of prescribed opioid class who get addicted to opioids
zeta = 0.25;      % H->R : rate at which Opioid abusers start treatment
delta = 0.09;     % R->S : people who finish their treatment and then go back to susceptible class
mu = 1-delta/2;       % R->H : rate at which users in treatment fall back into drug use
sigma = 1-delta-mu; % R->H : rate at which people in treatment fall back into use themselves.

% NOTE: sigma+delta+mu = 1.0;
% NOTE: eps+gamma = 1.0;

%
% Components of vectors
%
S = sol(1);         % Susceptible Class
G = sol(2);         % Prescribed Opioid Class
R = sol(3);         % People in Treatment
H = 1 - S - G - R;  % Abusing Opioid Class

Lambda = d*(S+R+G) + dSTAR*H;

%
% Stochastic Piece ("white noise")
%
if stochastic_flag == 1
    white_noise = awgn(x1,1,'measured');
else
    white_noise = 0;
end


%
% ODES (RHS)
%
dS = Lambda + delta*R - alpha*S - beta*S*H + eps*G - d*S;
dG = alpha*S - gamma*G - eps*G - d*G;
dR = zeta*H - mu*R*H - delta*R - d*R - sigma*R;


%
% Vector to be evaluated
%
dvdt = [dS dG dR]';


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Phase_Planes(x1,x2,x3,x4)

figure(1)
plot(x1,x2,'r-','LineWidth',3); hold on;
plot(x1,x3,'b-','LineWidth',3); hold on;
plot(x1,x4,'k-','LineWidth',3); hold on;
xlabel('x1');
ylabel('x2,x3,x4');
legend('G vs. S','H vs. S','R vs. S');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Evolutions(t,x1,x2,x3,x4)

figure(2)
plot(t,x1,'m-','LineWidth',4); hold on;
plot(t,x2,'r-','LineWidth',4); hold on;
plot(t,x3,'b-','LineWidth',4); hold on;
plot(t,x4,'k-','LineWidth',4); hold on;
xlabel('time');
ylabel('populations');
legend('Susceptible', 'Prescribed', 'Opioid Abuse', 'Treatment');