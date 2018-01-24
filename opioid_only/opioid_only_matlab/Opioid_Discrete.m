%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code solves a model of a basic opioid addiction epidemic using
% discrete dynamics
%
%
% Author: Nick Battista
% Date: August 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Opioid_Discrete()


%
% Temporal information
%
Tstop = 25;
t=0:1:Tstop;

%
% initial conditions
%
S(1) = 0.95; 
G(1)= 0.05; 
R(1)= 0.00; 
H(1) = 0.0;


%
% Dynamical Coupling Parameters
%
alpha = 0.9;        % S->G : people who are prescribed prescription opioids
eps = 0.74;         % G->S : people who use their prescriptions and then go back to susceptible
beta = 0.006;       % S->H : people who get opioids from their relatives/friends/etc to abuse them
d = 0.00824;       %      : natural death rate
dSTAR = 0.00834 ;   %      : enhanced death rate for opioid abusers
gamma =(1-eps);     % G->H : percent of prescribed opioid class who get addicted to opioids
zeta = 0.75;        % H->R : rate at which Opioid abusers start treatment
delta = 0.09;       % R->S : people who finish their treatment and then go back to susceptible class
mu = 1-delta/2;     % R->H : rate at which users in treatment fall back into drug use
sigma = 1-delta-mu; % R->H : rate at which people in treatment fall back into use themselves.


%
% Time Evolution Step
%
for i=1:Tstop 
   Lambda = d*(S(i)+G(i)+R(i)) + dSTAR*H(i); 
   S(i+1) = S(i) + Lambda + delta*R(i) - alpha*S(i) - beta*S(i)*H(i) + eps*G(i) - d*S(i);
   G(i+1) = G(i) + alpha*S(i) - gamma*G(i) - eps*G(i) - d*G(i);
   R(i+1) = R(i) + zeta*H(i) - mu*R(i)*H(i) - delta*R(i) - d*R(i) - sigma*R(i);
   H(i+1) = 1-S(i+1)-G(i+1)-R(i+1);
end


%
% Plotting solutions
%
%plot_Phase_Planes(S,G,H,R);
plot_Time_Evolutions(t,S,G,H,R)


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

