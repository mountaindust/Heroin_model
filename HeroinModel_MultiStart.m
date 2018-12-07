%File name: HeroinModel_MultiStart.m
%to avoid error with "manymins" 
clear all; 
clc;



% the parameter vector we will approximate 
%        z =[alpha  beta_A    beta_P   theta_1   epsilon  gamma   theta_2   sigma_A    zeta    theta_3  sigma_H    nu       P0      A0        H0       R0  ]
LowerBounds=[0.01    0.00001  0.00001   0.00001    0.8   0.00001  0.00002   0.00001   0.00001  0.00004  0.00001  0.00001   0.01   0.00001   0.00001  0.00001];
UpperBounds=[0.7      0.1      0.1       0.1        8      0.1      0.2        2       2        0.4        2        2      0.8      0.2       0.1      0.3];

% For now: alpha guess from opioid paper since higher in TN; beta_A guess from opioid paper; 
% beta_P guess from opioid paper; theta_1 complete guess, assume smaller than beta_A, beta_P; 
% epsilon 0.8-8 from opioid paper; gamma guess from opioid paper; theta_2
% guess twice as large as theta_1; sigma_A guess small because
% "successful recovereds"; zeta will be smaller than opioid paper because less "relapse" but still
% put opioid largest value in; theta_3 guess four times as large as
% theta_1; sigma_H guessed same as sigma_A; nu guess same as zeta; P0
% guess; A0 guess; H0 guess; R0 guess 


xstart=0.5*(LowerBounds + UpperBounds); % initial guesses for parameters (starting in the middle of each of the ranges)

% % % % %  Optimization Function fmincon  % % % % %
%x0 is xstart, objective is our specific model
problem=createOptimProblem('fmincon','x0', xstart,'objective',@HeroinModel_ODE45...
         ,'lb',LowerBounds,'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

%number of times want to do 
numstartpoints=4;

ms=MultiStart('Display', 'iter'); % Define a multistart problem

%manymins is a vector of solution objects (obtain multiple solutions), runs the multistart 
[x,fval, exitflag, output, manymins]=run(ms,problem,numstartpoints);

global ModelParameters

for i=1: length(manymins)
    ModelParameters(i,:)=manymins(i).X;
end

for i=1: length(manymins) 
    fval(i)=manymins(i).Fval;
end 
fval=fval';

for i=1: length(manymins) 
    EF(i)=manymins(i).Exitflag;
end
EF=EF';

beep on;
beep
beep
beep
beep
beep
beep

