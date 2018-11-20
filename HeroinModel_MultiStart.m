%to avoid error with "manymins" 
clear all; 
clc;

% the parameter vector we will approximate
% z=[alpha, beta, xi, theta_1, epsilon, mu, mu_A, mu_H, gamma, theta_2,
% sigma_A, zeta, theta_3, sigma_H, nu]


% The vector of the initilized parameters 
LowerBounds=[0.1  0.001  0.5   0.0001  0.8  0.001  0.001  0.002  0.0001  0.4  0.1  0.1  0.0001  0.7  0.05 0.1 0.01 0.01 .01];
%z0=[0.5 0.3 0.8 0.000008 0.00027 0.0000003 0.0000006 0.0000001 0.000001 0.2815 0.1];
UpperBounds=[0.5  0.005  0.9   0.001   2    0.009  0.009  0.01    0.001  0.9  0.5  0.5   0.001  1.0   0.5 0.3 0.2  0.3  0.3];

xstart=0.5*(LowerBounds + UpperBounds); % initial guesses for parameters (starting in the middle of each of the ranges)

% % % % %  Optimization Function fmincon  % % % % %
%x0 is xstart, objective is our specific model
problem=createOptimProblem('fmincon','x0', xstart,'objective',@HeroinModel_ODE45...
         ,'lb',LowerBounds,'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

%number of times want to do 
numstartpoints=3;

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

