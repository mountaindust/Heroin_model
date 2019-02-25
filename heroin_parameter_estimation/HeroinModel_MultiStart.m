% File name: HeroinModel_MultiStart.m (used to be in Heroin_model folder)

% Avoids error with "manymins," if needed:

clear all; 
clc;

% The parameter vector z we will approximate 
%        z =[alpha  beta_A    beta_P     theta_1   epsilon    gamma   theta_2    sigma     zeta    theta_3      nu       P0      A0        H0        R0  ]
%WORKED/VERY CLOSE WITH THESE BOUNDS for z0=[0.15  0.00094  0.00266   0.0001   2.5   0.00744   0.0002   0.5  0.05   0.0004   0.05 0.1 0.0057 0.0013 0.009];
%LowerBounds=[0.1    0.0008     0.001     0.00001      1      0.001     0.0001   0.1      0.01      0.0001    0.01      0.01    0.001     0.001   0.001];
%UpperBounds=[0.2    0.00099    0.003     0.0002       3      0.009     0.0003    0.9       0.2      0.0009     0.1       0.3    0.009     0.009      0.2];

LowerBounds=[0.01    0.00001    0.0001     0.00001      0.8     0.001     0.00001   0.01      0.01    0.0001    0.01     0.001   0.0001     0.00001  0.00001];
UpperBounds=[0.7       0.1      0.009      0.1          4       0.1         0.3      2         1       0.6       1       0.5     0.2         0.1      0.1];


% Bound choices in past:
% alpha: guess from opioid paper since higher in TN
% beta_A: guess from opioid paper
% beta_P: guess from opioid paper
% theta_1: complete guess, assume smaller than beta_A, beta_P
% epsilon: 0.8-8 from opioid paper
% gamma: guess from opioid paper
% theta_2: guess twice as large as theta_1
% sigma: guess small because those in a stable recovery state
% zeta: will be smaller than opioid paper because less "relapse" but still put opioid largest value in
% theta_3: guess four times as large as theta_1;
% nu: guess same as zeta
% P0: guess
% A0: guess 
% H0: guess
% R0: guess 

% Initial guesses for parameters (starting in the middle of each of the ranges)
xstart=0.5*(LowerBounds + UpperBounds); 

% % % % %  Optimization Function fmincon  % % % % %
% x0 is xstart, objective is what we are trying to minimize which comes from 
% value = HeroinModel_ODE45(z) = f(x) as output
problem=createOptimProblem('fmincon','x0', xstart,'objective',@HeroinModel_ODE45,...
         'lb',LowerBounds,'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

% Number of times I want to run optimization scheme
numstartpoints=2;

% Define a multistart problem; results are reported after each local solver run, in addition to the final summary
ms=MultiStart('Display', 'iter'); 

% Manymins is a vector of solutions containing the distinct local minima found during the run;
%  runs MultiStart on numstartpoints to find a solution or multiple local solutions to problem
[x,fval, exitflag, output, solutions]=run(ms,problem,numstartpoints);

global ModelParameters

%each component of manymins stores a vector of the parameters that were used each time step,
%ModelParameters puts them into a matrix 
for i=1: length(solutions)
    ModelParameters(i,:)=solutions(i).X;
end

%each time step has an fval; takes Fval values that are stored in manymins and
%creates vector out of them
for i=1: length(solutions) 
    fval(i)=solutions(i).Fval;
end 
fval=fval';

%each time step has an exitflag value; takes Exitflag values that are stored in
%manymins and creates vector out of them
for i=1: length(solutions) 
    EF(i)=solutions(i).Exitflag;
end
EF=EF';

beep on;
beep
beep
beep
beep
beep
beep




% Want value=f(x) to be small value when run MultiStart  

