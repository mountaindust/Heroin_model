% File name: HeroinModel.m (used to be in Heroin_model folder)

% Defining function with inputs: time, classes of people, and parameter vector z to estimate

function dy = f(t,y,z)

% Parameters

alpha=z(1); 

beta_A=z(2); 
 
beta_P=z(3);
 
theta_1=z(4);
 
epsilon=z(5);
 
mu=0.00868;  
 
mu_A=0.00775;   
 
mu_H=0.0271;
 
gamma=z(6);   
 
% Assume twice as likely for P individual to use heroin than an S individual 
theta_2=2*z(4); 
 
sigma=z(7);
 
zeta=z(8);
 
% Assume four times as likely for A individual to use heroin than an S individual 
theta_3=4*z(4);
 
nu=z(9);

omega=0.0000000001;

%Assume opioid-only project number (national level)
P0=0.05;

%Assume opioid-only project number (national level)
A0=0.0062;

%Assume 1/10th the number of heroin users compared to opioid addicts
H0=0.00062;

%Assume 1/10th of those who are addicted go into a stable recovery
R0=0.00062;


S0=1-0.05-0.0062-0.00062-0.00062; 


% ODEs 

dy(1) = -alpha*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+mu_A)*y(3)+(mu+mu_H)*y(4);
dy(2) = alpha*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2);
dy(3) = gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-mu*y(3)-mu_A*y(3);
dy(4) = theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4);
dy(5) = zeta*y(3)+nu*y(4)-sigma*y(5)*y(3)/(y(3)+y(4)+omega)-sigma*y(5)*y(4)/(y(3)+y(4)+omega)-mu*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
dy(6) = alpha*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
dy(7) = gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
dy(8) = theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega);

  dy=dy';  
 
end








           