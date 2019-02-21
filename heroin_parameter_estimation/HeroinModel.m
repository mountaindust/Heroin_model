% File name: HeroinModel.m (used to be in Heroin_model folder)

% Defining function with inputs: time, classes of people, and parameter vector z to estimate

function dy = HeroinModel(t,y,z)

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
 
theta_2=z(7); 
 
sigma=z(8);
 
zeta=z(9);
 
theta_3=z(10);
 
nu=z(11);

omega=0.0000000001;

P0=z(12);

A0=z(13);

H0=z(14);

R0=z(15); 

S0=1-z(12)-z(13)-z(14)-z(15); 


% ODEs 

dy(1) = -alpha*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+mu_A)*y(3)+(mu+mu_H)*y(4);
dy(2) = alpha*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2);
dy(3) = gamma*y(2)+(sigma*y(5)*y(3))/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-mu*y(3)-mu_A*y(3);
dy(4) = theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+(sigma*y(5)*y(4))/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4);
dy(5) = zeta*y(3)+nu*y(4)-(sigma*y(5)*y(3))/(y(3)+y(4)+omega)-(sigma*y(5)*y(4))/(y(3)+y(4)+omega)-mu*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
dy(6) = alpha*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
dy(7) = gamma*y(2)+(sigma*y(5)*y(3))/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
dy(8) = theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+(sigma*y(5)*y(4))/(y(3)+y(4)+omega);
% Transpose column vector solution into row vector 
  dy=dy';  
 
end








           
