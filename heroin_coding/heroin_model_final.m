%%%This code does the same things as HeroinModel_ODE45_Testing (put in parameters, get solution and plots out)%%%
% Defining function with inputs: time, classes of people, and parameter vector z to estimate

function dy = f(t,y)

% Parameters--values are realistic guesses for now until estimate parameters
    
alpha=0.15; 

beta_A= 0.00094; 
 
beta_P=0.00266;
 
theta_1=0.0003;
 
epsilon= 3.0;
 
mu=0.00868;  
 
mu_A=0.00775;   
 
mu_H=0.0271;
 
gamma=0.00744;   
 
theta_2=0.0003; 
 
sigma=1.5;
 
zeta=0.25;
 
theta_3=3;
 
nu=0.15;

omega=0.0000000001;




% ODEs 

dy(1) = -alpha*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+mu_A)*y(3)+(mu+mu_H)*y(4);
dy(2) = alpha*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2);
dy(3) = gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-mu*y(3)-mu_A*y(3);
dy(4) = theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4);
dy(5) = zeta*y(3)+nu*y(4)-sigma*y(5)*y(3)/(y(3)+y(4)+omega)-sigma*y(5)*y(4)/(y(3)+y(4)+omega)-mu*y(5);

  dy=dy';  
 
end








           
