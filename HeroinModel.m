
function dy = f(t,y,z)

% Parameters

alpha=z(1); 
 
beta=z(2); 
 
xi=z(3);
 
theta_1=z(4);
 
epsilon=z(5);
 
mu=z(6);  
 
mu_A=z(7);   
 
mu_H=z(8);
 
gamma=z(9);   
 
theta_2=z(10); 
 
sigma_A=z(11);
 
zeta=z(12);
 
theta_3=z(13);
 
sigma_H=z(14);
 
nu=z(15);


%ODE's 

dy(1) = -alpha*y(1)-beta*(1-xi)*y(1)*y(3)-beta*xi*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+mu_A)*y(3)+(mu+mu_H)*y(4);
dy(2) = alpha*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2);
dy(3) = gamma*y(2)+sigma_A*y(5)+beta*(1-xi)*y(1)*y(3)+beta*xi*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-mu*y(3)-mu_A*y(3);
dy(4) = theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma_H*y(5)-nu*y(4)-(mu+mu_H)*y(4);
dy(5) = zeta*y(3)+nu*y(4)-sigma_A*y(5)-sigma_H*y(5)-mu*y(5);
%ODE for new cases of prescription opioid use over time 
dy(6) = alpha*y(1);      
  dy=dy';  
 
end








           