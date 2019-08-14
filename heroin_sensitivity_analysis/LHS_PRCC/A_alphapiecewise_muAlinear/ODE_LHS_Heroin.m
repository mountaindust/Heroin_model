function dydt=LV_ODE_LHS_Heroin(t,y,LHSmatrix,x)
%% PARAMETERS %%
Parameter_settings_LHS_Heroin;

m= LHSmatrix(x,1);
beta_A= LHSmatrix(x,2);
beta_P= LHSmatrix(x,3);
theta_1= LHSmatrix(x,4);
epsilon= LHSmatrix(x,5);
gamma= LHSmatrix(x,6);
sigma= LHSmatrix(x,7);
mu= LHSmatrix(x,8);
mu_H= LHSmatrix(x,9);
theta_2= LHSmatrix(x,10);
zeta= LHSmatrix(x,11);
theta_3= LHSmatrix(x,12);
nu= LHSmatrix(x,13);
omega= LHSmatrix(x,14);
b= LHSmatrix(x,15);
c=LHSmatrix(x,16);
d=LHSmatrix(x,17);
e=LHSmatrix(x,18);
P0=LHSmatrix(x,19);
A0=LHSmatrix(x,20);
H0=LHSmatrix(x,21);
R0=LHSmatrix(x,22);
S0=1-LHSmatrix(x,19)-LHSmatrix(x,20)-LHSmatrix(x,21)-LHSmatrix(x,22);

%% ODE model 
    
    dydt=[-a(t,m,b,c)*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+muA(t,d,e))*y(3)+(mu+mu_H)*y(4)
         a(t,m,b,c)*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2)
         gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-(mu+muA(t,d,e))*y(3)
         theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4)
         zeta*y(3)+nu*y(4)-sigma*y(5)*y(3)/(y(3)+y(4)+omega)-sigma*y(5)*y(4)/(y(3)+y(4)+omega)-mu*y(5)];
    
     
     
end


function alpha = a(t,m,b,c)
if  t<=3.25 
    alpha = m*t+b;
else 
    alpha = m*3.25+b-c*3.25+c*t;
end
end


function mu_A = muA(t,d,e)
    mu_A = d*t+e;
end
