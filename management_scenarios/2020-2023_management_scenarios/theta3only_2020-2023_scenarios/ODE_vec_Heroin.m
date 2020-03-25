function dydt=LV_ODE_vec_Heroin(t,y,vec_matrix,x)
%% PARAMETERS %%
Parameter_settings_vec_Heroin;
%global LHSmatrix

%m= vec_matrix(x,1);
beta_A= vec_matrix(x,1);
beta_P= vec_matrix(x,2);
theta_1= vec_matrix(x,3);
epsilon= vec_matrix(x,4);
gamma= vec_matrix(x,5);
sigma= vec_matrix(x,6);
mu= vec_matrix(x,7);
mu_H= vec_matrix(x,8);
theta_2= vec_matrix(x,9);
zeta= vec_matrix(x,10);
theta_3= vec_matrix(x,11);
nu= vec_matrix(x,12);
omega= vec_matrix(x,13);
g=vec_matrix(x,14);
h=vec_matrix(x,15);
%b= vec_matrix(x,15);
%c=vec_matrix(x,16);
%d=vec_matrix(x,17);
%e=vec_matrix(x,18);

P0=0.0585;
A0=0.0037;
H0=0.00597;
R0=0.00751;
S0=1-P0-A0-H0-R0;
%% Distemper model for shelter
%have included a 5th equation that keeps track of the total exposed, which
%will be the metric we are interested in

%  dydt = [b - beta_s*y(1)*y(3) - delta*y(1) - a_s*y(1)
%         beta_s*y(1)*y(3) + beta_v*y(4)*y(3) - a_e*y(2) - alpha*y(2)
%         alpha*y(2) - d*y(3)
%         delta*y(1) - beta_v*y(4)*y(3) - a_v*y(4)
%         beta_s*y(1)*y(3) + beta_v*y(4)*y(3)];

    
    
   dydt=[-a(t,g)*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+muA(t,h))*y(3)+(mu+mu_H)*y(4)
         a(t,g)*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2)
         gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-(mu+muA(t,h))*y(3)
         theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4)
         zeta*y(3)+nu*y(4)-sigma*y(5)*y(3)/(y(3)+y(4)+omega)-sigma*y(5)*y(4)/(y(3)+y(4)+omega)-mu*y(5)];
    
     
end


function alpha = a(t,g)
    alpha = -0.00559565027929907*3.25+0.270110337915851+0.0269690987063522*3.25-0.0269690987063522*7+g*t;
end


function mu_A = muA(t,h)
    mu_A = 0.000977482526657751*7+0.00883138792481281+h*t;
end


