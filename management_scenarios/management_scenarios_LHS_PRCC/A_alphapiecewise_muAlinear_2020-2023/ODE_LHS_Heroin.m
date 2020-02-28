function dydt=LV_ODE_LHS_Heroin(t,y,LHSmatrix,x)
%% PARAMETERS %%
%Parameter_settings_LHS_Heroin;

beta_A= LHSmatrix(x,1);
beta_P= LHSmatrix(x,2);
theta_1= LHSmatrix(x,3);
epsilon= LHSmatrix(x,4);
gamma= LHSmatrix(x,5);
sigma= LHSmatrix(x,6);
mu= LHSmatrix(x,7);
mu_H= LHSmatrix(x,8);
theta_2= LHSmatrix(x,9);
zeta= LHSmatrix(x,10);
theta_3= LHSmatrix(x,11);
nu= LHSmatrix(x,12);
omega= LHSmatrix(x,13);
g=LHSmatrix(x,14);
h=LHSmatrix(x,15);

%% ODE model 
    
    dydt=[-a(t,g)*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+muA(t,h))*y(3)+(mu+mu_H)*y(4)
         a(t,g)*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2)
         gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-(mu+muA(t,h))*y(3)
         theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4)
         zeta*y(3)+nu*y(4)-sigma*y(5)*y(3)/(y(3)+y(4)+omega)-sigma*y(5)*y(4)/(y(3)+y(4)+omega)-mu*y(5)];
    
     
     
end


function alpha = a(t,g)
    alpha = -0.00559565027929907*3.25+0.270110337915851+0.0269690987063522*3.25-0.0269690987063522*7+g*t;
end

%shift time by 7 because forced to start at 0
function mu_A = muA(t,h)
    mu_A = 0.000977482526657751*7+0.00883138792481281+h*t;
end
