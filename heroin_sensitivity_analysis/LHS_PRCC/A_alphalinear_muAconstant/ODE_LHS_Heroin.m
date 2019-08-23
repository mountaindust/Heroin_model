function dydt=LV_ODE_LHS_Heroin(t,y,LHSmatrix,x)
%% PARAMETERS %%
Parameter_settings_LHS_Heroin;
%global LHSmatrix

m= LHSmatrix(x,1);
beta_A= LHSmatrix(x,2);
beta_P= LHSmatrix(x,3);
theta_1= LHSmatrix(x,4);
epsilon= LHSmatrix(x,5);
gamma= LHSmatrix(x,6);
sigma= LHSmatrix(x,7);
mu= LHSmatrix(x,8);
mu_A= LHSmatrix(x,9);
mu_H= LHSmatrix(x,10);
theta_2= LHSmatrix(x,11);
zeta= LHSmatrix(x,12);
theta_3= LHSmatrix(x,13);
nu= LHSmatrix(x,14);
omega= LHSmatrix(x,15);
b= LHSmatrix(x,16);
P0=LHSmatrix(x,17);
A0=LHSmatrix(x,18);
H0=LHSmatrix(x,19);
R0=LHSmatrix(x,20);
S0=1-P0-A0-H0-R0;


     
%% Distemper model for shelter
%have included a 5th equation that keeps track of the total exposed, which
%will be the metric we are interested in

%  dydt = [b - beta_s*y(1)*y(3) - delta*y(1) - a_s*y(1)
%         beta_s*y(1)*y(3) + beta_v*y(4)*y(3) - a_e*y(2) - alpha*y(2)
%         alpha*y(2) - d*y(3)
%         delta*y(1) - beta_v*y(4)*y(3) - a_v*y(4)
%         beta_s*y(1)*y(3) + beta_v*y(4)*y(3)];
    
    
    
    dydt=[-(m*t+b)*y(1)-beta_A*y(1)*y(3)-beta_P*y(1)*y(2)-theta_1*y(1)*y(4)+epsilon*y(2)+mu*(y(2)+y(5))+(mu+mu_A)*y(3)+(mu+mu_H)*y(4)
         (m*t+b)*y(1)-epsilon*y(2)-gamma*y(2)-theta_2*y(2)*y(4)-mu*y(2)
         gamma*y(2)+sigma*y(5)*y(3)/(y(3)+y(4)+omega)+beta_A*y(1)*y(3)+beta_P*y(1)*y(2)-zeta*y(3)-theta_3*y(3)*y(4)-(mu+mu_A)*y(3)
         theta_1*y(1)*y(4)+theta_2*y(2)*y(4)+theta_3*y(3)*y(4)+sigma*y(5)*y(4)/(y(3)+y(4)+omega)-nu*y(4)-(mu+mu_H)*y(4)
         zeta*y(3)+nu*y(4)-sigma*y(5)*y(3)/(y(3)+y(4)+omega)-sigma*y(5)*y(4)/(y(3)+y(4)+omega)-mu*y(5)];
    
     
     
end

%LEAVING OUT %% FOR NOW from Christina's Code 

%%function val = gamma(epi1b,epi2b,hospitalizations)

% function for controlling the rate at which people seek medical attention
%In theory, as more incidents are noticed, media/health officials get
%involved and start informing the general public which in turn makes people
%more aware of an outbreak and makes them more likely to seek medical
%attention.

%I beleive this is a function of detected infections (ie, hospiutalized)
%rather than a function of dead indiviudals.
% maxval = 1;
% minval = .25;
%exponential
%val = .75*exp(-epi*hospitalizations);
%first sigmoidal try
%val = maxval/(1+((maxval-minval)/minval)*exp(-epi*hospitalizations));
%second sigmoidal try
%%val =(1-exp(-hospitalizations/epi1b))/(1+epi2b*exp(-hospitalizations/epi1b));

%%val = max(val,0);

%%end
    
   