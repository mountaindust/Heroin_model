
clear all;
close all;

%% PRCC

%% Sample size N
 
%Total # of parameters values to test from each parameter interval (i.e. number of uniform intervals)
nsample = 1000; 

%% LHS MATRIX  %%

%Using baseline parameter settings 
Parameter_settings_LHS_Heroin; 

%xmin=0;
%xmax=xmin;


%The following uses LHS_Call_Distemper to set the vector of each parameter to be
%tested. LHS_Call_Distemper follows a method by Stein (cited in the LHS_Call file)
%Each parameter value will be a vector that will be arranged into a matrix
%(the latin hypercube matrix) shown below


%from LHS_Call_Distemper:
%s=LHS_Call_Distemper(xmin,xmean,xmax,xsd,nsample,distrib,logscale)
%none of the xmean or xsd values are nonzero because for the uniform
%distribution option you only need a max and min
% 

%Choosing the points from the +/-50% baseline parameter intervals using a uniform
%distribution (min, mean, max, std dev, number samples, distribution used)
m_LHS=LHS_Call_Heroin(m-(m/2),0,m+(m/2),0,nsample,'unif');
beta_A_LHS=LHS_Call_Heroin(beta_A-(beta_A/2),0,beta_A+(beta_A/2),0,nsample,'unif');
beta_P_LHS=LHS_Call_Heroin(beta_P-(beta_P/2),0,beta_P+(beta_P/2),0,nsample,'unif');
theta_1_LHS=LHS_Call_Heroin(theta_1-(theta_1/2),0,theta_1+(theta_1/2),0,nsample,'unif');
epsilon_LHS=LHS_Call_Heroin(epsilon-(epsilon/2),0,epsilon+(epsilon/2),0,nsample,'unif');
gamma_LHS=LHS_Call_Heroin(gamma-(gamma/2),0,gamma+(gamma/2),0,nsample,'unif');
sigma_LHS=LHS_Call_Heroin(sigma-(sigma/2),0,sigma+(sigma/2),0,nsample,'unif');
mu_LHS=LHS_Call_Heroin(mu-(mu/2),0,mu+(mu/2),0,nsample,'unif');
mu_A_LHS=LHS_Call_Heroin(mu_A-(mu_A/2),0,mu_A+(mu_A/2),0,nsample,'unif');
mu_H_LHS=LHS_Call_Heroin(mu_H-(mu_H/2),0,mu_H+(mu_H/2),0,nsample,'unif');
theta_2_LHS=LHS_Call_Heroin(theta_2-(theta_2/2),0,theta_2+(theta_2/2),0,nsample,'unif');
zeta_LHS=LHS_Call_Heroin(zeta-(zeta/2),0,zeta+(zeta/2),0,nsample,'unif');
theta_3_LHS=LHS_Call_Heroin(theta_3-(theta_3/2),0,theta_3+(theta_3/2),0,nsample,'unif');
nu_LHS=LHS_Call_Heroin(nu-(nu/2),0,nu+(nu/2),0,nsample,'unif');
omega_LHS=LHS_Call_Heroin(omega-(omega/2),0,omega+(omega/2),0,nsample,'unif');
b_LHS=LHS_Call_Heroin(b-(b/2),0,b+(b/2),0,nsample,'unif');


 
 

%% LHS MATRIX and PARAMETER LABELS
 %storing vectors from above in matrix form 
  LHSmatrix  = [m_LHS,beta_A_LHS,beta_P_LHS,theta_1_LHS,epsilon_LHS,gamma_LHS,sigma_LHS,mu_LHS,mu_A_LHS,mu_H_LHS,theta_2_LHS,zeta_LHS,theta_3_LHS,nu_LHS,omega_LHS,b_LHS];

 
for x=1:nsample %Run solution x times choosing different values
    f=@ODE_LHS_Heroin;
%     x;
%     
%      LHSmatrix(x,:);
     
     [t,y] = ode15s(@(t,y)f(t,y,LHSmatrix,x),tspan,y0,[]); 
     %[t,y]=ode15s(@ODE_LHS_Heroin,tspan,y0,[],LHSmatrix);
     W = [t y]; % [time y]
    
   
    %% Save the outputs at ALL time points [tspan]
     %first column of W is are the values of t 
     S_lhs(:,x)=W(:,2);
     P_lhs(:,x)=W(:,3);
     A_lhs(:,x)=W(:,4);
     H_lhs(:,x)=W(:,5);
     R_lhs(:,x)=W(:,6);
     
    %% Save only the outputs at the time points of interest [time_points]:
    %%More efficient    

%     N_lhs(:,x)=A(time_points+1,2);
%     S_lhs(:,x)=A(time_points+1,3);
%     I_lhs(:,x)=A(time_points+1,4);

end
 
%% Save the workspace
 
save LV_Model_LHS_Heroin.mat;

%% CALCULATE PRCC


  alpha = 1e-3;
  %this uses total exposures as the metric. Change second input argument if want
  %to test something else
 [prcc sign sign_label]=PRCC_Heroin(LHSmatrix,H_lhs,time_points,PRCC_var,alpha); %PRCC_var and time_points set in parameter file


%% Scatter plots

  PRCC_PLOT_Heroin(LHSmatrix,H_lhs,time_points,PRCC_var,'Total Heroin Addicts');
  
  
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 %% Monotonicity plots
 

%fixing all but running LHS for one at a time

%m
       
m_LHS1=LHS_Call_Heroin(m-(m/2),0,m+(m/2),0,nsample,'unif');
beta_A_LHS1=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS1=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS1=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS1=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS1=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS1=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS1=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS1=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS1=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS1=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS1=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS1=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS1=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS1=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS1=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%beta_A

m_LHS2=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS2=LHS_Call_Heroin(beta_A-(beta_A/2),0,beta_A+(beta_A/2),0,nsample,'unif');
beta_P_LHS2=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS2=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS2=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS2=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS2=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS2=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS2=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS2=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS2=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS2=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS2=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS2=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS2=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS2=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');


%beta_P

m_LHS3=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS3=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS3=LHS_Call_Heroin(beta_P-(beta_P/2),0,beta_P+(beta_P/2),0,nsample,'unif');
theta_1_LHS3=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS3=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS3=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS3=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS3=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS3=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS3=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS3=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS3=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS3=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS3=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS3=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS3=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');


%theta_1

m_LHS4=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS4=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS4=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS4=LHS_Call_Heroin(theta_1-(theta_1/2),0,theta_1+(theta_1/2),0,nsample,'unif');
epsilon_LHS4=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS4=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS4=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS4=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS4=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS4=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS4=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS4=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS4=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS4=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS4=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS4=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');





%epsilon

m_LHS5=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS5=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS5=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS5=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS5=LHS_Call_Heroin(epsilon-(epsilon/2),0,epsilon+(epsilon/2),0,nsample,'unif');
gamma_LHS5=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS5=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS5=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS5=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS5=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS5=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS5=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS5=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS5=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS5=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS5=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');






%gamma

m_LHS6=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS6=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS6=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS6=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS6=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS6=LHS_Call_Heroin(gamma-(gamma/2),0,gamma+(gamma/2),0,nsample,'unif');
sigma_LHS6=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS6=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS6=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS6=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS6=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS6=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS6=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS6=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS6=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS6=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');



%sigma

m_LHS7=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS7=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS7=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS7=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS7=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS7=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS7=LHS_Call_Heroin(sigma-(sigma/2),0,sigma+(sigma/2),0,nsample,'unif');
mu_LHS7=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS7=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS7=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS7=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS7=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS7=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS7=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS7=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS7=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');



%mu

m_LHS8=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS8=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS8=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS8=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS8=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS8=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS8=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS8=LHS_Call_Heroin(mu-(mu/2),0,mu+(mu/2),0,nsample,'unif');
mu_A_LHS8=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS8=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS8=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS8=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS8=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS8=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS8=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS8=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%mu_A

m_LHS9=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS9=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS9=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS9=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS9=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS9=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS9=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS9=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS9=LHS_Call_Heroin(mu_A-(mu_A/2),0,mu_A+(mu_A/2),0,nsample,'unif');
mu_H_LHS9=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS9=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS9=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS9=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS9=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS9=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS9=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%mu_H

m_LHS10=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS10=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS10=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS10=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS10=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS10=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS10=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS10=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS10=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS10=LHS_Call_Heroin(mu_H-(mu_H/2),0,mu_H+(mu_H/2),0,nsample,'unif');
theta_2_LHS10=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS10=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS10=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS10=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS10=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS10=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%theta_2
m_LHS11=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS11=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS11=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS11=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS11=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS11=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS11=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS11=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS11=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS11=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS11=LHS_Call_Heroin(theta_2-(theta_2/2),0,theta_2+(theta_2/2),0,nsample,'unif');
zeta_LHS11=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS11=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS11=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS11=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS11=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');


%zeta

m_LHS12=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS12=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS12=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS12=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS12=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS12=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS12=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS12=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS12=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS12=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS12=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS12=LHS_Call_Heroin(zeta-(zeta/2),0,zeta+(zeta/2),0,nsample,'unif');
theta_3_LHS12=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS12=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS12=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS12=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');


%theta_3

m_LHS13=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS13=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS13=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS13=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS13=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS13=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS13=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS13=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS13=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS13=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS13=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS13=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS13=LHS_Call_Heroin(theta_3-(theta_3/2),0,theta_3+(theta_3/2),0,nsample,'unif');
nu_LHS13=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS13=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS13=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%nu 

m_LHS14=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS14=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS14=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS14=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS14=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS14=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS14=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS14=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS14=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS14=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS14=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS14=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS14=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS14=LHS_Call_Heroin(nu-(nu/2),0,nu+(nu/2),0,nsample,'unif');
omega_LHS14=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS14=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%omega

m_LHS15=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS15=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS15=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS15=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS15=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS15=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS15=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS15=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS15=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS15=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS15=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS15=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS15=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS15=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS15=LHS_Call_Heroin(omega-(omega/2),0,omega+(omega/2),0,nsample,'unif');
b_LHS15=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');

%b

m_LHS16=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS16=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS16=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS16=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS16=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS16=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS16=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS16=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS16=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS16=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS16=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS16=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS16=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS16=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS16=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS16=LHS_Call_Heroin(b-(b/2),0,b+(b/2),0,nsample,'unif');




%% LHS MATRIX and PARAMETER LABELS

  
LHSmatrix1 =[m_LHS1 beta_A_LHS1 beta_P_LHS1 theta_1_LHS1 epsilon_LHS1 gamma_LHS1 sigma_LHS1 mu_LHS1 mu_A_LHS1 mu_H_LHS1 theta_2_LHS1 zeta_LHS1 theta_3_LHS1   nu_LHS1 omega_LHS1 b_LHS1 ];
LHSmatrix2 =[m_LHS2 beta_A_LHS2 beta_P_LHS2 theta_1_LHS2 epsilon_LHS2 gamma_LHS2 sigma_LHS2 mu_LHS2 mu_A_LHS2 mu_H_LHS2 theta_2_LHS2 zeta_LHS2 theta_3_LHS2 nu_LHS2 omega_LHS2 b_LHS2 ];
LHSmatrix3 =[m_LHS3 beta_A_LHS3 beta_P_LHS3 theta_1_LHS3 epsilon_LHS3 gamma_LHS3 sigma_LHS3 mu_LHS3 mu_A_LHS3 mu_H_LHS3 theta_2_LHS3 zeta_LHS3 theta_3_LHS3 nu_LHS3 omega_LHS3 b_LHS3 ];
LHSmatrix4 =[m_LHS4 beta_A_LHS4 beta_P_LHS4 theta_1_LHS4 epsilon_LHS4 gamma_LHS4 sigma_LHS4 mu_LHS4 mu_A_LHS4 mu_H_LHS4 theta_2_LHS4 zeta_LHS4 theta_3_LHS4 nu_LHS4 omega_LHS4 b_LHS4 ];
LHSmatrix5 =[m_LHS5 beta_A_LHS5 beta_P_LHS5 theta_1_LHS5 epsilon_LHS5 gamma_LHS5 sigma_LHS5 mu_LHS5 mu_A_LHS5 mu_H_LHS5 theta_2_LHS5 zeta_LHS5 theta_3_LHS5 nu_LHS5 omega_LHS5 b_LHS5 ]; 
LHSmatrix6 =[m_LHS6 beta_A_LHS6 beta_P_LHS6 theta_1_LHS6 epsilon_LHS6 gamma_LHS6 sigma_LHS6 mu_LHS6 mu_A_LHS6 mu_H_LHS6 theta_2_LHS6 zeta_LHS6 theta_3_LHS6 nu_LHS6 omega_LHS6 b_LHS6 ];
LHSmatrix7 =[m_LHS7 beta_A_LHS7 beta_P_LHS7 theta_1_LHS7 epsilon_LHS7 gamma_LHS7 sigma_LHS7 mu_LHS7 mu_A_LHS7 mu_H_LHS7 theta_2_LHS7 zeta_LHS7 theta_3_LHS7 nu_LHS7 omega_LHS7 b_LHS7 ];
LHSmatrix8 =[m_LHS8 beta_A_LHS8 beta_P_LHS8 theta_1_LHS8 epsilon_LHS8 gamma_LHS8 sigma_LHS8 mu_LHS8 mu_A_LHS8 mu_H_LHS8 theta_2_LHS8 zeta_LHS8 theta_3_LHS8 nu_LHS8 omega_LHS8 b_LHS8 ];
LHSmatrix9 =[m_LHS9 beta_A_LHS9 beta_P_LHS9 theta_1_LHS9 epsilon_LHS9 gamma_LHS9 sigma_LHS9 mu_LHS9 mu_A_LHS9 mu_H_LHS9 theta_2_LHS9 zeta_LHS9 theta_3_LHS9 nu_LHS9 omega_LHS9 b_LHS9 ];
LHSmatrix10 =[m_LHS10 beta_A_LHS10 beta_P_LHS10 theta_1_LHS10 epsilon_LHS10 gamma_LHS10 sigma_LHS10 mu_LHS10 mu_A_LHS10 mu_H_LHS10 theta_2_LHS10 zeta_LHS10 theta_3_LHS10 nu_LHS10 omega_LHS10 b_LHS10 ];
LHSmatrix11 =[m_LHS11 beta_A_LHS11 beta_P_LHS11 theta_1_LHS11 epsilon_LHS11 gamma_LHS11 sigma_LHS11 mu_LHS11 mu_A_LHS11 mu_H_LHS11 theta_2_LHS11 zeta_LHS11 theta_3_LHS11 nu_LHS11 omega_LHS11 b_LHS11 ];
LHSmatrix12 =[m_LHS12 beta_A_LHS12 beta_P_LHS12 theta_1_LHS12 epsilon_LHS12 gamma_LHS12 sigma_LHS12 mu_LHS12 mu_A_LHS12 mu_H_LHS12 theta_2_LHS12 zeta_LHS12 theta_3_LHS12 nu_LHS12 omega_LHS12 b_LHS12 ];
LHSmatrix13 =[m_LHS13 beta_A_LHS13 beta_P_LHS13 theta_1_LHS13 epsilon_LHS13 gamma_LHS13 sigma_LHS13 mu_LHS13 mu_A_LHS13 mu_H_LHS13 theta_2_LHS13 zeta_LHS13 theta_3_LHS13 nu_LHS13 omega_LHS13 b_LHS13 ];
LHSmatrix14 =[m_LHS14 beta_A_LHS14 beta_P_LHS14 theta_1_LHS14 epsilon_LHS14 gamma_LHS14 sigma_LHS14 mu_LHS14 mu_A_LHS14 mu_H_LHS14 theta_2_LHS14 zeta_LHS14 theta_3_LHS14 nu_LHS14 omega_LHS14 b_LHS14 ];
LHSmatrix15 =[m_LHS15 beta_A_LHS15 beta_P_LHS15 theta_1_LHS15 epsilon_LHS15 gamma_LHS15 sigma_LHS15 mu_LHS15 mu_A_LHS15 mu_H_LHS15 theta_2_LHS15 zeta_LHS15 theta_3_LHS15 nu_LHS15 omega_LHS15 b_LHS15 ]; 
LHSmatrix16 =[m_LHS16 beta_A_LHS16 beta_P_LHS16 theta_1_LHS16 epsilon_LHS16 gamma_LHS16 sigma_LHS16 mu_LHS16 mu_A_LHS16 mu_H_LHS16 theta_2_LHS16 zeta_LHS16 theta_3_LHS16 nu_LHS16 omega_LHS16 b_LHS16 ];

  
for x=1:nsample %Run solution x times choosing different values, represents each row of the matrix that's going to go through the ODE solver to produce the plots 
    f=@ODE_LHS_Heroin;
    x;
     
     %values to test
    LHSmatrix1(x,:);
    LHSmatrix2(x,:);
    LHSmatrix3(x,:);
    LHSmatrix4(x,:);
    LHSmatrix5(x,:);
    LHSmatrix6(x,:);
    LHSmatrix7(x,:);
    LHSmatrix8(x,:);
    LHSmatrix9(x,:);
    LHSmatrix10(x,:);
    LHSmatrix11(x,:);
    LHSmatrix12(x,:);
    LHSmatrix13(x,:);
    LHSmatrix14(x,:);
    LHSmatrix15(x,:);
    LHSmatrix16(x,:);
 
    
    %run each with only 1 parameter varying in each (each ODE run nsample
    %times, 1 time for each row in LHSmatrix# to produce output for each parameter in a
    %certain range)
    [t,y1] = ode15s(@(t,y1)f(t,y1,LHSmatrix1,x),tspan,y0,[]); 
    [t,y2] = ode15s(@(t,y2)f(t,y2,LHSmatrix2,x),tspan,y0,[]); 
    [t,y3] = ode15s(@(t,y3)f(t,y3,LHSmatrix3,x),tspan,y0,[]);
    [t,y4] = ode15s(@(t,y4)f(t,y4,LHSmatrix4,x),tspan,y0,[]); 
    [t,y5] = ode15s(@(t,y5)f(t,y5,LHSmatrix5,x),tspan,y0,[]); 
    [t,y6] = ode15s(@(t,y6)f(t,y6,LHSmatrix6,x),tspan,y0,[]);
    [t,y7] = ode15s(@(t,y7)f(t,y7,LHSmatrix7,x),tspan,y0,[]); 
    [t,y8] = ode15s(@(t,y8)f(t,y8,LHSmatrix8,x),tspan,y0,[]);
    [t,y9] = ode15s(@(t,y9)f(t,y9,LHSmatrix9,x),tspan,y0,[]);
    [t,y10] = ode15s(@(t,y10)f(t,y10,LHSmatrix10,x),tspan,y0,[]);
    [t,y11] = ode15s(@(t,y11)f(t,y11,LHSmatrix11,x),tspan,y0,[]);
    [t,y12] = ode15s(@(t,y12)f(t,y12,LHSmatrix12,x),tspan,y0,[]);
    [t,y13] = ode15s(@(t,y13)f(t,y13,LHSmatrix13,x),tspan,y0,[]);
    [t,y14] = ode15s(@(t,y14)f(t,y14,LHSmatrix14,x),tspan,y0,[]);
    [t,y15] = ode15s(@(t,y15)f(t,y15,LHSmatrix15,x),tspan,y0,[]);
    [t,y16] = ode15s(@(t,y16)f(t,y16,LHSmatrix16,x),tspan,y0,[]);
   
 
    
     %store results of each [t,y1] = ode15s(@(t,y1)f(t,y1,LHSmatrix1,x),tspan,y0,[]);
     
     W1 = [t y1]; % [time y]
     W2 = [t y2];
     W3 = [t y3];
     W4 = [t y4]; 
     W5 = [t y5];
     W6 = [t y6];
     W7 = [t y7];
     W8 = [t y8];
     W9 = [t y9];
     W10 = [t y10]; 
     W11 = [t y11]; 
     W12 = [t y12]; 
     W13 = [t y13]; 
     W14 = [t y14]; 
     W15 = [t y15]; 
     W16 = [t y16]; 
   
     
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %store final value of each class (time_points+1 is final time since first
  %column is IC's), start with 2nd column because 1st just counts classes
     S_lhs1(:,x)=W1(time_points+1,2);
     P_lhs1(:,x)=W1(time_points+1,3);
     A_lhs1(:,x)=W1(time_points+1,4);
     H_lhs1(:,x)=W1(time_points+1,5);
     R_lhs1(:,x)=W1(time_points+1,6);

     
     S_lhs2(:,x)=W2(time_points+1,2);
     P_lhs2(:,x)=W2(time_points+1,3);
     A_lhs2(:,x)=W2(time_points+1,4);
     H_lhs2(:,x)=W2(time_points+1,5);
     R_lhs2(:,x)=W2(time_points+1,6);
   
     
     S_lhs3(:,x)=W3(time_points+1,2);
     P_lhs3(:,x)=W3(time_points+1,3);
     A_lhs3(:,x)=W3(time_points+1,4);
     H_lhs3(:,x)=W3(time_points+1,5);
     R_lhs3(:,x)=W3(time_points+1,6);
     
     
     S_lhs4(:,x)=W4(time_points+1,2);
     P_lhs4(:,x)=W4(time_points+1,3);
     A_lhs4(:,x)=W4(time_points+1,4);
     H_lhs4(:,x)=W4(time_points+1,5);
     R_lhs4(:,x)=W4(time_points+1,6);
    
     S_lhs5(:,x)=W5(time_points+1,2);
     P_lhs5(:,x)=W5(time_points+1,3);
     A_lhs5(:,x)=W5(time_points+1,4);
     H_lhs5(:,x)=W5(time_points+1,5);
     R_lhs5(:,x)=W5(time_points+1,6);
  
     
     S_lhs6(:,x)=W6(time_points+1,2);
     P_lhs6(:,x)=W6(time_points+1,3);
     A_lhs6(:,x)=W6(time_points+1,4);
     H_lhs6(:,x)=W6(time_points+1,5);
     R_lhs6(:,x)=W6(time_points+1,6);
  
     
     S_lhs7(:,x)=W7(time_points+1,2);
     P_lhs7(:,x)=W7(time_points+1,3);
     A_lhs7(:,x)=W7(time_points+1,4);
     H_lhs7(:,x)=W7(time_points+1,5);
     R_lhs7(:,x)=W7(time_points+1,6);
    
     
     S_lhs8(:,x)=W8(time_points+1,2);
     P_lhs8(:,x)=W8(time_points+1,3);
     A_lhs8(:,x)=W8(time_points+1,4);
     H_lhs8(:,x)=W8(time_points+1,5);
     R_lhs8(:,x)=W8(time_points+1,6);
   
     
     S_lhs9(:,x)=W9(time_points+1,2);
     P_lhs9(:,x)=W9(time_points+1,3);
     A_lhs9(:,x)=W9(time_points+1,4);
     H_lhs9(:,x)=W9(time_points+1,5);
     R_lhs9(:,x)=W9(time_points+1,6);
  
     
     S_lhs10(:,x)=W10(time_points+1,2);
     P_lhs10(:,x)=W10(time_points+1,3);
     A_lhs10(:,x)=W10(time_points+1,4);
     H_lhs10(:,x)=W10(time_points+1,5);
     R_lhs10(:,x)=W10(time_points+1,6);
   
     S_lhs11(:,x)=W11(time_points+1,2);
     P_lhs11(:,x)=W11(time_points+1,3);
     A_lhs11(:,x)=W11(time_points+1,4);
     H_lhs11(:,x)=W11(time_points+1,5);
     R_lhs11(:,x)=W11(time_points+1,6);
    
     
     S_lhs12(:,x)=W12(time_points+1,2);
     P_lhs12(:,x)=W12(time_points+1,3);
     A_lhs12(:,x)=W12(time_points+1,4);
     H_lhs12(:,x)=W12(time_points+1,5);
     R_lhs12(:,x)=W12(time_points+1,6);

     
     S_lhs13(:,x)=W13(time_points+1,2);
     P_lhs13(:,x)=W13(time_points+1,3);
     A_lhs13(:,x)=W13(time_points+1,4);
     H_lhs13(:,x)=W13(time_points+1,5);
     R_lhs13(:,x)=W13(time_points+1,6);
    
     
     S_lhs14(:,x)=W14(time_points+1,2);
     P_lhs14(:,x)=W14(time_points+1,3);
     A_lhs14(:,x)=W14(time_points+1,4);
     H_lhs14(:,x)=W14(time_points+1,5);
     R_lhs14(:,x)=W14(time_points+1,6);
 
     
     S_lhs15(:,x)=W15(time_points+1,2);
     P_lhs15(:,x)=W15(time_points+1,3);
     A_lhs15(:,x)=W15(time_points+1,4);
     H_lhs15(:,x)=W15(time_points+1,5);
     R_lhs15(:,x)=W15(time_points+1,6);
 
     
     
     S_lhs16(:,x)=W16(time_points+1,2);
     P_lhs16(:,x)=W16(time_points+1,3);
     A_lhs16(:,x)=W16(time_points+1,4);
     H_lhs16(:,x)=W16(time_points+1,5);
     R_lhs16(:,x)=W16(time_points+1,6);
    
    
end
%% Save the workspace
  
 save LV_Model_LHS_Heroin.mat;
 
%  %% CALCULATE PRCC 

%only stored last time point so 1:length(time_points)=1 makes sense
[prcc1 sign1 sign_label1] = PRCC_Heroin(LHSmatrix1,H_lhs1,1:length(time_points),PRCC_var,alpha);
[prcc2 sign2 sign_label2] = PRCC_Heroin(LHSmatrix2,H_lhs2,1:length(time_points),PRCC_var,alpha);
[prcc3 sign3 sign_label3] = PRCC_Heroin(LHSmatrix3,H_lhs3,1:length(time_points),PRCC_var,alpha);
[prcc4 sign4 sign_label4] = PRCC_Heroin(LHSmatrix4,H_lhs4,1:length(time_points),PRCC_var,alpha);
[prcc5 sign5 sign_label5] = PRCC_Heroin(LHSmatrix5,H_lhs5,1:length(time_points),PRCC_var,alpha);
[prcc6 sign6 sign_label6] = PRCC_Heroin(LHSmatrix6,H_lhs6,1:length(time_points),PRCC_var,alpha);
[prcc7 sign7 sign_label7] = PRCC_Heroin(LHSmatrix7,H_lhs7,1:length(time_points),PRCC_var,alpha);
[prcc8 sign8 sign_label8] = PRCC_Heroin(LHSmatrix8,H_lhs8,1:length(time_points),PRCC_var,alpha);
[prcc9 sign9 sign_label9] = PRCC_Heroin(LHSmatrix9,H_lhs9,1:length(time_points),PRCC_var,alpha);
[prcc10 sign10 sign_label10] = PRCC_Heroin(LHSmatrix10,H_lhs10,1:length(time_points),PRCC_var,alpha);
[prcc11 sign11 sign_label11] = PRCC_Heroin(LHSmatrix11,H_lhs11,1:length(time_points),PRCC_var,alpha);
[prcc12 sign12 sign_label12] = PRCC_Heroin(LHSmatrix12,H_lhs12,1:length(time_points),PRCC_var,alpha);
[prcc13 sign13 sign_label13] = PRCC_Heroin(LHSmatrix13,H_lhs13,1:length(time_points),PRCC_var,alpha);
[prcc14 sign14 sign_label14] = PRCC_Heroin(LHSmatrix14,H_lhs14,1:length(time_points),PRCC_var,alpha);
[prcc15 sign15 sign_label15] = PRCC_Heroin(LHSmatrix15,H_lhs15,1:length(time_points),PRCC_var,alpha);
[prcc16 sign16 sign_label16] = PRCC_Heroin(LHSmatrix16,H_lhs16,1:length(time_points),PRCC_var,alpha);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Monotonicity curves for opioid addicted individuals at a particular time point (last time step)
figure(1);
 subplot(331)
 %take first column of LHSmatrix1 because those are the values m is varying over to produce the output in A_lhs1
 plot(LHSmatrix1(:,1),A_lhs1,'o') 
 xlabel('m')
 ylabel('Opioid addicted')
 subplot(332)
 plot(LHSmatrix2(:,2),A_lhs2,'o')
 xlabel('beta_A')
 ylabel('Opioid addicted')
 subplot(333)
 plot(LHSmatrix3(:,3),A_lhs3,'o')
 xlabel('beta_P')
 ylabel('Opioid addicted')
subplot(334)
plot(LHSmatrix4(:,4),A_lhs4,'o')
xlabel('theta_1')
ylabel('Opioid addicted')
subplot(335)
plot(LHSmatrix5(:,5),A_lhs5,'o')
xlabel('epsilon')
ylabel('Opioid addicted')
subplot(336)
plot(LHSmatrix6(:,6),A_lhs6,'o')
xlabel('gamma')
ylabel('Opioid addicted')
subplot(337)
plot(LHSmatrix7(:,7),A_lhs7,'o')
xlabel('sigma')
ylabel('Opioid addicted')
subplot(338)
plot(LHSmatrix8(:,8),A_lhs8,'o')
xlabel('mu')
ylabel('Opioid addicted')
subplot(339)
plot(LHSmatrix9(:,9),A_lhs9,'o')
xlabel('mu_A')
ylabel('Opioid addicted') 

figure(2);
subplot(331)
 plot(LHSmatrix10(:,10),A_lhs10,'o')
 xlabel('mu_H')
 ylabel('Opioid addicted')
 subplot(332)
 plot(LHSmatrix11(:,11),A_lhs11,'o')
 xlabel('theta_2')
 ylabel('Opioid addicted')
 subplot(333)
 plot(LHSmatrix12(:,12),A_lhs12,'o')
 xlabel('zeta')
 ylabel('Opioid addicted')
 subplot(334)
 plot(LHSmatrix13(:,13),A_lhs13,'o')
 xlabel('theta_3')
 ylabel('Opioid addicted')
 subplot(335)
 plot(LHSmatrix14(:,14),A_lhs14,'o')
 xlabel('nu')
 ylabel('Opioid addicted')
 subplot(337)
 plot(LHSmatrix15(:,15),A_lhs15,'o')
 xlabel('omega')
 ylabel('Opioid addicted')
 subplot(338)
 plot(LHSmatrix16(:,16),A_lhs16,'o')
 xlabel('b')
 ylabel('Opioid addicted')

 
 %Monotonicity curves for heroin addicted individuals at last time step 
 
 figure(3);
 subplot(331)
 plot(LHSmatrix1(:,1),H_lhs1,'o')
 xlabel('m')
 ylabel('Heroin addicted')
 subplot(332)
 plot(LHSmatrix2(:,2),H_lhs2,'o')
 xlabel('beta_A')
 ylabel('Heroin addicted')
 subplot(333)
 plot(LHSmatrix3(:,3),H_lhs3,'o')
 xlabel('beta_P')
 ylabel('Heroin addicted')
subplot(334)
plot(LHSmatrix4(:,4),H_lhs4,'o')
xlabel('theta_1')
ylabel('Heroin addicted')
subplot(335)
plot(LHSmatrix5(:,5),H_lhs5,'o')
xlabel('epsilon')
ylabel('Heroin addicted')
subplot(336)
plot(LHSmatrix6(:,6),H_lhs6,'o')
xlabel('gamma')
ylabel('Heroin addicted')
subplot(337)
plot(LHSmatrix7(:,7),H_lhs7,'o')
xlabel('sigma')
ylabel('Heroin addicted')
subplot(338)
plot(LHSmatrix8(:,8),H_lhs8,'o')
xlabel('mu')
ylabel('Heroin addicted')
subplot(339)
plot(LHSmatrix9(:,9),H_lhs9,'o')
xlabel('mu_A')
ylabel('Heroin addicted') 

figure(4);
subplot(331)
 plot(LHSmatrix10(:,10),H_lhs10,'o')
 xlabel('mu_H')
 ylabel('Heroin addicted')
 subplot(332)
 plot(LHSmatrix11(:,11),H_lhs11,'o')
 xlabel('theta_2')
 ylabel('Heroin addicted')
 subplot(333)
 plot(LHSmatrix12(:,12),H_lhs12,'o')
 xlabel('zeta')
 ylabel('Heroin addicted')
 subplot(334)
 plot(LHSmatrix13(:,13),H_lhs13,'o')
 xlabel('theta_3')
 ylabel('Heroin addicted')
 subplot(335)
 plot(LHSmatrix14(:,14),H_lhs14,'o')
 xlabel('nu')
 ylabel('Heroin addicted')
 subplot(337)
 plot(LHSmatrix15(:,15),H_lhs15,'o')
 %set(gca,'yticklabel',num2str(get(gca,'ytick')','%.9f'))
 xlabel('omega')
 ylabel('Heroin addicted') 
 subplot(338)
 plot(LHSmatrix16(:,16),H_lhs16,'o')
 xlabel('b')
 ylabel('Heroin addicted')


 
 figure(5);
 subplot(331)
 plot(LHSmatrix1(:,1),S_lhs1,'o')
 xlabel('m')
 ylabel('Susceptible Individuals')
 subplot(332)
 plot(LHSmatrix2(:,2),S_lhs2,'o')
 xlabel('beta_A')
 ylabel('Susceptible Individuals')
 subplot(333)
 plot(LHSmatrix3(:,3),S_lhs3,'o')
 xlabel('beta_P')
 ylabel('Susceptible Individuals')
subplot(334)
plot(LHSmatrix4(:,4),S_lhs4,'o')
xlabel('theta_1')
ylabel('Susceptible Individuals')
subplot(335)
plot(LHSmatrix5(:,5),S_lhs5,'o')
xlabel('epsilon')
ylabel('Susceptible Individuals')
subplot(336)
plot(LHSmatrix6(:,6),S_lhs6,'o')
xlabel('gamma')
ylabel('Susceptible Individuals')
subplot(337)
plot(LHSmatrix7(:,7),S_lhs7,'o')
xlabel('sigma')
ylabel('Susceptible Individuals')
subplot(338)
plot(LHSmatrix8(:,8),S_lhs8,'o')
xlabel('mu')
ylabel('Susceptible Individuals')
subplot(339)
plot(LHSmatrix9(:,9),S_lhs9,'o')
xlabel('mu_A')
ylabel('Susceptible Individuals')

figure(6);
subplot(331)
 plot(LHSmatrix10(:,10),S_lhs10,'o')
 xlabel('mu_H')
 ylabel('Susceptible Individuals')
 subplot(332)
 plot(LHSmatrix11(:,11),S_lhs11,'o')
 xlabel('theta_2')
 ylabel('Susceptible Individuals')
 subplot(333)
 plot(LHSmatrix12(:,12),S_lhs12,'o')
 xlabel('zeta')
 ylabel('Susceptible Individuals')
 subplot(334)
 plot(LHSmatrix13(:,13),S_lhs13,'o')
 xlabel('theta_3')
 ylabel('Susceptible Individuals')
 subplot(335)
 plot(LHSmatrix14(:,14),S_lhs14,'o')
 xlabel('nu')
 ylabel('Susceptible Individuals')
 subplot(337)
 plot(LHSmatrix15(:,15),S_lhs15,'o')
 %set(gca,'yticklabel',num2str(get(gca,'ytick')','%.9f'))
 xlabel('omega')
 ylabel('Susceptible Individuals')
 subplot(338)
 plot(LHSmatrix16(:,16),S_lhs16,'o')
 xlabel('b')
 ylabel('Susceptible Individuals')

 
 
 
 figure(7);
 subplot(331)
 plot(LHSmatrix1(:,1),P_lhs1,'o')
 xlabel('m')
 ylabel('Prescription opioid users')
 subplot(332)
 plot(LHSmatrix2(:,2),P_lhs2,'o')
 xlabel('beta_A')
 ylabel('Prescription opioid users')
 subplot(333)
 plot(LHSmatrix3(:,3),P_lhs3,'o')
 xlabel('beta_P')
 ylabel('Prescription opioid users')
subplot(334)
plot(LHSmatrix4(:,4),P_lhs4,'o')
xlabel('theta_1')
ylabel('Prescription opioid users')
subplot(335)
plot(LHSmatrix5(:,5),P_lhs5,'o')
xlabel('epsilon')
ylabel('Prescription opioid users')
subplot(336)
plot(LHSmatrix6(:,6),P_lhs6,'o')
xlabel('gamma')
ylabel('Prescription opioid users')
subplot(337)
plot(LHSmatrix7(:,7),P_lhs7,'o')
xlabel('sigma')
ylabel('Prescription opioid users')
subplot(338)
plot(LHSmatrix8(:,8),P_lhs8,'o')
xlabel('mu')
ylabel('Prescription opioid users')
subplot(339)
plot(LHSmatrix9(:,9),P_lhs9,'o')
xlabel('mu_A')
ylabel('Prescription opioid users') 

figure(8);
subplot(331)
 plot(LHSmatrix10(:,10),P_lhs10,'o')
 xlabel('mu_H')
 ylabel('Prescription opioid users')
 subplot(332)
 plot(LHSmatrix11(:,11),P_lhs11,'o')
 xlabel('theta_2')
 ylabel('Prescription opioid users')
 subplot(333)
 plot(LHSmatrix12(:,12),P_lhs12,'o')
 xlabel('zeta')
 ylabel('Prescription opioid users')
 subplot(334)
 plot(LHSmatrix13(:,13),P_lhs13,'o')
 xlabel('theta_3')
 ylabel('Prescription opioid users')
 subplot(335)
 plot(LHSmatrix14(:,14),P_lhs14,'o')
 xlabel('nu')
 ylabel('Prescription opioid users')
 subplot(337)
 plot(LHSmatrix15(:,15),P_lhs15,'o')
 %set(gca,'yticklabel',num2str(get(gca,'ytick')','%.9f'))
 xlabel('omega')
 ylabel('Prescription opioid users') 
 subplot(338)
 plot(LHSmatrix16(:,16),P_lhs16,'o')
 xlabel('b')
 ylabel('Prescription opioid users')
 
 
 
 
 
 figure(9);
 subplot(331)
 plot(LHSmatrix1(:,1),R_lhs1,'o')
 xlabel('m')
 ylabel('Stable Recovereds')
 subplot(332)
 plot(LHSmatrix2(:,2),R_lhs2,'o')
 xlabel('beta_A')
 ylabel('Stable Recovereds')
 subplot(333)
 plot(LHSmatrix3(:,3),R_lhs3,'o')
 xlabel('beta_P')
 ylabel('Stable Recovereds')
subplot(334)
plot(LHSmatrix4(:,4),R_lhs4,'o')
xlabel('theta_1')
ylabel('Stable Recovereds')
subplot(335)
plot(LHSmatrix5(:,5),R_lhs5,'o')
xlabel('epsilon')
ylabel('Stable Recovereds')
subplot(336)
plot(LHSmatrix6(:,6),R_lhs6,'o')
xlabel('gamma')
ylabel('Stable Recovereds')
subplot(337)
plot(LHSmatrix7(:,7),R_lhs7,'o')
xlabel('sigma')
ylabel('Stable Recovereds')
subplot(338)
plot(LHSmatrix8(:,8),R_lhs8,'o')
xlabel('mu')
ylabel('Stable Recovereds')
subplot(339)
plot(LHSmatrix9(:,9),R_lhs9,'o')
xlabel('mu_A')
ylabel('Stable Recovereds')

figure(10);
subplot(331)
 plot(LHSmatrix10(:,10),R_lhs10,'o')
 xlabel('mu_H')
 ylabel('Stable Recovereds')
 subplot(332)
 plot(LHSmatrix11(:,11),R_lhs11,'o')
 xlabel('theta_2')
 ylabel('Stable Recovereds')
 subplot(333)
 plot(LHSmatrix12(:,12),R_lhs12,'o')
 xlabel('zeta')
 ylabel('Stable Recovereds')
 subplot(334)
 plot(LHSmatrix13(:,13),R_lhs13,'o')
 xlabel('theta_3')
 ylabel('Stable Recovereds')
 subplot(335)
 plot(LHSmatrix14(:,14),R_lhs14,'o')
 xlabel('nu')
 ylabel('Stable Recovereds')
 subplot(337)
 plot(LHSmatrix15(:,15),R_lhs15,'o')
 %set(gca,'yticklabel',num2str(get(gca,'ytick')','%.9f'))
 xlabel('omega')
 ylabel('Stable Recovereds')
 subplot(338)
 plot(LHSmatrix16(:,16),R_lhs16,'o')
 xlabel('b')
 ylabel('Stable Recovereds')

 
 
 
 