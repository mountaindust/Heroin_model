%H_alphapiecewise_muAlinear
clear all;
close all;

%% PRCC

%% Sample size N
 
%Total # of parameters values to test, one from each parameter interval (i.e. number of uniform intervals)
nsample = 800; 

%% LHS MATRIX  %%

%Using baseline parameter settings 
Parameter_settings_LHS_Heroin; 



%Set the vector of each parameter to be
%tested. LHS_Call_Heroin follows a method by Stein (cited in the LHS_Call
%file).
%Each parameter value will be a vector that will be arranged into a matrix
%(the latin hypercube matrix) shown below.


%From LHS_Call_Heroin:
%s=LHS_Call_Heroin(xmin,xmean,xmax,xsd,nsample,distrib,logscale)
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
mu_H_LHS=LHS_Call_Heroin(mu_H-(mu_H/2),0,mu_H+(mu_H/2),0,nsample,'unif');
theta_2_LHS=LHS_Call_Heroin(theta_2-(theta_2/2),0,theta_2+(theta_2/2),0,nsample,'unif');
zeta_LHS=LHS_Call_Heroin(zeta-(zeta/2),0,zeta+(zeta/2),0,nsample,'unif');
theta_3_LHS=LHS_Call_Heroin(theta_3-(theta_3/2),0,theta_3+(theta_3/2),0,nsample,'unif');
nu_LHS=LHS_Call_Heroin(nu-(nu/2),0,nu+(nu/2),0,nsample,'unif');
omega_LHS=LHS_Call_Heroin(omega-(omega/2),0,omega+(omega/2),0,nsample,'unif');
b_LHS=LHS_Call_Heroin(b-(b/2),0,b+(b/2),0,nsample,'unif');
c_LHS=LHS_Call_Heroin(c-(c/2),0,c+(c/2),0,nsample,'unif');
d_LHS=LHS_Call_Heroin(d-(d/2),0,d+(d/2),0,nsample,'unif');
e_LHS=LHS_Call_Heroin(e-(e/2),0,e+(e/2),0,nsample,'unif');
P0_LHS=LHS_Call_Heroin(P0-(P0/2),0,P0+(P0/2),0,nsample,'unif');
A0_LHS=LHS_Call_Heroin(A0-(A0/2),0,A0+(A0/2),0,nsample,'unif');
H0_LHS=LHS_Call_Heroin(H0-(H0/2),0,H0+(H0/2),0,nsample,'unif');
R0_LHS=LHS_Call_Heroin(R0-(R0/2),0,R0+(R0/2),0,nsample,'unif');

 

%% LHS MATRIX and PARAMETER LABELS
 %storing vectors from above in matrix form 
  LHSmatrix  = [m_LHS,beta_A_LHS,beta_P_LHS,theta_1_LHS,epsilon_LHS,gamma_LHS,sigma_LHS,mu_LHS,mu_H_LHS,theta_2_LHS,zeta_LHS,theta_3_LHS,nu_LHS,omega_LHS,b_LHS,c_LHS,d_LHS,e_LHS,P0_LHS,A0_LHS,H0_LHS,R0_LHS];

 
for x=1:nsample %Run solution nsample times choosing different values
    f=@ODE_LHS_Heroin;
%     x;
%     
%      LHSmatrix(x,:);
     y0=[1-LHSmatrix(x,19)-LHSmatrix(x,20)-LHSmatrix(x,21)-LHSmatrix(x,22),LHSmatrix(x,19),LHSmatrix(x,20),LHSmatrix(x,21),LHSmatrix(x,22)];
     [t,y] = ode15s(@(t,y)f(t,y,LHSmatrix,x),tspan,y0,[]); 
     %[t,y]=ode15s(@ODE_LHS_Heroin,tspan,y0,[],LHSmatrix);
     W = [t y]; % [time y]
    
   
    %% Save the outputs at ALL time points [tspan]
     %first column of W are the values of t 
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

  PRCC_PLOT_Heroin(LHSmatrix,H_lhs,time_points,PRCC_var,'Total Heroin/Fentanyl Addicts');
  
  
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 %% Monotonicity plots
 

%fixing all except running LHS for one parameter at a time

%m
       
m_LHS1=LHS_Call_Heroin(m-(m/2),0,m+(m/2),0,nsample,'unif');
beta_A_LHS1=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS1=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS1=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS1=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS1=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS1=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS1=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS1=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS1=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS1=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS1=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS1=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS1=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS1=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS1=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS1=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS1=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS1=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS1=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS1=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS1=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');



%beta_A

m_LHS2=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS2=LHS_Call_Heroin(beta_A-(beta_A/2),0,beta_A+(beta_A/2),0,nsample,'unif');
beta_P_LHS2=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS2=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS2=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS2=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS2=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS2=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS2=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS2=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS2=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS2=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS2=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS2=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS2=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS2=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS2=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS2=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS2=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS2=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS2=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS2=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');


%beta_P

m_LHS3=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS3=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS3=LHS_Call_Heroin(beta_P-(beta_P/2),0,beta_P+(beta_P/2),0,nsample,'unif');
theta_1_LHS3=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS3=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS3=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS3=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS3=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS3=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS3=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS3=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS3=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS3=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS3=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS3=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS3=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS3=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS3=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS3=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS3=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS3=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS3=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%theta_1

m_LHS4=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS4=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS4=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS4=LHS_Call_Heroin(theta_1-(theta_1/2),0,theta_1+(theta_1/2),0,nsample,'unif');
epsilon_LHS4=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS4=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS4=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS4=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS4=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS4=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS4=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS4=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS4=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS4=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS4=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS4=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS4=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS4=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS4=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS4=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS4=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS4=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');



%epsilon

m_LHS5=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS5=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS5=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS5=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS5=LHS_Call_Heroin(epsilon-(epsilon/2),0,epsilon+(epsilon/2),0,nsample,'unif');
gamma_LHS5=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS5=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS5=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS5=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS5=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS5=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS5=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS5=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS5=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS5=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS5=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS5=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS5=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS5=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS5=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS5=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS5=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');




%gamma

m_LHS6=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS6=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS6=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS6=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS6=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS6=LHS_Call_Heroin(gamma-(gamma/2),0,gamma+(gamma/2),0,nsample,'unif');
sigma_LHS6=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS6=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS6=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS6=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS6=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS6=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS6=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS6=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS6=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS6=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS6=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS6=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS6=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS6=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS6=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS6=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%sigma

m_LHS7=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS7=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS7=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS7=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS7=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS7=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS7=LHS_Call_Heroin(sigma-(sigma/2),0,sigma+(sigma/2),0,nsample,'unif');
mu_LHS7=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS7=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS7=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS7=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS7=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS7=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS7=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS7=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS7=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS7=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS7=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS7=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS7=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS7=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS7=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');



%mu

m_LHS8=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS8=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS8=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS8=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS8=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS8=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS8=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS8=LHS_Call_Heroin(mu-(mu/2),0,mu+(mu/2),0,nsample,'unif');
mu_H_LHS8=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS8=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS8=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS8=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS8=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS8=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS8=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS8=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS8=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS8=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS8=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS8=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS8=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS8=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');



%mu_H

m_LHS9=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS9=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS9=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS9=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS9=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS9=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS9=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS9=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS9=LHS_Call_Heroin(mu_H-(mu_H/2),0,mu_H+(mu_H/2),0,nsample,'unif');
theta_2_LHS9=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS9=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS9=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS9=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS9=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS9=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS9=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS9=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS9=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS9=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS9=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS9=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS9=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%theta_2

m_LHS10=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS10=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS10=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS10=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS10=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS10=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS10=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS10=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS10=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS10=LHS_Call_Heroin(theta_2-(theta_2/2),0,theta_2+(theta_2/2),0,nsample,'unif');
zeta_LHS10=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS10=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS10=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS10=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS10=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS10=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS10=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS10=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS10=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS10=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS10=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS10=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');


%zeta

m_LHS11=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS11=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS11=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS11=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS11=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS11=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS11=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS11=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS11=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS11=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS11=LHS_Call_Heroin(zeta-(zeta/2),0,zeta+(zeta/2),0,nsample,'unif');
theta_3_LHS11=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS11=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS11=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS11=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS11=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS11=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS11=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS11=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS11=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS11=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS11=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%theta_3

m_LHS12=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS12=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS12=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS12=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS12=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS12=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS12=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS12=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS12=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS12=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS12=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS12=LHS_Call_Heroin(theta_3-(theta_3/2),0,theta_3+(theta_3/2),0,nsample,'unif');
nu_LHS12=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS12=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS12=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS12=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS12=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS12=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS12=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS12=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS12=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS12=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%nu 

m_LHS13=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS13=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS13=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS13=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS13=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS13=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS13=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS13=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS13=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS13=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS13=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS13=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS13=LHS_Call_Heroin(nu-(nu/2),0,nu+(nu/2),0,nsample,'unif');
omega_LHS13=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS13=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS13=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS13=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS13=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS13=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS13=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS13=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS13=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%omega

m_LHS14=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS14=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS14=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS14=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS14=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS14=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS14=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS14=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS14=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS14=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS14=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS14=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS14=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS14=LHS_Call_Heroin(omega-(omega/2),0,omega+(omega/2),0,nsample,'unif');
b_LHS14=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS14=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS14=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS14=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS14=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS14=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS14=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS14=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');


%b

m_LHS15=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS15=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS15=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS15=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS15=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS15=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS15=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS15=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS15=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS15=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS15=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS15=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS15=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS15=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS15=LHS_Call_Heroin(b-(b/2),0,b+(b/2),0,nsample,'unif');
c_LHS15=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS15=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS15=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS15=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS15=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS15=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS15=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');


%c

m_LHS16=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS16=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS16=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS16=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS16=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS16=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS16=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS16=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS16=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS16=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS16=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS16=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS16=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS16=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS16=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS16=LHS_Call_Heroin(c-(c/2),0,c+(c/2),0,nsample,'unif');
d_LHS16=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS16=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS16=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS16=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS16=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS16=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

%d

m_LHS17=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS17=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS17=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS17=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS17=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS17=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS17=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS17=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS17=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS17=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS17=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS17=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS17=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS17=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS17=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS17=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS17=LHS_Call_Heroin(d-(d/2),0,d+(d/2),0,nsample,'unif');
e_LHS17=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');
P0_LHS17=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS17=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS17=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS17=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');




%e

m_LHS18=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS18=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS18=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS18=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS18=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS18=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS18=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS18=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS18=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS18=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS18=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS18=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS18=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS18=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS18=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS18=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
d_LHS18=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS18=LHS_Call_Heroin(e-(e/2),0,e+(e/2),0,nsample,'unif');
P0_LHS18=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS18=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS18=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS18=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');


%P0

m_LHS19=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS19=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS19=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS19=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS19=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS19=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS19=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS19=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS19=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS19=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS19=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS19=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS19=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS19=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS19=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS19=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
P0_LHS19=LHS_Call_Heroin(P0-(P0/2),0,P0+(P0/2),0,nsample,'unif');
A0_LHS19=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS19=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS19=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');
d_LHS19=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS19=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');

%A0

m_LHS20=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS20=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS20=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS20=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS20=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS20=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS20=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS20=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS20=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS20=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS20=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS20=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS20=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS20=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS20=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS20=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
P0_LHS20=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS20=LHS_Call_Heroin(A0-(A0/2),0,A0+(A0/2),0,nsample,'unif');
H0_LHS20=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS20=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');
d_LHS20=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS20=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');

%H0

m_LHS21=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS21=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS21=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS21=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS21=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS21=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS21=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS21=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS21=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS21=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS21=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS21=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS21=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS21=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS21=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS21=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
P0_LHS21=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS21=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS21=LHS_Call_Heroin(H0-(H0/2),0,H0+(H0/2),0,nsample,'unif');
R0_LHS21=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');
d_LHS21=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS21=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');

%R0

m_LHS22=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS22=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS22=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS22=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS22=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS22=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS22=LHS_Call_Heroin(sigma-0,0,sigma+0,0,nsample,'unif');
mu_LHS22=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_H_LHS22=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS22=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS22=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS22=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS22=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS22=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS22=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS22=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
P0_LHS22=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS22=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS22=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS22=LHS_Call_Heroin(R0-(R0/2),0,R0+(R0/2),0,nsample,'unif');
d_LHS22=LHS_Call_Heroin(d-0,0,d+0,0,nsample,'unif');
e_LHS22=LHS_Call_Heroin(e-0,0,e+0,0,nsample,'unif');






%% LHS MATRIX and PARAMETER LABELS
%taking parameter blocks above and putting them in matrix form 
  
LHSmatrix1 =[m_LHS1 beta_A_LHS1 beta_P_LHS1 theta_1_LHS1 epsilon_LHS1 gamma_LHS1 sigma_LHS1 mu_LHS1 mu_H_LHS1 theta_2_LHS1 zeta_LHS1 theta_3_LHS1 nu_LHS1 omega_LHS1 b_LHS1 c_LHS1 d_LHS1 e_LHS1 P0_LHS1 A0_LHS1 H0_LHS1 R0_LHS1];
LHSmatrix2 =[m_LHS2 beta_A_LHS2 beta_P_LHS2 theta_1_LHS2 epsilon_LHS2 gamma_LHS2 sigma_LHS2 mu_LHS2 mu_H_LHS2 theta_2_LHS2 zeta_LHS2 theta_3_LHS2 nu_LHS2 omega_LHS2 b_LHS2 c_LHS2 d_LHS2 e_LHS2 P0_LHS2 A0_LHS2 H0_LHS2 R0_LHS2];
LHSmatrix3 =[m_LHS3 beta_A_LHS3 beta_P_LHS3 theta_1_LHS3 epsilon_LHS3 gamma_LHS3 sigma_LHS3 mu_LHS3 mu_H_LHS3 theta_2_LHS3 zeta_LHS3 theta_3_LHS3 nu_LHS3 omega_LHS3 b_LHS3 c_LHS3 d_LHS3 e_LHS3 P0_LHS3 A0_LHS3 H0_LHS3 R0_LHS3];
LHSmatrix4 =[m_LHS4 beta_A_LHS4 beta_P_LHS4 theta_1_LHS4 epsilon_LHS4 gamma_LHS4 sigma_LHS4 mu_LHS4 mu_H_LHS4 theta_2_LHS4 zeta_LHS4 theta_3_LHS4 nu_LHS4 omega_LHS4 b_LHS4 c_LHS4 d_LHS4 e_LHS4 P0_LHS4 A0_LHS4 H0_LHS4 R0_LHS4];
LHSmatrix5 =[m_LHS5 beta_A_LHS5 beta_P_LHS5 theta_1_LHS5 epsilon_LHS5 gamma_LHS5 sigma_LHS5 mu_LHS5 mu_H_LHS5 theta_2_LHS5 zeta_LHS5 theta_3_LHS5 nu_LHS5 omega_LHS5 b_LHS5 c_LHS5 d_LHS5 e_LHS5 P0_LHS5 A0_LHS5 H0_LHS5 R0_LHS5]; 
LHSmatrix6 =[m_LHS6 beta_A_LHS6 beta_P_LHS6 theta_1_LHS6 epsilon_LHS6 gamma_LHS6 sigma_LHS6 mu_LHS6 mu_H_LHS6 theta_2_LHS6 zeta_LHS6 theta_3_LHS6 nu_LHS6 omega_LHS6 b_LHS6 c_LHS6 d_LHS6 e_LHS6 P0_LHS6 A0_LHS6 H0_LHS6 R0_LHS6];
LHSmatrix7 =[m_LHS7 beta_A_LHS7 beta_P_LHS7 theta_1_LHS7 epsilon_LHS7 gamma_LHS7 sigma_LHS7 mu_LHS7 mu_H_LHS7 theta_2_LHS7 zeta_LHS7 theta_3_LHS7 nu_LHS7 omega_LHS7 b_LHS7 c_LHS7 d_LHS7 e_LHS7 P0_LHS7 A0_LHS7 H0_LHS7 R0_LHS7];
LHSmatrix8 =[m_LHS8 beta_A_LHS8 beta_P_LHS8 theta_1_LHS8 epsilon_LHS8 gamma_LHS8 sigma_LHS8 mu_LHS8 mu_H_LHS8 theta_2_LHS8 zeta_LHS8 theta_3_LHS8 nu_LHS8 omega_LHS8 b_LHS8 c_LHS8 d_LHS8 e_LHS8 P0_LHS8 A0_LHS8 H0_LHS8 R0_LHS8];
LHSmatrix9 =[m_LHS9 beta_A_LHS9 beta_P_LHS9 theta_1_LHS9 epsilon_LHS9 gamma_LHS9 sigma_LHS9 mu_LHS9 mu_H_LHS9 theta_2_LHS9 zeta_LHS9 theta_3_LHS9 nu_LHS9 omega_LHS9 b_LHS9 c_LHS9 d_LHS9 e_LHS9 P0_LHS9 A0_LHS9 H0_LHS9 R0_LHS9];
LHSmatrix10 =[m_LHS10 beta_A_LHS10 beta_P_LHS10 theta_1_LHS10 epsilon_LHS10 gamma_LHS10 sigma_LHS10 mu_LHS10 mu_H_LHS10 theta_2_LHS10 zeta_LHS10 theta_3_LHS10 nu_LHS10 omega_LHS10 b_LHS10 c_LHS10 d_LHS10 e_LHS10 P0_LHS10 A0_LHS10 H0_LHS10 R0_LHS10];
LHSmatrix11 =[m_LHS11 beta_A_LHS11 beta_P_LHS11 theta_1_LHS11 epsilon_LHS11 gamma_LHS11 sigma_LHS11 mu_LHS11 mu_H_LHS11 theta_2_LHS11 zeta_LHS11 theta_3_LHS11 nu_LHS11 omega_LHS11 b_LHS11 c_LHS11 d_LHS11 e_LHS11 P0_LHS11 A0_LHS11 H0_LHS11 R0_LHS11];
LHSmatrix12 =[m_LHS12 beta_A_LHS12 beta_P_LHS12 theta_1_LHS12 epsilon_LHS12 gamma_LHS12 sigma_LHS12 mu_LHS12 mu_H_LHS12 theta_2_LHS12 zeta_LHS12 theta_3_LHS12 nu_LHS12 omega_LHS12 b_LHS12 c_LHS12 d_LHS12 e_LHS12 P0_LHS12 A0_LHS12 H0_LHS12 R0_LHS12];
LHSmatrix13 =[m_LHS13 beta_A_LHS13 beta_P_LHS13 theta_1_LHS13 epsilon_LHS13 gamma_LHS13 sigma_LHS13 mu_LHS13 mu_H_LHS13 theta_2_LHS13 zeta_LHS13 theta_3_LHS13 nu_LHS13 omega_LHS13 b_LHS13 c_LHS13 d_LHS13 e_LHS13 P0_LHS13 A0_LHS13 H0_LHS13 R0_LHS13];
LHSmatrix14 =[m_LHS14 beta_A_LHS14 beta_P_LHS14 theta_1_LHS14 epsilon_LHS14 gamma_LHS14 sigma_LHS14 mu_LHS14 mu_H_LHS14 theta_2_LHS14 zeta_LHS14 theta_3_LHS14 nu_LHS14 omega_LHS14 b_LHS14 c_LHS14 d_LHS14 e_LHS14 P0_LHS14 A0_LHS14 H0_LHS14 R0_LHS14];
LHSmatrix15 =[m_LHS15 beta_A_LHS15 beta_P_LHS15 theta_1_LHS15 epsilon_LHS15 gamma_LHS15 sigma_LHS15 mu_LHS15 mu_H_LHS15 theta_2_LHS15 zeta_LHS15 theta_3_LHS15 nu_LHS15 omega_LHS15 b_LHS15 c_LHS15 d_LHS15 e_LHS15 P0_LHS15 A0_LHS15 H0_LHS15 R0_LHS15]; 
LHSmatrix16 =[m_LHS16 beta_A_LHS16 beta_P_LHS16 theta_1_LHS16 epsilon_LHS16 gamma_LHS16 sigma_LHS16 mu_LHS16 mu_H_LHS16 theta_2_LHS16 zeta_LHS16 theta_3_LHS16 nu_LHS16 omega_LHS16 b_LHS16 c_LHS16 d_LHS16 e_LHS16 P0_LHS16 A0_LHS16 H0_LHS16 R0_LHS16];
LHSmatrix17 =[m_LHS17 beta_A_LHS17 beta_P_LHS17 theta_1_LHS17 epsilon_LHS17 gamma_LHS17 sigma_LHS17 mu_LHS17 mu_H_LHS17 theta_2_LHS17 zeta_LHS17 theta_3_LHS17 nu_LHS17 omega_LHS17 b_LHS17 c_LHS17 d_LHS17 e_LHS17 P0_LHS17 A0_LHS17 H0_LHS17 R0_LHS17];
LHSmatrix18 =[m_LHS18 beta_A_LHS18 beta_P_LHS18 theta_1_LHS18 epsilon_LHS18 gamma_LHS18 sigma_LHS18 mu_LHS18 mu_H_LHS18 theta_2_LHS18 zeta_LHS18 theta_3_LHS18 nu_LHS18 omega_LHS18 b_LHS18 c_LHS18 d_LHS18 e_LHS18 P0_LHS18 A0_LHS18 H0_LHS18 R0_LHS18];
LHSmatrix19 =[m_LHS19 beta_A_LHS19 beta_P_LHS19 theta_1_LHS19 epsilon_LHS19 gamma_LHS19 sigma_LHS19 mu_LHS19 mu_H_LHS19 theta_2_LHS19 zeta_LHS19 theta_3_LHS19 nu_LHS19 omega_LHS19 b_LHS19 c_LHS19 d_LHS19 e_LHS19 P0_LHS19 A0_LHS19 H0_LHS19 R0_LHS19];
LHSmatrix20 =[m_LHS20 beta_A_LHS20 beta_P_LHS20 theta_1_LHS20 epsilon_LHS20 gamma_LHS20 sigma_LHS20 mu_LHS20 mu_H_LHS20 theta_2_LHS20 zeta_LHS20 theta_3_LHS20 nu_LHS20 omega_LHS20 b_LHS20 c_LHS20 d_LHS20 e_LHS20 P0_LHS20 A0_LHS20 H0_LHS20 R0_LHS20];
LHSmatrix21 =[m_LHS21 beta_A_LHS21 beta_P_LHS21 theta_1_LHS21 epsilon_LHS21 gamma_LHS21 sigma_LHS21 mu_LHS21 mu_H_LHS21 theta_2_LHS21 zeta_LHS21 theta_3_LHS21 nu_LHS21 omega_LHS21 b_LHS21 c_LHS21 d_LHS21 e_LHS21 P0_LHS21 A0_LHS21 H0_LHS21 R0_LHS21];
LHSmatrix22 =[m_LHS22 beta_A_LHS22 beta_P_LHS22 theta_1_LHS22 epsilon_LHS22 gamma_LHS22 sigma_LHS22 mu_LHS22 mu_H_LHS22 theta_2_LHS22 zeta_LHS22 theta_3_LHS22 nu_LHS22 omega_LHS22 b_LHS22 c_LHS22 d_LHS22 e_LHS22 P0_LHS22 A0_LHS22 H0_LHS22 R0_LHS22];
    
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
    LHSmatrix17(x,:);
    LHSmatrix18(x,:);
    LHSmatrix19(x,:);
    LHSmatrix20(x,:);
    LHSmatrix21(x,:);
    LHSmatrix22(x,:);

    
     %Initial conditions for each LHSMatrix# solve, last four entries of LHSMatrix# rows 
    k1=[1-LHSmatrix1(x,19)-LHSmatrix1(x,20)-LHSmatrix1(x,21)-LHSmatrix1(x,22),LHSmatrix1(x,19),LHSmatrix1(x,20),LHSmatrix1(x,21),LHSmatrix1(x,22)];
    k2=[1-LHSmatrix2(x,19)-LHSmatrix2(x,20)-LHSmatrix2(x,21)-LHSmatrix2(x,22),LHSmatrix2(x,19),LHSmatrix2(x,20),LHSmatrix2(x,21),LHSmatrix2(x,22)];
    k3=[1-LHSmatrix3(x,19)-LHSmatrix3(x,20)-LHSmatrix3(x,21)-LHSmatrix3(x,22),LHSmatrix3(x,19),LHSmatrix3(x,20),LHSmatrix3(x,21),LHSmatrix3(x,22)];
    k4=[1-LHSmatrix4(x,19)-LHSmatrix4(x,20)-LHSmatrix4(x,21)-LHSmatrix4(x,22),LHSmatrix4(x,19),LHSmatrix4(x,20),LHSmatrix4(x,21),LHSmatrix4(x,22)];
    k5=[1-LHSmatrix5(x,19)-LHSmatrix5(x,20)-LHSmatrix5(x,21)-LHSmatrix5(x,22),LHSmatrix5(x,19),LHSmatrix5(x,20),LHSmatrix5(x,21),LHSmatrix5(x,22)];
    k6=[1-LHSmatrix6(x,19)-LHSmatrix6(x,20)-LHSmatrix6(x,21)-LHSmatrix6(x,22),LHSmatrix6(x,19),LHSmatrix6(x,20),LHSmatrix6(x,21),LHSmatrix6(x,22)];
    k7=[1-LHSmatrix7(x,19)-LHSmatrix7(x,20)-LHSmatrix7(x,21)-LHSmatrix7(x,22),LHSmatrix7(x,19),LHSmatrix7(x,20),LHSmatrix7(x,21),LHSmatrix7(x,22)];
    k8=[1-LHSmatrix8(x,19)-LHSmatrix8(x,20)-LHSmatrix8(x,21)-LHSmatrix8(x,22),LHSmatrix8(x,19),LHSmatrix8(x,20),LHSmatrix8(x,21),LHSmatrix8(x,22)];
    k9=[1-LHSmatrix9(x,19)-LHSmatrix9(x,20)-LHSmatrix9(x,21)-LHSmatrix9(x,22),LHSmatrix9(x,19),LHSmatrix9(x,20),LHSmatrix9(x,21),LHSmatrix9(x,22)];
    k10=[1-LHSmatrix10(x,19)-LHSmatrix10(x,20)-LHSmatrix10(x,21)-LHSmatrix10(x,22),LHSmatrix10(x,19),LHSmatrix10(x,20),LHSmatrix10(x,21),LHSmatrix10(x,22)];
    k11=[1-LHSmatrix11(x,19)-LHSmatrix11(x,20)-LHSmatrix11(x,21)-LHSmatrix11(x,22),LHSmatrix11(x,19),LHSmatrix11(x,20),LHSmatrix11(x,21),LHSmatrix11(x,22)];
    k12=[1-LHSmatrix12(x,19)-LHSmatrix12(x,20)-LHSmatrix12(x,21)-LHSmatrix12(x,22),LHSmatrix12(x,19),LHSmatrix12(x,20),LHSmatrix12(x,21),LHSmatrix12(x,22)];
    k13=[1-LHSmatrix13(x,19)-LHSmatrix13(x,20)-LHSmatrix13(x,21)-LHSmatrix13(x,22),LHSmatrix13(x,19),LHSmatrix13(x,20),LHSmatrix13(x,21),LHSmatrix13(x,22)];
    k14=[1-LHSmatrix14(x,19)-LHSmatrix14(x,20)-LHSmatrix14(x,21)-LHSmatrix14(x,22),LHSmatrix14(x,19),LHSmatrix14(x,20),LHSmatrix14(x,21),LHSmatrix14(x,22)];
    k15=[1-LHSmatrix15(x,19)-LHSmatrix15(x,20)-LHSmatrix15(x,21)-LHSmatrix15(x,22),LHSmatrix15(x,19),LHSmatrix15(x,20),LHSmatrix15(x,21),LHSmatrix15(x,22)];
    k16=[1-LHSmatrix16(x,19)-LHSmatrix16(x,20)-LHSmatrix16(x,21)-LHSmatrix16(x,22),LHSmatrix16(x,19),LHSmatrix16(x,20),LHSmatrix16(x,21),LHSmatrix16(x,22)];
    k17=[1-LHSmatrix17(x,19)-LHSmatrix17(x,20)-LHSmatrix17(x,21)-LHSmatrix17(x,22),LHSmatrix17(x,19),LHSmatrix17(x,20),LHSmatrix17(x,21),LHSmatrix17(x,22)];
    k18=[1-LHSmatrix18(x,19)-LHSmatrix18(x,20)-LHSmatrix18(x,21)-LHSmatrix18(x,22),LHSmatrix18(x,19),LHSmatrix18(x,20),LHSmatrix18(x,21),LHSmatrix18(x,22)];
    k19=[1-LHSmatrix19(x,19)-LHSmatrix19(x,20)-LHSmatrix19(x,21)-LHSmatrix19(x,22),LHSmatrix19(x,19),LHSmatrix19(x,20),LHSmatrix19(x,21),LHSmatrix19(x,22)];
    k20=[1-LHSmatrix20(x,19)-LHSmatrix20(x,20)-LHSmatrix20(x,21)-LHSmatrix20(x,22),LHSmatrix20(x,19),LHSmatrix20(x,20),LHSmatrix20(x,21),LHSmatrix20(x,22)];
    k21=[1-LHSmatrix21(x,19)-LHSmatrix21(x,20)-LHSmatrix21(x,21)-LHSmatrix21(x,22),LHSmatrix21(x,19),LHSmatrix21(x,20),LHSmatrix21(x,21),LHSmatrix21(x,22)];
    k22=[1-LHSmatrix22(x,19)-LHSmatrix22(x,20)-LHSmatrix22(x,21)-LHSmatrix22(x,22),LHSmatrix22(x,19),LHSmatrix22(x,20),LHSmatrix22(x,21),LHSmatrix22(x,22)];
   
    
    
    %run each with only 1 parameter varying in each (each ODE run nsample
    %times, 1 time for each row in LHSmatrix# to produce output for each parameter in a
    %certain range)
    [t,y1] = ode15s(@(t,y1)f(t,y1,LHSmatrix1,x),tspan,k1,[]); 
    [t,y2] = ode15s(@(t,y2)f(t,y2,LHSmatrix2,x),tspan,k2,[]); 
    [t,y3] = ode15s(@(t,y3)f(t,y3,LHSmatrix3,x),tspan,k3,[]);
    [t,y4] = ode15s(@(t,y4)f(t,y4,LHSmatrix4,x),tspan,k4,[]); 
    [t,y5] = ode15s(@(t,y5)f(t,y5,LHSmatrix5,x),tspan,k5,[]); 
    [t,y6] = ode15s(@(t,y6)f(t,y6,LHSmatrix6,x),tspan,k6,[]);
    [t,y7] = ode15s(@(t,y7)f(t,y7,LHSmatrix7,x),tspan,k7,[]); 
    [t,y8] = ode15s(@(t,y8)f(t,y8,LHSmatrix8,x),tspan,k8,[]);
    [t,y9] = ode15s(@(t,y9)f(t,y9,LHSmatrix9,x),tspan,k9,[]);
    [t,y10] = ode15s(@(t,y10)f(t,y10,LHSmatrix10,x),tspan,k10,[]);
    [t,y11] = ode15s(@(t,y11)f(t,y11,LHSmatrix11,x),tspan,k11,[]);
    [t,y12] = ode15s(@(t,y12)f(t,y12,LHSmatrix12,x),tspan,k12,[]);
    [t,y13] = ode15s(@(t,y13)f(t,y13,LHSmatrix13,x),tspan,k13,[]);
    [t,y14] = ode15s(@(t,y14)f(t,y14,LHSmatrix14,x),tspan,k14,[]);
    [t,y15] = ode15s(@(t,y15)f(t,y15,LHSmatrix15,x),tspan,k15,[]);
    [t,y16] = ode15s(@(t,y16)f(t,y16,LHSmatrix16,x),tspan,k16,[]);
    [t,y17] = ode15s(@(t,y17)f(t,y17,LHSmatrix17,x),tspan,k17,[]);
    [t,y18] = ode15s(@(t,y18)f(t,y18,LHSmatrix18,x),tspan,k18,[]);
    [t,y19] = ode15s(@(t,y19)f(t,y19,LHSmatrix19,x),tspan,k19,[]);
    [t,y20] = ode15s(@(t,y20)f(t,y20,LHSmatrix20,x),tspan,k20,[]);
    [t,y21] = ode15s(@(t,y21)f(t,y21,LHSmatrix21,x),tspan,k21,[]);
    [t,y22] = ode15s(@(t,y22)f(t,y22,LHSmatrix22,x),tspan,k22,[]);

    
    
     %store results of each [t,y1] = ode15s(@(t,y1)f(t,y1,LHSmatrix1,x),tspan,y0,[]);
     %These get overwritten for each x of the for loop, so the final result
     %in the workspace shown is for the final x value (the last row of each LHS matrix) over the
     %entire time span 
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
     W17 = [t y17];
     W18 = [t y18];
     W19 = [t y19];
     W20 = [t y20];
     W21 = [t y21];
     W22 = [t y22];
 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %store final value of each class (time_points+1 is final time since first
  %column is IC's), start with 2nd column because 1st just for time values 
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
     
     S_lhs17(:,x)=W17(time_points+1,2);
     P_lhs17(:,x)=W17(time_points+1,3);
     A_lhs17(:,x)=W17(time_points+1,4);
     H_lhs17(:,x)=W17(time_points+1,5);
     R_lhs17(:,x)=W17(time_points+1,6);
     
     S_lhs18(:,x)=W18(time_points+1,2);
     P_lhs18(:,x)=W18(time_points+1,3);
     A_lhs18(:,x)=W18(time_points+1,4);
     H_lhs18(:,x)=W18(time_points+1,5);
     R_lhs18(:,x)=W18(time_points+1,6);
     
     S_lhs19(:,x)=W19(time_points+1,2);
     P_lhs19(:,x)=W19(time_points+1,3);
     A_lhs19(:,x)=W19(time_points+1,4);
     H_lhs19(:,x)=W19(time_points+1,5);
     R_lhs19(:,x)=W19(time_points+1,6);
    
     
     S_lhs20(:,x)=W20(time_points+1,2);
     P_lhs20(:,x)=W20(time_points+1,3);
     A_lhs20(:,x)=W20(time_points+1,4);
     H_lhs20(:,x)=W20(time_points+1,5);
     R_lhs20(:,x)=W20(time_points+1,6);
    
     
     S_lhs21(:,x)=W21(time_points+1,2);
     P_lhs21(:,x)=W21(time_points+1,3);
     A_lhs21(:,x)=W21(time_points+1,4);
     H_lhs21(:,x)=W21(time_points+1,5);
     R_lhs21(:,x)=W21(time_points+1,6);
     
     
     S_lhs22(:,x)=W22(time_points+1,2);
     P_lhs22(:,x)=W22(time_points+1,3);
     A_lhs22(:,x)=W22(time_points+1,4);
     H_lhs22(:,x)=W22(time_points+1,5);
     R_lhs22(:,x)=W22(time_points+1,6);
    
    
    
    
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
[prcc17 sign17 sign_label17] = PRCC_Heroin(LHSmatrix17,H_lhs17,1:length(time_points),PRCC_var,alpha);
[prcc18 sign18 sign_label18] = PRCC_Heroin(LHSmatrix18,H_lhs18,1:length(time_points),PRCC_var,alpha);
[prcc19 sign19 sign_label19] = PRCC_Heroin(LHSmatrix19,H_lhs19,1:length(time_points),PRCC_var,alpha);
[prcc20 sign20 sign_label20] = PRCC_Heroin(LHSmatrix20,H_lhs20,1:length(time_points),PRCC_var,alpha);
[prcc21 sign21 sign_label21] = PRCC_Heroin(LHSmatrix21,H_lhs21,1:length(time_points),PRCC_var,alpha);
[prcc22 sign22 sign_label22] = PRCC_Heroin(LHSmatrix22,H_lhs22,1:length(time_points),PRCC_var,alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Monotonicity curves 


 figure(1);
 subplot(421)
 plot(LHSmatrix1(:,1),S_lhs1,'o') %(i.e. first column of LHSmatrix1 is m varying, and S_lhs1 is ODE output from using those values (while all other parameters are fixed), so this is how S is affected
 xlabel('m')
 ylabel('S')
 subplot(422)
 plot(LHSmatrix2(:,2),S_lhs2,'o')
 xlabel('\beta_A')
 ylabel('S')
 subplot(423)
 plot(LHSmatrix3(:,3),S_lhs3,'o')
 xlabel('\beta_P')
 ylabel('S')
subplot(424)
plot(LHSmatrix4(:,4),S_lhs4,'o')
xlabel('\theta_1')
ylabel('S')
subplot(425)
plot(LHSmatrix5(:,5),S_lhs5,'o')
xlabel('\epsilon')
ylabel('S')
subplot(426)
plot(LHSmatrix6(:,6),S_lhs6,'o')
xlabel('\gamma')
ylabel('S')
subplot(427)
plot(LHSmatrix7(:,7),S_lhs7,'o')
xlabel('\sigma')
ylabel('S')
subplot(428)
plot(LHSmatrix8(:,8),S_lhs8,'o')
xlabel('\mu')
ylabel('S')

figure(2);
subplot(421)
plot(LHSmatrix9(:,9),S_lhs9,'o')
xlabel('\mu_H')
ylabel('S')
subplot(422)
 plot(LHSmatrix10(:,10),S_lhs10,'o')
 xlabel('\theta_2')
 ylabel('S')
 subplot(423)
 plot(LHSmatrix11(:,11),S_lhs11,'o')
 xlabel('\zeta')
 ylabel('S')
 subplot(424)
 plot(LHSmatrix12(:,12),S_lhs12,'o')
 xlabel('\theta_3')
 ylabel('S')
 subplot(425)
 plot(LHSmatrix13(:,13),S_lhs13,'o')
 xlabel('\nu')
 ylabel('S')
 subplot(426)
 plot(LHSmatrix14(:,14),S_lhs14,'o')
 xlabel('\omega')
 ylabel('S')
 subplot(427)
 plot(LHSmatrix15(:,15),S_lhs15,'o')
 xlabel('b')
 ylabel('S')
 subplot(428)
 plot(LHSmatrix16(:,16),S_lhs16,'o')
 xlabel('c')
 ylabel('S')
 
 figure(3);
 subplot(421)
 plot(LHSmatrix17(:,17),S_lhs17,'o')
 xlabel('d')
 ylabel('S')
 subplot(422)
 plot(LHSmatrix18(:,18),S_lhs18,'o')
 xlabel('e')
 ylabel('S')
 subplot(423)
 plot(LHSmatrix19(:,19),S_lhs19,'o')
 xlabel('P0')
 ylabel('S')
 subplot(424)
 plot(LHSmatrix20(:,20),S_lhs20,'o')
 xlabel('A0')
 ylabel('S')
 subplot(425)
 plot(LHSmatrix21(:,21),S_lhs21,'o')
 xlabel('H0')
 ylabel('S')
 subplot(426)
 plot(LHSmatrix22(:,22),S_lhs22,'o')
 xlabel('R0')
 ylabel('S')

 
 
 %Monotonicity curves for prescription opioids users at last time step 
 
  
 figure(4);
 subplot(421)
 plot(LHSmatrix1(:,1),P_lhs1,'o')
 xlabel('m')
 ylabel('P')
 subplot(422)
 plot(LHSmatrix2(:,2),P_lhs2,'o')
 xlabel('\beta_A')
 ylabel('P')
 subplot(423)
 plot(LHSmatrix3(:,3),P_lhs3,'o')
 xlabel('\beta_P')
 ylabel('P')
subplot(424)
plot(LHSmatrix4(:,4),P_lhs4,'o')
xlabel('\theta_1')
ylabel('P')
subplot(425)
plot(LHSmatrix5(:,5),P_lhs5,'o')
xlabel('\epsilon')
ylabel('P')
subplot(426)
plot(LHSmatrix6(:,6),P_lhs6,'o')
xlabel('\gamma')
ylabel('P')
subplot(427)
plot(LHSmatrix7(:,7),P_lhs7,'o')
xlabel('\sigma')
ylabel('P')
subplot(428)
plot(LHSmatrix8(:,8),P_lhs8,'o')
xlabel('\mu')
ylabel('P')


figure(5);
subplot(421)
plot(LHSmatrix9(:,9),P_lhs9,'o')
xlabel('\mu_H')
ylabel('P') 
subplot(422)
 plot(LHSmatrix10(:,10),P_lhs10,'o')
 xlabel('\theta_2')
 ylabel('P')
 subplot(423)
 plot(LHSmatrix11(:,11),P_lhs11,'o')
 xlabel('\zeta')
 ylabel('P')
 subplot(424)
 plot(LHSmatrix12(:,12),P_lhs12,'o')
 xlabel('\theta_3')
 ylabel('P')
 subplot(425)
 plot(LHSmatrix13(:,13),P_lhs13,'o')
 xlabel('\nu')
 ylabel('P')
 subplot(426)
 plot(LHSmatrix14(:,14),P_lhs14,'o')
 xlabel('\omega')
 ylabel('P')
 subplot(427)
 plot(LHSmatrix15(:,15),P_lhs15,'o')
 %set(gca,'yticklabel',num2str(get(gca,'ytick')','%.9f'))
 xlabel('b')
 ylabel('P') 
 subplot(428)
 plot(LHSmatrix16(:,16),P_lhs16,'o')
 xlabel('c')
 ylabel('P')
 
 figure(6);
 subplot(421)
 plot(LHSmatrix17(:,17),P_lhs17,'o')
 xlabel('d')
 ylabel('P')
 subplot(422)
 plot(LHSmatrix18(:,18),P_lhs18,'o')
 xlabel('e')
 ylabel('P')
 subplot(423)
 plot(LHSmatrix19(:,19),P_lhs19,'o')
 xlabel('P0')
 ylabel('P')
 subplot(424)
 plot(LHSmatrix20(:,20),P_lhs20,'o')
 xlabel('A0')
 ylabel('P')
 subplot(425)
 plot(LHSmatrix21(:,21),P_lhs21,'o')
 xlabel('H0')
 ylabel('P')
 subplot(426)
 plot(LHSmatrix22(:,22),P_lhs22,'o')
 xlabel('R0')
 ylabel('P')
 
 
 
 %Monotonicity curves for prescription addicted individuals at last time step 
 
 figure(7);
 subplot(421)
 %take first column of LHSmatrix1 because those are the values m is varying over to produce the output in A_lhs1
 plot(LHSmatrix1(:,1),A_lhs1,'o') 
 xlabel('m')
 ylabel('A')
 subplot(422)
 plot(LHSmatrix2(:,2),A_lhs2,'o')
 xlabel('\beta_A')
 ylabel('A')
 subplot(423)
 plot(LHSmatrix3(:,3),A_lhs3,'o')
 xlabel('\beta_P')
 ylabel('A')
subplot(424)
plot(LHSmatrix4(:,4),A_lhs4,'o')
xlabel('\theta_1')
ylabel('A')
subplot(425)
plot(LHSmatrix5(:,5),A_lhs5,'o')
xlabel('\epsilon')
ylabel('A')
subplot(426)
plot(LHSmatrix6(:,6),A_lhs6,'o')
xlabel('\gamma')
ylabel('A')
subplot(427)
plot(LHSmatrix7(:,7),A_lhs7,'o')
xlabel('\sigma')
ylabel('A')
subplot(428)
plot(LHSmatrix8(:,8),A_lhs8,'o')
xlabel('\mu')
ylabel('A')

figure(8);
subplot(421)
plot(LHSmatrix9(:,9),A_lhs9,'o')
xlabel('\mu_H')
ylabel('A') 
subplot(422)
 plot(LHSmatrix10(:,10),A_lhs10,'o')
 xlabel('\theta_2')
 ylabel('A')
 subplot(423)
 plot(LHSmatrix11(:,11),A_lhs11,'o')
 xlabel('\zeta')
 ylabel('A')
 subplot(424)
 plot(LHSmatrix12(:,12),A_lhs12,'o')
 xlabel('\theta_3')
 ylabel('A')
 subplot(425)
 plot(LHSmatrix13(:,13),A_lhs13,'o')
 xlabel('\nu')
 ylabel('A')
 subplot(426)
 plot(LHSmatrix14(:,14),A_lhs14,'o')
 xlabel('\omega')
 ylabel('A')
 subplot(427)
 plot(LHSmatrix15(:,15),A_lhs15,'o')
 xlabel('b')
 ylabel('A')
 subplot(428)
 plot(LHSmatrix16(:,16),A_lhs16,'o')
 xlabel('c')
 ylabel('A')
 
 figure(9);
 subplot(421)
 plot(LHSmatrix17(:,17),A_lhs17,'o')
 xlabel('d')
 ylabel('A')
 subplot(422)
 plot(LHSmatrix18(:,18),A_lhs18,'o')
 xlabel('e')
 ylabel('A')
 subplot(423)
 plot(LHSmatrix19(:,19),A_lhs19,'o')
 xlabel('P0')
 ylabel('A')
 subplot(424)
 plot(LHSmatrix20(:,20),A_lhs20,'o')
 xlabel('A0')
 ylabel('A')
 subplot(425)
 plot(LHSmatrix21(:,21),A_lhs21,'o')
 xlabel('H0')
 ylabel('A')
 subplot(426)
 plot(LHSmatrix22(:,22),A_lhs22,'o')
 xlabel('R0')
 ylabel('A')

 
 
 %Monotonicity curves for heroin addicted individuals at last time step 
 
 figure(10);
 subplot(421)
 plot(LHSmatrix1(:,1),H_lhs1,'o')
 xlabel('m')
 ylabel('H')
 subplot(422)
 plot(LHSmatrix2(:,2),H_lhs2,'o')
 xlabel('\beta_A')
 ylabel('H')
 subplot(423)
 plot(LHSmatrix3(:,3),H_lhs3,'o')
 xlabel('\beta_P')
 ylabel('H')
subplot(424)
plot(LHSmatrix4(:,4),H_lhs4,'o')
xlabel('\theta_1')
ylabel('H')
subplot(425)
plot(LHSmatrix5(:,5),H_lhs5,'o')
xlabel('\epsilon')
ylabel('H')
subplot(426)
plot(LHSmatrix6(:,6),H_lhs6,'o')
xlabel('\gamma')
ylabel('H')
subplot(427)
plot(LHSmatrix7(:,7),H_lhs7,'o')
xlabel('\sigma')
ylabel('H')
subplot(428)
plot(LHSmatrix8(:,8),H_lhs8,'o')
xlabel('\mu')
ylabel('H')

figure(11);
subplot(421)
plot(LHSmatrix9(:,9),H_lhs9,'o')
xlabel('\mu_H')
ylabel('H') 
subplot(422)
 plot(LHSmatrix10(:,10),H_lhs10,'o')
 xlabel('\theta_2')
 ylabel('H')
 subplot(423)
 plot(LHSmatrix11(:,11),H_lhs11,'o')
 xlabel('\zeta')
 ylabel('H')
 subplot(424)
 plot(LHSmatrix12(:,12),H_lhs12,'o')
 xlabel('\theta_3')
 ylabel('H')
 subplot(425)
 plot(LHSmatrix13(:,13),H_lhs13,'o')
 xlabel('\nu')
 ylabel('H')
 subplot(426)
 plot(LHSmatrix14(:,14),H_lhs14,'o')
 xlabel('\omega')
 ylabel('H')
 subplot(427)
 plot(LHSmatrix15(:,15),H_lhs15,'o')
 xlabel('b')
 ylabel('H') 
 subplot(428)
 plot(LHSmatrix16(:,16),H_lhs16,'o')
 xlabel('c')
 ylabel('H')
 
 figure(12);
 subplot(421)
 plot(LHSmatrix17(:,17),H_lhs17,'o')
 xlabel('d')
 ylabel('H')
 subplot(422)
 plot(LHSmatrix18(:,18),H_lhs18,'o')
 xlabel('e')
 ylabel('H')
 subplot(423)
 plot(LHSmatrix19(:,19),H_lhs19,'o')
 xlabel('P0')
 ylabel('H')
 subplot(424)
 plot(LHSmatrix20(:,20),H_lhs20,'o')
 xlabel('A0')
 ylabel('H')
 subplot(425)
 plot(LHSmatrix21(:,21),H_lhs21,'o')
 xlabel('H0')
 ylabel('H')
 subplot(426)
 plot(LHSmatrix22(:,22),H_lhs22,'o')
 xlabel('R0')
 ylabel('H')

 
 
 
%Monotonicity curves for stably recovered individuals at last time step 
  
 figure(13);
 subplot(421)
 plot(LHSmatrix1(:,1),R_lhs1,'o')
 xlabel('m')
 ylabel('R')
 subplot(422)
 plot(LHSmatrix2(:,2),R_lhs2,'o')
 xlabel('\beta_A')
 ylabel('R')
 subplot(423)
 plot(LHSmatrix3(:,3),R_lhs3,'o')
 xlabel('\beta_P')
 ylabel('R')
subplot(424)
plot(LHSmatrix4(:,4),R_lhs4,'o')
xlabel('\theta_1')
ylabel('R')
subplot(425)
plot(LHSmatrix5(:,5),R_lhs5,'o')
xlabel('\epsilon')
ylabel('R')
subplot(426)
plot(LHSmatrix6(:,6),R_lhs6,'o')
xlabel('\gamma')
ylabel('R')
subplot(427)
plot(LHSmatrix7(:,7),R_lhs7,'o')
xlabel('\sigma')
ylabel('R')
subplot(428)
plot(LHSmatrix8(:,8),R_lhs8,'o')
xlabel('\mu')
ylabel('R')

figure(14);
subplot(421)
plot(LHSmatrix9(:,9),R_lhs9,'o')
xlabel('\mu_H')
ylabel('R')
subplot(422)
 plot(LHSmatrix10(:,10),R_lhs10,'o')
 xlabel('\theta_2')
 ylabel('R')
 subplot(423)
 plot(LHSmatrix11(:,11),R_lhs11,'o')
 xlabel('\zeta')
 ylabel('R')
 subplot(424)
 plot(LHSmatrix12(:,12),R_lhs12,'o')
 xlabel('\theta_3')
 ylabel('R')
 subplot(425)
 plot(LHSmatrix13(:,13),R_lhs13,'o')
 xlabel('\nu')
 ylabel('R')
 subplot(426)
 plot(LHSmatrix14(:,14),R_lhs14,'o')
 xlabel('\omega')
 ylabel('R')
 subplot(427)
 plot(LHSmatrix15(:,15),R_lhs15,'o')
 xlabel('b')
 ylabel('R')
 subplot(428)
 plot(LHSmatrix16(:,16),R_lhs16,'o')
 xlabel('c')
 ylabel('R')
 
 figure(15);
 subplot(421)
 plot(LHSmatrix17(:,17),R_lhs17,'o')
 xlabel('d')
 ylabel('R')
 subplot(422)
 plot(LHSmatrix18(:,18),R_lhs18,'o')
 xlabel('e')
 ylabel('R')
 subplot(423)
 plot(LHSmatrix19(:,19),R_lhs19,'o')
 xlabel('P0')
 ylabel('R')
 subplot(424)
 plot(LHSmatrix20(:,20),R_lhs20,'o')
 xlabel('A0')
 ylabel('R')
 subplot(425)
 plot(LHSmatrix21(:,21),R_lhs21,'o')
 xlabel('H0')
 ylabel('R')
 subplot(426)
 plot(LHSmatrix22(:,22),R_lhs22,'o')
 xlabel('R0')
 ylabel('R')

 