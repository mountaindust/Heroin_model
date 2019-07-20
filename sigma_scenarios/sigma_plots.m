%sigma varying
clf;
clear all;

Parameter_settings_LHS_Heroin; 
nsample = 500; 

%sigma

m_LHS=LHS_Call_Heroin(m-0,0,m+0,0,nsample,'unif');
beta_A_LHS=LHS_Call_Heroin(beta_A-0,0,beta_A+0,0,nsample,'unif');
beta_P_LHS=LHS_Call_Heroin(beta_P-0,0,beta_P+0,0,nsample,'unif');
theta_1_LHS=LHS_Call_Heroin(theta_1-0,0,theta_1+0,0,nsample,'unif');
epsilon_LHS=LHS_Call_Heroin(epsilon-0,0,epsilon+0,0,nsample,'unif');
gamma_LHS=LHS_Call_Heroin(gamma-0,0,gamma+0,0,nsample,'unif');
sigma_LHS=LHS_Call_Heroin(sigma*0.1,0,sigma,0,nsample,'unif');
mu_LHS=LHS_Call_Heroin(mu-0,0,mu+0,0,nsample,'unif');
mu_A_LHS=LHS_Call_Heroin(mu_A-0,0,mu_A+0,0,nsample,'unif');
mu_H_LHS=LHS_Call_Heroin(mu_H-0,0,mu_H+0,0,nsample,'unif');
theta_2_LHS=LHS_Call_Heroin(theta_2-0,0,theta_2+0,0,nsample,'unif');
zeta_LHS=LHS_Call_Heroin(zeta-0,0,zeta+0,0,nsample,'unif');
theta_3_LHS=LHS_Call_Heroin(theta_3-0,0,theta_3+0,0,nsample,'unif');
nu_LHS=LHS_Call_Heroin(nu-0,0,nu+0,0,nsample,'unif');
omega_LHS=LHS_Call_Heroin(omega-0,0,omega+0,0,nsample,'unif');
b_LHS=LHS_Call_Heroin(b-0,0,b+0,0,nsample,'unif');
c_LHS=LHS_Call_Heroin(c-0,0,c+0,0,nsample,'unif');
P0_LHS=LHS_Call_Heroin(P0-0,0,P0+0,0,nsample,'unif');
A0_LHS=LHS_Call_Heroin(A0-0,0,A0+0,0,nsample,'unif');
H0_LHS=LHS_Call_Heroin(H0-0,0,H0+0,0,nsample,'unif');
R0_LHS=LHS_Call_Heroin(R0-0,0,R0+0,0,nsample,'unif');

LHSmatrix  = [m_LHS,beta_A_LHS,beta_P_LHS,theta_1_LHS,epsilon_LHS,gamma_LHS,sigma_LHS,mu_LHS,mu_A_LHS,mu_H_LHS,theta_2_LHS,zeta_LHS,theta_3_LHS,nu_LHS,omega_LHS,b_LHS,c_LHS,P0_LHS,A0_LHS,H0_LHS,R0_LHS];


for x=1:nsample %Run solution x times choosing different values, represents each row of the matrix that's going to go through the ODE solver to produce the plots 
    f=@ODE_LHS_Heroin;
    x;
     
     %values to test
    LHSmatrix(x,:);
   
    
    
    
    %run each with only 1 parameter varying in each (each ODE run nsample
    %times, 1 time for each row in LHSmatrix# to produce output for each parameter in a
    %certain range)
    [t,y] = ode15s(@(t,y)f(t,y,LHSmatrix,x),tspan,y0,[]); 
   
 
    
    
     %store results of each [t,y1] = ode15s(@(t,y1)f(t,y1,LHSmatrix1,x),tspan,y0,[]);
     %These get overwritten for each x of the for loop, so the final result
     %in the workspace shown is for the final x value (the last row of each LHS matrix) over the
     %entire time span 
    W = [t y]; % [time y]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %store final value of each class (time_points+1 is final time since first
  %column is IC's), start with 2nd column because 1st just for time values 
    S_lhs(:,x)=W(time_points+1,2);
    P_lhs(:,x)=W(time_points+1,3);
    A_lhs(:,x)=W(time_points+1,4);
    H_lhs(:,x)=W(time_points+1,5);
    R_lhs(:,x)=W(time_points+1,6);

    
end

figure(1);
 subplot(221)
 plot(LHSmatrix(:,7),A_lhs,'o') %(i.e. first column of LHSmatrix1 is m varying, and S_lhs1 is ODE output from using those values (while all other parameters are fixed), so this is how S is affected
 xlabel('\sigma')
 ylabel('A')
 subplot(222)
 plot(LHSmatrix(:,7),H_lhs,'o')
 xlabel('\sigma')
 ylabel('H')
 

figure(2);
 subplot(221)
 plot(LHSmatrix(:,7),A_lhs,'o') %(i.e. first column of LHSmatrix1 is m varying, and S_lhs1 is ODE output from using those values (while all other parameters are fixed), so this is how S is affected
 set ( gca, 'xdir', 'reverse' )
 xlabel('\sigma')
 ylabel('A')
 subplot(222)
 plot(LHSmatrix(:,7),H_lhs,'o')
 set ( gca, 'xdir', 'reverse' )
 xlabel('\sigma')
 ylabel('H')

  
 figure(3);
 plot((0.0283-LHSmatrix(:,7))*100/0.0283,(0.004418807691222-A_lhs)*100/0.004418807691222,'o') %(i.e. first column of LHSmatrix1 is m varying, and S_lhs1 is ODE output from using those values (while all other parameters are fixed), so this is how S is affected
 hold on
 plot((0.0283-LHSmatrix(:,7))*100/0.0283,(0.002197350239427-H_lhs)*100/0.002197350239427,'o')
 hold off
 xlabel('Percent reduction in \sigma')
 ylabel('Percent reduction in A or H')
 legend({'Percent reduction in A', 'Percent reduction in H'},'FontSize', 10)
 set(gca,'XTick',0:10:100);
 xlim([0 90])
 ylim([0 85])
 
 
 


