
Parameter_settings_vec_Heroin;

%Number of intervals to uniformly sample zeta 
N=1000;

%Create vectors of parameters (with zeta varying) to run model with
beta_A_vec=linspace(beta_A,beta_A,N);
beta_P_vec=linspace(beta_P,beta_P,N);
theta_1_vec=linspace(0,theta_1,N);
epsilon_vec=linspace(epsilon,epsilon,N);
gamma_vec=linspace(gamma,gamma,N);
sigma_vec=linspace(sigma,sigma,N);
mu_vec=linspace(mu,mu,N);
mu_H_vec=linspace(mu_H,mu_H,N);
theta_2_vec=linspace(0,theta_2,N);
zeta_vec=linspace(zeta,zeta,N);
theta_3_vec=linspace(0,theta_3,N);
nu_vec=linspace(nu,nu,N);
omega_vec=linspace(omega,omega,N);
g_vec=linspace(g,g,N);
h_vec=linspace(h,h,N);
%b_vec=linspace(b,b,N);
%c_vec=linspace(c,c,N);
%d_vec=linspace(d,d,N);
%e_vec=linspace(e,e,N);

%Matrix of parameter sets
%vec_matrix  = [m_vec,beta_A_vec,beta_P_vec,theta_1_vec,epsilon_vec,gamma_vec,sigma_vec,mu_vec,mu_A_vec,mu_H_vec,theta_2_vec,zeta_vec,theta_3_vec,nu_vec,omega_vec,b_vec,c_vec];
%vec_matrix  = [m_vec;beta_A_vec;beta_P_vec;theta_1_vec;epsilon_vec;gamma_vec;sigma_vec;mu_vec;mu_H_vec;theta_2_vec;zeta_vec;theta_3_vec;nu_vec;omega_vec;b_vec;c_vec;d_vec;e_vec]';
vec_matrix  = [beta_A_vec;beta_P_vec;theta_1_vec;epsilon_vec;gamma_vec;sigma_vec;mu_vec;mu_H_vec;theta_2_vec;zeta_vec;theta_3_vec;nu_vec;omega_vec;g_vec;h_vec]';


%Run ODE with rows of vec_matrix 
for x=1:N;
    f=@ODE_vec_Heroin;

    %could use ode45 if issues with ode15s (so doesn't match but okay
    %because close enough) per 12/17/2019 meeting
    [t,y] = ode15s(@(t,y)f(t,y,vec_matrix,x),tspan,y0,[]); 

    W = [t y];
    
    S_lhs(:,x)=W(time_points+1,2);
    P_lhs(:,x)=W(time_points+1,3);
    A_lhs(:,x)=W(time_points+1,4);
    H_lhs(:,x)=W(time_points+1,5);
    R_lhs(:,x)=W(time_points+1,6);
    
end


 %number of prescription opioid addict overdoses at the beginning of 2023
 C=0.018606213191390.*A_lhs(1,:);
 %number of heroin/fentanyl overdoses at the beginning of 2023
 D=mu_H.*H_lhs(1,:);
 
 %number of prescription opioid addict overdoses at any time point
 %Aoverdoses=muA(t,h).*A_lhs(1,:);
 %number of heroin/fentanyl overdoses at any time point
 %Hoverdoses=mu_H.*H_lhs(1,:);
 
 %Z=muA(t,d,e).*y(:,3);
 
 %Percent change in theta_1: taking baseline theta_1 and subtracting new
 %theta_1
 %value and dividing by baseline. Then plotting final time output of A from baseline
 %theta_1-final time output of A from new theta_1 value divided by baseline final time output. Same
 %for H. 
 
 figure(1);
 plot((0.222457489109919-vec_matrix(:,3))*100/0.222457489109919,(A_lhs(1,:)-0.001681716135399)*100/0.001681716135399,'LineWidth',2) 
 hold on
 plot((0.222457489109919-vec_matrix(:,3))*100/0.222457489109919,(0.013798429999561-H_lhs(1,:))*100/0.013798429999561,'LineWidth',2)
 xlabel('Percent decrease in \theta_1')
 ylabel('Percent change in A or H at final time')
 legend({'Percent increase in A at final time', 'Percent reduction in H at final time'},'FontSize', 16)
 set(gca,'XTick',-10:10:110);
 set(gca,'FontSize',16)
 xlim([0 100])
 ylim([-5 110])
 
 %Would need to edit to include theta_2 and theta_3
 figure(2);
 plot(vec_matrix(:,3),A_lhs(1,:),'LineWidth',2) 
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\theta_1')
 ylabel('A at final time')
 set(gca,'FontSize',16)
 xlim([0 0.3])

 %Would need to edit to include theta_2 and theta_3
 figure(3);
 plot(vec_matrix(:,3),H_lhs(1,:),'-r','LineWidth',2)
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\theta_1')
 ylabel('H at final time')
 set(gca,'FontSize',16)
 xlim([0 0.3])
 
 
%  %Would need to edit to include theta_2 and theta_3
%  figure(4);
%  plot(vec_matrix(:,4),mu_H.*H_lhs(1,:),'-r','LineWidth',2)
%  %set ( gca, 'xdir', 'reverse' )
%  xlabel('\theta_1')
%  ylabel('Overdoses from H at final time')
%  set(gca,'FontSize',16)
%  xlim([0 0.3])

%  %Would need to edit to include theta_2 and theta_3
%  figure(5);
%  plot(vec_matrix(:,4),mu_H.*A_lhs(1,:),'-r','LineWidth',2)
%  %set ( gca, 'xdir', 'reverse' )
%  xlabel('\theta_1')
%  ylabel('Overdoses from H at final time')
%  set(gca,'FontSize',16)
%  xlim([0 0.3])
%  
 

%whatever entry of each vector I am interested in (anywhere from 1 to
 %1000) 
 entry=500
 
 disp('percent decrease of theta_1')
 (0.222457489109919-vec_matrix(entry,3))*100/0.222457489109919
 disp('value of theta_1')
 theta_1_vec(1,entry)
 
 disp('percent decrease of theta_2')
 (0.236479520411597-vec_matrix(entry,9))*100/0.236479520411597
 disp('value of theta_2')
 theta_2_vec(1,entry)
 
 disp('percent decrease of theta_3')
 (19.7264083013258-vec_matrix(entry,11))*100/19.7264083013258
 disp('value of theta_3')
 theta_3_vec(1,entry)
 
 disp('value of A with these theta values')
 A_lhs(1,entry)
 disp('percent increase of A from baseline 2023 value with these new theta values')
 (A_lhs(1,entry)-0.001681716135399)*100/0.001681716135399
 
 disp('value of H with these theta values')
 H_lhs(1,entry)
 disp('percent decrease of H from baseline 2023 value with these new theta values')
 (0.013798429999561-H_lhs(1,entry))*100/0.013798429999561
 
 disp('value of A overdoses at 2023 with these theta values')
 C(1,entry)
 disp('percent increase of A overdoses at 2023 with these theta values')
 (C(1,entry)-3.129036894264230e-05)*100/3.129036894264230e-05
 
 disp('value of H overdoses at 2023 with these theta values')
 D(1,entry)
 disp('percent decrease of H overdoses at 2023 with these theta values')
 (6.430068379795345e-04-D(1,entry))*100/6.430068379795345e-04
 

 
 
 
 