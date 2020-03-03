
Parameter_settings_vec_Heroin;

%Number of intervals to uniformly sample zeta 
N=1000;

%Create vectors of parameters (with zeta varying) to run model with
%m_vec=linspace(m,m,N);
beta_A_vec=linspace(beta_A,beta_A,N);
beta_P_vec=linspace(beta_P,beta_P,N);
theta_1_vec=linspace(theta_1,theta_1,N);
epsilon_vec=linspace(epsilon,epsilon,N);
gamma_vec=linspace(gamma,gamma,N);
sigma_vec=linspace(sigma,sigma,N);
mu_vec=linspace(mu,mu,N);
mu_H_vec=linspace(mu_H,mu_H,N);
theta_2_vec=linspace(theta_2,theta_2,N);
zeta_vec=linspace(zeta,1.5*zeta,N);
theta_3_vec=linspace(theta_3,theta_3,N);
nu_vec=linspace(nu,1.5*nu,N);
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

    %using ode45 because issues with ode15s (so doesn't match but okay
    %because close enough) per 12/17/2019 meeting
    [t,y] = ode45(@(t,y)f(t,y,vec_matrix,x),tspan,y0,[]); 

    W = [t y];
    
    S_lhs(:,x)=W(time_points+1,2);
    P_lhs(:,x)=W(time_points+1,3);
    A_lhs(:,x)=W(time_points+1,4);
    H_lhs(:,x)=W(time_points+1,5);
    R_lhs(:,x)=W(time_points+1,6);
    
    
end  
    
 %Percent change in zeta: taking new zeta and subtracting baseline zeta
 %value and dividing by baseline. Then plotting final time output of A from baseline
 %zeta-final time output of A from new zeta value divided by baseline final time output. Same
 %for H. 
 
 figure(1);
 plot((vec_matrix(:,10)-0.198182427387906)*100/0.198182427387906,(0.001681923193256-A_lhs(1,:))*100/0.001681923193256,'LineWidth',2) 
 hold on
 plot((vec_matrix(:,10)-0.198182427387906)*100/0.198182427387906,(0.013796340958938-H_lhs(1,:))*100/0.013796340958938,'LineWidth',2)
 xlabel('Percent increase in \zeta and \nu')
 ylabel('Percent change in A or H at final time')
 legend({'Percent reduction in A at final time', 'Percent reduction in H at final time'},'FontSize', 16)
 set(gca,'XTick',-10:10:110);
 set(gca,'FontSize',16)
 xlim([0 100])
 ylim([-5 110])
 
 
 %Would need to edit to include change in nu, as well
 figure(2);
 plot(vec_matrix(:,10),A_lhs(1,:),'LineWidth',2) 
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\zeta')
 ylabel('A at final time')
 set(gca,'FontSize',16)
 xlim([0.2 0.4])

 %Would need to edit to include change in nu, as well
 figure(3);
 plot(vec_matrix(:,10),H_lhs(1,:),'-r','LineWidth',2)
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\zeta')
 ylabel('H at final time')
 set(gca,'FontSize',16)
 xlim([0.2 0.4])
 
 
 %number of prescription opioid addict overdoses at the final time
 W=0.018606213191390.*A_lhs(1,:);
 %number of heroin/fentanyl overdoses at the final time
 Z=mu_H.*H_lhs(1,:);
 
 
 
 figure(4);
 plot((vec_matrix(:,10)-0.198182427387906)*100/0.198182427387906,(3.129422150525879e-05-W)*100/3.129422150525879e-05,'LineWidth',2) 
 hold on
 plot((vec_matrix(:,10)-0.198182427387906)*100/0.198182427387906,(6.429094886864968e-04-Z)*100/6.429094886864968e-04,'LineWidth',2)
 xlabel('Percent increase in \zeta and \nu')
 ylabel('Percent change in A or H overdoses at final time')
 legend({'Percent reduction in A overdoses at final time', 'Percent reduction in H overdoses at final time'},'FontSize', 16)
 set(gca,'XTick',-10:10:110);
 set(gca,'FontSize',16)
 xlim([0 100])
 ylim([-5 110])
 
 
 %Would need to edit to include change in nu, as well
 figure(5);
 plot(vec_matrix(:,10),W,'LineWidth',2) 
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\zeta')
 ylabel('A overdoses at final time')
 set(gca,'FontSize',16)
 xlim([0.2 0.4])

 %Would need to edit to include change in nu, as well
 figure(6);
 plot(vec_matrix(:,10),Z,'-r','LineWidth',2)
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\zeta')
 ylabel('H overdoses at final time')
 set(gca,'FontSize',16)
 xlim([0.2 0.4])
 
 