%FIND ERROR FOR HOW CHANGING ZETA AFFECT A PLOT!!!

Parameter_settings_vec_Heroin;

%Number of intervals to uniformly sample zeta 
N=100;

%Create vectors of parameters (with zeta varying) to run model with
m_vec=linspace(m,m,N);
beta_A_vec=linspace(beta_A,beta_A,N);
beta_P_vec=linspace(beta_P,beta_P,N);
theta_1_vec=linspace(theta_1,theta_1,N);
epsilon_vec=linspace(epsilon, epsilon,N);
gamma_vec=linspace(gamma,gamma,N);
sigma_vec=linspace(sigma,sigma,N);
mu_vec=linspace(mu,mu,N);
mu_H_vec=linspace(mu_H,mu_H,N);
theta_2_vec=linspace(theta_2,theta_2,N);
zeta_vec=linspace(0.1*zeta,zeta,N);
theta_3_vec=linspace(theta_3,theta_3,N);
nu_vec=linspace(nu,nu,N);
omega_vec=linspace(omega,omega,N);
b_vec=linspace(b,b,N);
c_vec=linspace(c,c,N);
d_vec=linspace(d,d,N);
e_vec=linspace(e,e,N);

%Matrix of parameter sets
%vec_matrix  = [m_vec,beta_A_vec,beta_P_vec,theta_1_vec,epsilon_vec,gamma_vec,sigma_vec,mu_vec,mu_A_vec,mu_H_vec,theta_2_vec,zeta_vec,theta_3_vec,nu_vec,omega_vec,b_vec,c_vec];
vec_matrix  = [m_vec;beta_A_vec;beta_P_vec;theta_1_vec;epsilon_vec;gamma_vec;sigma_vec;mu_vec;mu_H_vec;theta_2_vec;zeta_vec;theta_3_vec;nu_vec;omega_vec;b_vec;c_vec;d_vec;e_vec]';


%Run ODE with rows of vec_matrix 
for x=1:N;
    f=@ODE_vec_Heroin;


    [t,y] = ode15s(@(t,y)f(t,y,vec_matrix,x),tspan,y0,[]); 

    W = [t y];
    
    S_lhs(:,x)=W(time_points+1,2);
    P_lhs(:,x)=W(time_points+1,3);
    A_lhs(:,x)=W(time_points+1,4);
    H_lhs(:,x)=W(time_points+1,5);
    R_lhs(:,x)=W(time_points+1,6);
    
    
end  
    
 %Percent change in zeta: taking baseline zeta and subtracting new zeta
 %value and dividing by baseline. Then plotting final time output of A from baseline
 %zeta-final time output of A from new zeta value divided by baseline final time output. Same
 %for H. 
 figure(1);
 hold all
 plot((0.198182427387906-vec_matrix(:,11))*100/0.198182427387906,(0.004329295236935-A_lhs(1,:))*100/0.004329295236935,'LineWidth',2) 
 hold on
 plot((0.198182427387906-vec_matrix(:,11))*100/0.198182427387906,(0.004298941991162-H_lhs(1,:))*100/0.004298941991162,'LineWidth',2)
 xlabel('Percent reduction in \zeta')
 ylabel('Percent reduction in A or H')
 %legend({'Percent reduction in A', 'Percent reduction in H'},'FontSize', 16)
 set(gca,'XTick',0:10:100);
 set(gca,'FontSize',16)
 xlim([0 60])
 ylim([0 100])
 
 figure(2);
 plot(vec_matrix(:,11),A_lhs(1,:),'LineWidth',2) 
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\zeta')
 ylabel('A')
 set(gca,'FontSize',16)
 xlim([0 0.2])

 figure(3);
 plot(vec_matrix(:,11),H_lhs(1,:),'-r','LineWidth',2)
 %set ( gca, 'xdir', 'reverse' )
 xlabel('\zeta')
 ylabel('H')
 set(gca,'FontSize',16)
 xlim([0 0.2])
 
 