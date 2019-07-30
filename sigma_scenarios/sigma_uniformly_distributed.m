%NEED TO FIX!!!

%Baseline parameter values
m=-0.00483;
beta_A=0.0044;
beta_P=0.000469;
theta_1=0.000502;
epsilon=2.49;
gamma=0.00146;
sigma=0.0283;
mu=0.00868;
mu_A=0.00870;
mu_H=0.0507;
theta_2=0.148;
zeta=0.318;
theta_3=2.38;
nu=0.0482;
omega=0.0000000001;
b=0.283;
c=-0.0313;

%Baseline initial conditions 
P0=0.095;
A0=0.00647;
H0=0.000843;
R0=0.0584;
S0=1-P0-A0-H0-R0;

%Number of intervals to uniformly sample sigma 
N=100;

%Create vectors of parameters (with sigma varying) to run model with
m_vec=linspace(m,m,N);
beta_A_vec=linspace(beta_A,beta_A,N);
beta_P_vec=linspace(beta_P,beta_P,N);
theta_1_vec=linspace(theta_1,theta_1,N);
epsilon_vec=linspace(epsilon, epsilon,N);
gamma_vec=linspace(gamma,gamma,N);
sigma_vec=linspace(0.1*sigma,sigma,N);
mu_vec=linspace(mu,mu,N);
mu_A_vec=linspace(mu_A,mu_A,N);
mu_H_vec=linspace(mu_H,mu_H,N);
theta_2_vec=linspace(theta_2,theta_2,N);
zeta_vec=linspace(zeta,zeta,N);
theta_3_vec=linspace(theta_3,theta_3,N);
nu_vec=linspace(nu,nu,N);
omega_vec=linspace(omega,omega,N);
b_vec=linspace(b,b,N);
c_vec=linspace(c,c,N);
P0_vec=linspace(P0,P0,N);
A0_vec=linspace(A0,A0,N);
H0_vec=linspace(H0,H0,N);
R0_vec=linspace(R0,R0,N);

%Matrix of parameter sets
vec_matrix  = [m_vec;beta_A_vec;beta_P_vec;theta_1_vec;epsilon_vec;gamma_vec;sigma_vec;mu_vec;mu_A_vec;mu_H_vec;theta_2_vec;zeta_vec;theta_3_vec;nu_vec;omega_vec;b_vec;c_vec;P0_vec;A0_vec;H0_vec;R0_vec]';


%Run ODE with rows of vec_matrix 
for x=1:N;
    tspan=linspace(0,6,N);
    y0 = [S0,P0,A0,H0,R0]; 
    time_points = 6; 
    f=@ODE_LHS_Heroin;
    x;
    vec_matrix(x,:);
    
    [t,y] = ode15s(@(t,y)f(t,y,vec_matrix,x),tspan,y0,[]); 

    W = [t y];
    
    S_lhs(:,x)=W(time_points+1,2);
    P_lhs(:,x)=W(time_points+1,3);
    A_lhs(:,x)=W(time_points+1,4);
    H_lhs(:,x)=W(time_points+1,5);
    R_lhs(:,x)=W(time_points+1,6);
    
    
end  
    
 figure(1);
 plot((0.0283-vec_matrix(:,7))*100/0.0283,(0.004418807691222-A_lhs(1,:))*100/0.004418807691222,'o') 
 hold on
 plot((0.0283-vec_matrix(:,7))*100/0.0283,(0.002197350239427-H_lhs(1,:))*100/0.002197350239427,'o')
 hold off
 xlabel('Percent reduction in \sigma')
 ylabel('Percent reduction in A or H')
 legend({'Percent reduction in A', 'Percent reduction in H'},'FontSize', 10)
 set(gca,'XTick',0:10:100);
 xlim([0 90])
 ylim([0 85])
 
 
 