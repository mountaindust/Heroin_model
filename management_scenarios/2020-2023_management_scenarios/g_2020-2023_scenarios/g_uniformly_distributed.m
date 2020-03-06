
Parameter_settings_vec_Heroin;

%Number of intervals to uniformly sample sigma 
N=1000;

%Create vectors of parameters (with sigma varying) to run model with
beta_A_vec=linspace(beta_A,beta_A,N);
beta_P_vec=linspace(beta_P,beta_P,N);
theta_1_vec=linspace(theta_1,theta_1,N);
epsilon_vec=linspace(epsilon,epsilon,N);
gamma_vec=linspace(gamma,gamma,N);
sigma_vec=linspace(sigma,sigma,N);
mu_vec=linspace(mu,mu,N);
mu_H_vec=linspace(mu_H,mu_H,N);
theta_2_vec=linspace(theta_2,theta_2,N);
zeta_vec=linspace(zeta,zeta,N);
theta_3_vec=linspace(theta_3,theta_3,N);
nu_vec=linspace(nu,nu,N);
omega_vec=linspace(omega,omega,N);
g_vec=linspace(-0.038616505394549,g,N); %50% alpha value at 2023 compared to baseline case requires g to go to -0.038616505394549
h_vec=linspace(h,h,N); %to make muA decrease and go toward 0 like mu_H
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
    
 
%whatever entry of each vector I am interested in (anywhere from 1 to
 %1000) 
 entry=500
 
 disp('percent decrease of sigma')
 (0.101518004918260-vec_matrix(entry,6))*100/0.101518004918260
 disp('value of sigma')
 sigma_vec(1,entry)
 
 disp('value of A with this sigma value')
 A_lhs(1,entry)
 disp('percent decrease of A from baseline 2023 value with this new sigma value')
 (0.001681716135399-A_lhs(1,entry))*100/0.001681716135399
 
 disp('value of H with this new sigma value')
 H_lhs(1,entry)
 disp('percent decrease of H from baseline 2023 value with this new sigma value')
 (0.013798429999561-H_lhs(1,entry))*100/0.013798429999561
 

 
 
