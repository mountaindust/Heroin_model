%File name: heroin_multistart_final_testing.m

%Parameters
%alpha=0.3; 
%beta_A=0.0094; 
%beta_P=0.00266; 
%theta_1=0.0003;
%epsilon=2.0;
%mu=0.00868; 
%mu_A=0.00775;   
%mu_H=0.0271;
%gamma=0.00744;   
%theta_2=0.0005; 
%sigma=0.7;
%zeta=0.1;
%theta_3=0.005; 
%nu=0.05;
%omega=0.0000000001;

%alpha=0.2; 
%beta_A=0.00094; 
%beta_P=0.00266; 
%theta_1=0.0003;
%epsilon=1.5;
%mu=0.00868; 
%mu_A=0.00775;   
%mu_H=0.0271;
%gamma=0.00744;   
%theta_2=0.0006; 
%sigma=0.7;
%zeta=0.25;
%theta_3=0.0009; 
%nu=0.1;
%omega=0.0000000001;



alpha=0.2; 
beta_A=0.00094; 
beta_P=0.00266; 
theta_1=0.0003;
epsilon=1.5;
mu=0.00868; 
mu_A=0.00775;   
mu_H=0.0271;
gamma=0.00744;   
theta_2=3*theta_1; 
sigma=0.7;
zeta=0.25;
theta_3=16*theta_1; 
nu=0.1;
omega=0.0000000001;

pars=[alpha,beta_A,beta_P,theta_1,epsilon,0.00868,0.00775,0.0271,gamma,theta_2,sigma,zeta,theta_3,nu,0.0000000001];

% Final time and N+# is # of equally spaced points from 0 to N 
N = 5;
tspan=linspace(0,N,N+1);

% Initial conditions 
%S0=1-0.13-0.01-0.001-0.0003; 
%P0=0.13;
%A0=0.01;
%H0=0.001;
%R0=0.0003; 
%X0=0;
%L0=0;
%M0=0;
%initials = [S0;P0;A0;H0;R0;X0;L0;M0];

S0=1-0.0538-0.0022-0.00074-0.000091;
P0=0.0538;
A0=0.0022;
H0=0.00074;
R0=0.000091;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];



[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);
  
  
  % Making sure S+P+A+H+R=1
  for i=1:N+1
      sum(i)=y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5);
  end
  
 
 % Yearly output from the model as a proportion of individuals in P at 
 % some point during the yearfor 2013-final year, Data1 is a column vector
% Data1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  
 %for 2013-2017
 Data1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 
 % Yearly output from the model as a proportion of individuals in A at 
 % some point during the yearfor 2013-final year, Data2 is a column vector
% Data2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7);  
 %for 2014 and 2015 
 Data2=y(2:3,3)+y(3:4,7)-y(2:3,7); 

    
 % Yearly output from the model as a proportion of individuals in H at 
 % some point during the yearfor 2013-final year, Data3 is a column vector 
% Data3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8); 
 %for 2014 and 2015 
 Data3=y(2:4,4)+y(3:5,8)-y(2:4,8);

 
  
 % ODE solutions plotted separately 
 figure(1)
         
           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
          
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
           
           subplot(2,2,3);plot(t,y(:,4) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Stably Recovered Individuals')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
           
          
               
                 
 % ODE Solutions plotted all together
 figure(2)
           plot(t,y(:,3),'r-','LineWidth',1);
           hold all
           plot(t,y(:,4),'g-','LineWidth',1); 
           hold all
           plot(t,y(:,5),'m-','LineWidth',1); 
           xlabel('time')
           ylabel('Size of Populations');
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           legend('A','H','R')
           xlim([0 , N])
       
  figure(3) 
    subplot(2,2,1);plot(t,y(:,1),'r-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Susceptibles')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
           
           subplot(2,2,2);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Prescription Users')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           %set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
           
 %{      
 %Data points from Data1 and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 scatter(t(1:end-1), Data1)
 %plot(t(1:end-1), Data1)
 %set(gca, 'xtick', [1 2 3 4 5])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2013','2014','2015','2016','2017'})
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('Data points interested in for P')

 
 %Data points from Data2 and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 scatter(t(1:end-1), Data2)
 %plot(t(1:end-1), Data2)
 %set(gca, 'xtick', [ 1 2 ])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2015' '2016'})
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('Data points interested in for A')

  
 % Data points from Data3 and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 scatter(t(1:end-1), Data3)
 %plot(t(1:end-1), Data3)
 %set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2014', '2015', '2016'})
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('Data points interested in for H')
  %}

           
function f = HeroinModel(t,y,pars)
f=zeros(8,1);
f(1)=-pars(1)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=pars(1)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);


% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
f(6) = pars(1)*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));


end





 