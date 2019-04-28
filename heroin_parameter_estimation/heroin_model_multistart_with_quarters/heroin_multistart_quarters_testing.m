%File name: heroin_multistart_final_testing.m

%Parameters
%slope of alpha 
m=-.0123;
beta_A=0.000273; 
beta_P=0.000777; 
theta_1=0.00001;
epsilon=3;
mu=0.00868; 
mu_A=0.00870;      
mu_H=0.0507;
gamma=0.00744;
theta_2=3*theta_1; 
sigma=0.00744;
zeta=0.0214;
theta_3=16*theta_1; 
nu=0.0155;
omega=0.0000000001;
%y-intercept of alpha 
b=0.3; 

%{
% For R_0 checking:
alpha=0.2; 
beta_A=0.000273; 
beta_P=0; 
theta_1=0.0003;
epsilon=1.5;
mu=0.00868; 
mu_A=0.00775;   
mu_H=0.0271;
gamma=0;   
theta_2=3*theta_1; 
sigma=0.7;
zeta=0.25;
theta_3=16*theta_1; 
nu=0.1;
omega=0.0000000001;
%}

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Final time and N+# is # of equally spaced points from 0 to N 
N = 24;
tspan=linspace(0,N,N+1);

% Initial Conditions
P0=0.2;
A0=0.00760;
H0=0.00121;
R0=0.00443;
S0=1-P0-A0-H0-R0;
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
  
  alpha=m*t+b;
  
  % Making sure S+P+A+H+R=1
  for i=1:N+1
      total(i)=y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5);
  end
  
 
Estim1=[sum(y(1:4,2)+y(2:5,6)-y(1:4,6)); sum(y(5:8,2)+y(6:9,6)-y(5:8,6)); sum(y(9:12,2)+y(10:13,6)-y(9:12,6));...
         sum(y(13:16,2)+y(14:17,6)-y(13:16,6)); sum(y(17:20,2)+y(18:21,6)-y(17:20,6))];
 

Estim2=[sum(y(1:4,3)+y(2:5,7)-y(1:4,7)); sum(y(5:8,3)+y(6:9,7)-y(5:8,7)); sum(y(9:12,3)+y(10:13,7)-y(9:12,7));...
         sum(y(13:16,3)+y(14:17,7)-y(13:16,7)); sum(y(17:20,3)+y(18:21,7)-y(17:20,7))];
    
     
Estim3=[sum(y(5:8,4)+y(6:9,8)-y(5:8,8)); sum(y(9:12,4)+y(10:13,8)-y(9:12,8));...
         sum(y(13:16,4)+y(14:17,8)-y(13:16,8))];
 
 

 
Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);

 
 % ODE solutions plotted separately shown all together
 figure(1)
         
           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           set(gca, 'xtick', [ 0 4 8 12 16 20 24 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018','2019'})
 
 
          
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           set(gca, 'xtick', [ 0 4 8 12 16 20 24 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})
 
           
           subplot(2,2,3);plot(t,y(:,4) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           set(gca, 'xtick', [ 0 4 8 12 16 20 24 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})
 
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Stably Recovered Individuals')
           set(gca, 'xtick', [ 0 4 8 12 16 20 24 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})  
               
                 
 % ODE Solutions for P, A, H plotted all together
 figure(2)
           plot(t,y(:,3),'r-','LineWidth',1);
           hold all
           plot(t,y(:,4),'g-','LineWidth',1); 
           hold all
           plot(t,y(:,5),'m-','LineWidth',1); 
           xlabel('Year')
           ylabel('Size of Populations');
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
           set(gca, 'fontsize',10)
           xtickangle(90)
           legend('A','H','R')
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 

 
 figure(3)
           hold all
           plot(t,y(:,1))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Susceptible individuals')
           legend('Proportion of susceptible individuals')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
           
 figure(4)
           hold all
           plot(t,y(:,2))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Prescription Users')
           legend('Proportion of prescription users')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
           
           
 figure(5)
           hold all
           plot(t,y(:,3))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Prescription opioid addicts')
           legend('Proportion of prescription opioid addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
                
 figure(6)
           hold all
           plot(t,y(:,4))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Heroin addicts')
           legend('Proportion of heroin addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
      
 figure(7)
           hold all
           plot(t,y(:,5))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Stably recovered individuals')
           legend('Proportion of stably recovered individuals')
          set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
              
           
 figure(8)
 hold all
 scatter(0:1:4, Estim1,'o')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('ODE solution')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
 

 
 % Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 scatter(0:1:4, Estim2,'o')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('ODE solution')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017'})
 



 % Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 scatter(1:1:3, Estim3,'o')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('ODE solution')
 set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 

 figure(11)
 hold all
 plot(t(1:end-1),Estim4, 'o')
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in P each quarter')
 legend('ODE solution')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 23])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
 
           
function f = HeroinModel(t,y,pars)
f=zeros(8,1);
f(1)=-(pars(1)*t+pars(16))*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=(pars(1)*t+pars(16))*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);


% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
f(6) = (pars(1)*t+pars(16))*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));


end




 