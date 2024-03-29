%File name: heroin_multistart_final_testing.m

%Parameters
%slope of alpha 
m=-.0332;
beta_A=0.000273; 
beta_P=0.000777; 
theta_1=0.0001;
epsilon=3.1419;
mu=0.00868; 
mu_A=0.00870;      
mu_H=0.0507;
gamma=0.0001;
theta_2=3*theta_1; 
sigma=0.0001;
zeta=0.0214;
theta_3=16*theta_1; 
nu=0.0155;
omega=0.0000000001;
%y-intercept of alpha 
b=0.6348;

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
N = 5;
tspan=linspace(0,N,N+1);

% Initial Conditions
P0=0.0710;
A0=0.0075;
H0=0.0017;
R0=0.4928;
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
      sum(i)=y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5);
  end
  
 
 % Yearly output from the model as a proportion of individuals in P at 
 % some point during the year for 2013-2017, Data1 is a column vector
 %Data1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 
 % Yearly output from the model as a proportion of individuals in A at 
 % some point during the year for 2013-2017, Data2 is a column vector
 %Data2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 

 % Yearly output from the model as a proportion of individuals in H at 
 % some point during the year for 2014-2016, Data3 is a column vector 
 %Data3=y(2:4,4)+y(3:5,8)-y(2:4,8);

 
 
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 %Data1=[0.399384466780766;0.476721593432771;0.469765087997695;0.449866817308604;0.427737148099884];
 
 % Actual Data for years 2013-2017
 Data1=[1825910./5517176; 1805325./5559006; 1800613./5602117; 1744766./5651993; 1620955./5708586];
 

 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 %Data2=[0.00709261474856600;0.0106675766257930;0.0130723402928730;0.0148284456654410;0.0162165492691030];
 
 % Actual Data for years 2013-2017 
 Data2=[43418./5517176; 42928./5559006; 42816./5602117; 37464./5651993; 34805./5708586];


 Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8); 
 %Data3=[0.00116527288223448;0.00120952017524577;0.00118883157707289];
 
 % Actual Data for years 2014-2016
 Data3=[7560./5559006; 7560./5602117; 10260./5651993];
 
 Diff1=Estim1-Data1;
 Diff2=Estim2-Data2;
 Diff3=Estim3-Data3;
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)

  
 % ODE solutions plotted separately shown all together
 figure(1)
         
           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           xlim([0 , N])
          
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           xlim([0 , N])
           
           subplot(2,2,3);plot(t,y(:,4) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           xlim([0 , N])
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Stably Recovered Individuals')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           xlim([0 , N])
           
          
               
                 
 % ODE Solutions for P, A, H plotted all together
 figure(2)
           plot(t,y(:,3),'r-','LineWidth',1);
           hold all
           plot(t,y(:,4),'g-','LineWidth',1); 
           hold all
           plot(t,y(:,5),'m-','LineWidth',1); 
           xlabel('Year')
           ylabel('Size of Populations');
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           legend('A','H','R')
           xlim([0 , N])
   
  % ODE Solutions for S and P plotted together
  figure(3) 
           subplot(2,2,1);plot(t,y(:,1),'r-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Susceptibles')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           xlim([0 , N])
           
           subplot(2,2,2);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Prescription Users')
           %set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           xlim([0 , N])
           
 figure(4)
           hold all
           plot(t,y(:,1))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Susceptible individuals')
           legend('Proportion of susceptible individuals')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           
 figure(5)
           hold all
           plot(t,y(:,2))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Prescription Users')
           legend('Proportion of prescription users')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           
           
 figure(6)
           hold all
           plot(t,y(:,3))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Prescription opioid addicts')
           legend('Proportion of prescription opioid addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           
                
 figure(7)
           hold all
           plot(t,y(:,4))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Heroin addicts')
           legend('Proportion of heroin addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           
      
 figure(8)
           hold all
           plot(t,y(:,5))
           set(gca, 'fontsize',10)
           xlabel('Year')
           ylabel('Stably recovered individuals')
           legend('Proportion of stably recovered individuals')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
           
                   
 % Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 plot(t(1:end-1),Estim1, 'o')
 plot(t(1:end-1), Data1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
 
 % Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 plot(t(1:end-1),Estim2, 'o')
 plot(t(1:end-1), Data2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [0 1 2 3 4])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017'})
 
 
 % Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(11)
 hold all
 plot(t(2:4),Estim3, 'o')
 plot(t(2:4), Data3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 
            
    

           
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




 