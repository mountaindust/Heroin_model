% File name: HeroinModel_ODE45_Testing.m (used to be in Heroin_model folder)

% Later: if plotting Estim/Data points, put in: 
% function value = HeroinModel_ODE45_Testing(z)


% Final time 
N = 25;
T = N;
tspan=linspace(0,T,N+1);
% Later, if need at all with Estim/Data points: global value 
 

% Estimated parameter values from "HeroinModel_MultiStart.m"
% z0=[0.15 0.00094 0.00266 0.0001 3.25 0.00744 0.0002 0.5 0.05 0.0004 0.05 0.1 0.0057 0.0013 0.009];

% For testing, selected parameter values based on opioid paper/something realistic, 
% then simulated data, put data into HeroinModel_ODE45.m file Data1, Data2, Data3 vectors, 
% then ran HeroinModel_MultiStart.m to see if got back these parameter values
z0=[0.15 0.00094 0.00266 0.0001 3.25 0.00744 0.0002 0.5 0.05 0.0004 0.05 0.1 0.0057 0.0013 0.009];

z=z0;

%Parameters
alpha=z(1); 

beta_A=z(2); 
 
beta_P=z(3);
 
theta_1=z(4);
 
epsilon=z(5);
 
mu=0.00868;  
 
mu_A=0.00775;   
 
mu_H=0.0271;
 
gamma=z(6);   
 
theta_2=z(7); 
 
sigma=z(8);
 
zeta=z(9);
 
theta_3=z(10);
 
nu=z(11);

P0=z(12);

A0=z(13);

H0=z(14);

R0=z(15);

omega=.00001;

% Initial conditions 
S0=1-z(12)-z(13)-z(14)-z(15); 
P0=z(12);
A0=z(13);
H0=z(14);
R0=z(15); 
X0=0;
L0=0;
M0=0;
initials = [S0,P0,A0,H0,R0,X0,L0,M0];



[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);
  
  
 % Later: if plotting Estim/Data points, copy/paste "%%COMPARING MODEL
 % ESTIMATES TO DATA%%" from HeroinModel_ODE45.m file.
 % Note: if want to display Estim# points explicitly in command window, 
 % write the following in the code, for example:
 % Estim1 

 % ODE solutions plotted separately 
 figure(1)

           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           
           subplot(2,2,3);plot(t,y(:,4) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
          
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Recovered Individuals')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
               
                 
 % ODE Solutions all plotted together
 figure(2)
           plot(t,y(:,2),'b-','LineWidth',1);
           hold all
           plot(t,y(:,3),'r-','LineWidth',1);
           hold all
           plot(t,y(:,4),'g-','LineWidth',1); 
           hold all
           plot(t,y(:,5),'m-','LineWidth',1); 
           xlabel('time')
           ylabel('Size of Populations');
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           legend('P','A','H','R')
           xlim([0 , T])
           ylim([0 , 0.1])
           legend('P','A','H','R')
           
           
%%% Later: plots with Estim/Data points if needed 
 %{
 %Data points from Estim1, Data1 (in future: plot model on top?) 
 figure(3)
 hold all
 scatter(1:1:5, Estim1,'filled')
 scatter(1:1:5, Data1, 'filled')
 set(gca, 'xtick', [1 2 3 4 5])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017'})
 xlabel('Year')
 ylabel('Proportion of prescription users')
 legend('Proportion in prescription users simulated','Proportion of prescription users data')

 
 %Data points from Estim2, Data2 (in future: plot model on top?) 
 figure(4)
 hold all
 scatter(1:1:1, Estim2,'filled')
 scatter(1:1:1, Data2, 'filled')
 set(gca, 'xtick', [ 1 2 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2015' '2016'})
 xlabel('Year')
 ylabel('Proportion of opioid addicts')
 legend('Proportion of opioid addicts simulated','Proportion of opioid addicts data' )
 
  
 %Data points from Estim3, Data3 (in future: plot model on top?) 
 figure(5)
 hold all
 scatter(1:1:1, Estim3,'filled')
 scatter(1:1:1, Data3, 'filled')
 set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 xlabel('Year')
 ylabel('Proportion of heroin/fentanyl users')
 legend('Proportion of heroin/fentanyl users simulated','Proportion of heroin/fentanyl users data' )
 

 %Plots of ODES used in Estim1-Estim3
 figure(6)
          
           subplot(2,2,1);plot(t,y(:,6),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Proportion entering P')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           
           
           subplot(2,2,2);plot(t,y(:,7) ,' y-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Proportion entering A')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           
           
          
           subplot(2,2,3);plot(t,y(:,8) ,' m-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Proportion entering H')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
    
           
           %}


 
 
 