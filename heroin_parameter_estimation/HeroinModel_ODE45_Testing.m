% File name: HeroinModel_ODE45_Testing.m (used to be in Heroin_model folder)

% Later: if plotting Estim/Data points, put in:
%function value = HeroinModel_ODE45_Testing(z)


% Final time 
N = 25;
T = N;
tspan=linspace(0,T,N+1);
% Later, if need at all with Estim/Data points: global value 
 

% Estimated parameter values from "HeroinModel_MultiStart.m"
%z =[alpha beta_A  beta_P   theta_1 epsilon  gamma   theta_2  sigma zeta   theta_3   nu  P0    A0    H0     R0  ]
z0=[0.15  0.00094  0.00266   0.0001   2.5   0.00744   0.0002   0.5  0.05   0.0004   0.05 0.1 0.0057 0.0013 0.009];

% For testing, selected parameter values based on opioid paper/something realistic, 
% then simulated data, put data into HeroinModel_ODE45.m file Data1, Data2, Data3 vectors, 
% then ran HeroinModel_MultiStart.m to see if got back these parameter values
%z0=[0.2 0.3 0.2 0.2 0.85 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.2 0.1 0.1 0.3];

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


  S=y(:,1)';
  P=y(:,2)';
  A=y(:,3)';
  H=y(:,4)';
  R=y(:,5)';
  X=y(:,6)';
  L=y(:,7)';
  M=y(:,8)';
  
  
 % Later: if plotting Estim/Data points, copy/paste "%%COMPARING MODEL
 % ESTIMATES TO DATA%%" from HeroinModel_ODE45.m file.
 % Note: if want to display Estim# points explicitly in command window, 
 % write the following in the code, for example:
 % Estim1 
 
 
 % Yearly output from the model as a proportion of P individuals for
 % 2013-final year, Estim1 is a row vector
 Data1=zeros(1,25);
 % For 2013:
 Data1(1)=P0+y(1,6);  
 % For 2014-final year:
    for i=2:25
       Data1(i)= y(i-1,2)+y(i,6)-y(i-1,6);
    end
    
 
 % Yearly output from the model as a proportion of A individuals for
 % 2013-final year, Estim2 is a row vector
 Data2=zeros(1,25);  
 % For 2013:
 Data2(1)=A0+y(1,7);  
 % For 2014-final year:
    for i=2:25
       Data2(i)= y(i-1,3)+y(i,7)-y(i-1,7);
    end
    
    
 % Yearly output from the model as a proportion of H individuals for
 % 2013-final year, Estim3 is a row vector 
 Data3=zeros(1,25);  
 % For 2013:
 Data3(1)=H0+y(1,8);  
 % For 2014-final year:
    for i=2:25
       Data3(i)= y(i-1,4)+y(i,8)-y(i-1,8);
    end
   
 
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
           
           
 
 %Data points from Estim1
 figure(3)
 hold all
 %scatter(1:1:25, Estim1,'filled')
 scatter(1:1:N, Data1, 'filled')
 %set(gca, 'xtick', [1 2 3 4 5])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2013','2014','2015','2016','2017'})
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 %legend('Proportion in prescription users simulated','Proportion of prescription users data')

 
 %Data points from Estim2
 figure(4)
 hold all
 %scatter(1:1:25, Estim2,'filled')
 scatter(1:1:N, Data2, 'filled')
 %set(gca, 'xtick', [ 1 2 ])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2015' '2016'})
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 %legend('Proportion of opioid addicts simulated','Proportion of opioid addicts data' )
 
  
 %Data points from Estim3, Data3 (in future: plot model on top?) 
 figure(5)
 hold all
 %scatter(1:1:25, Estim3,'filled')
 scatter(1:1:N, Data3, 'filled')
% set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2014', '2015', '2016'})
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 %legend('Proportion of heroin/fentanyl users simulated','Proportion of heroin/fentanyl users data' )
 
 
 
 %%% Later: plots with extra ODEs if needed 
%{
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


 
 
 