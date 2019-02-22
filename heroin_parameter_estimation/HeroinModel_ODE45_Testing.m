% File name: HeroinModel_ODE45_Testing.m (used to be in Heroin_model folder)

% Later: if want function value put in:
%function value = HeroinModel_ODE45_Testing(z)


% Final time 
N = 25;
% # of equally spaced points from 0 to N 
tspan=linspace(0,N,N+1);

% Later: if want function value put in: 
%global value 
 
% Estimated parameter values from "HeroinModel_MultiStart.m" (or parameter
% values used for testing)
%z =[alpha beta_A  beta_P   theta_1 epsilon  gamma   theta_2  sigma  zeta   theta_3   nu     P0       A0      H0     R0  ]

%z0=[0.2  0.00094  0.00266   0.0003    1.5    0.00744   0.0006  0.7   0.25    0.009   0.1    0.05   0.0062   0.0026  0.0006];
%z0=[0.199898262823408,0.00119461029644510,0.00575064364967549,0.000713076997491307,1.50066832051824,0.00466755968403883,0.00550270808269071,0.710660286077390,0.247820314304348,0.0299600600512476,0.0999524189000383,0.0501091691248067,0.00620843889669422,0.00260306913866477,0.000552026009579225];

%z0=[0.4 0.002 0.007 0.06 2.2 0.03 0.1 0.99 0.09 0.2 0.3 0.2 0.03 0.0003 0.0001];
%z0=[0.399847330011085,0.00185318387071895,0.00612436772706274,0.0930041222760489,2.19597306249974,0.0306124131142860,0.250862360061604,1.00022276811523,0.108720635939099,0.240070511978861,0.367229181276194,0.200101464894243,0.0294703373972747,0.000283848333990047,0.000588889097266240];

%z0=[0.4 0.002 0.007 0.06 2.2 0.03 0.1 0.99 0.09 0.2 0.3 0.3 0.03 0.0003 0.001];
%z0=[0.399742161527240,0.00178358977815605,0.00549114597752625,0.0826408509368144,2.19557128379611,0.0310749089090507,0.193609088002418,1.03204397859285,0.0996508670400122,0.233228789911452,0.342423492635243,0.300266107091212,0.0297862607645851,0.000289534406510186,0.00109989310808790];

z0=[0.3 0.00094 0.00266 0.0003 2 0.00744 0.0005 0.7 0.1 0.005 0.05 0.13 0.01 0.001 0.0003];
%z0=[0.299399217572142,0.000656061593271862,0.00542021884103976,0.0481875542428343,1.98402418305403,0.00459277587485556,0.247136325502546,0.692239069236622,0.111094641850517,0.194567054510244,0.139621617363135,0.130785262291572,0.00991135366170681,0.000922012938305255,0.000523903325876204];

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

omega=0.0000000001;

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
  
  
  % Making sure S+P+A+H+R=1
  for i=1:N+1
      sum(i)=y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5);
  end
  
 

 % Yearly output from the model as a proportion of P individuals for
 % 2013-final year, Data1 is a column vector
 Data1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  
 Data1=Data1';
% Same as Data1
% Check1=zeros(1,length(y)-1);  
 % For 2013-year N-1:
    %for i=1:length(y)-1
     %  Check1(i)= y(i,2)+y(i+1,6)-y(i,6);
  %  end 
    
    
 % Yearly output from the model as a proportion of A individuals for
 % 2013-final year, Data2 is a column vector
 Data2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7);  
 Data2=Data2';
    
 %Same as Data2
 % Check2=zeros(1,length(y)-1);  
 % For 2013-year N-1:
    %for i=1:length(y)-1
      % Check2(i)= y(i,3)+y(i+1,7)-y(i,7);
   % end   
    
    
 % Yearly output from the model as a proportion of H individuals for
 % 2013-final year, Data3 is a column vector 
 Data3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);   
 Data3=Data3';
 % For 2013 to year N-1:
    
   
 % Same as Data3
 % Check3=zeros(1,length(y)-1);  
 % For 2013-year N-1:
   % for i=1:length(y)-1
     % Check3(i)= y(i,4)+y(i+1,8)-y(i,8);
    %end   
    
    
 % ODE solutions plotted separately 
 figure(1)
         
           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
          % set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
          % set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
          % set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
         %  set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
           
           subplot(2,2,3);plot(t,y(:,4) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
          % set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
          % set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Recovered Individuals')
          % set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
          % set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , N])
               
                 
 % ODE Solutions all plotted together
 figure(2)
           plot(t,y(:,3),'r-','LineWidth',1);
           hold all
           plot(t,y(:,4),'g-','LineWidth',1); 
           hold all
           plot(t,y(:,5),'m-','LineWidth',1); 
           xlabel('time')
           ylabel('Size of Populations');
          % set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
          % set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           legend('A','H','R')
           xlim([0 , N])
           legend('A','H','R')
           
          
 %Data points from Data1 and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 %scatter(1:1:25, Estim1,'filled')
 %scatter(1:1:N, Data1, 'filled')
 %plot(t,y(:,6),'b-','LineWidth',1)
 scatter(t(1:end-1), Data1)
 %plot(t(1:end-1), Data1)
 plot(t(1:end-1), Data1)
 %plot(t(1:end-1), Check1)
 %set(gca, 'xtick', [1 2 3 4 5])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2013','2014','2015','2016','2017'})
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('Data points interested in', 'ODE solution')
 %legend('Proportion in prescription users simulated','Proportion of prescription users data')

 
 %Data points from Data2 and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 %scatter(1:1:N, Estim2,'filled')
 %scatter(0:1:N-1, Data2, 'filled')
 %plot(t(1:end-1), Check2)
 scatter(t(1:end-1), Data2)
 %plot(t(1:end-1), Data2)
 plot(t(1:end-1), Data2)
 %plot(t(1:end-1), Check2)
 %set(gca, 'xtick', [ 1 2 ])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2015' '2016'})
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('Data points interested in', 'ODE solution')
 %legend('Proportion of opioid addicts simulated','Proportion of opioid addicts data' )
 
  
 % Data points from Data3 and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 scatter(t(1:end-1), Data3)
 %scatter(0:1:N-1, Data3)
 plot(t(1:end-1), Data3)
 %plot(t(1:end-1), Check3)
 %set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 %set(gca,'xticklabel',{'2014', '2015', '2016'})
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('Data points interested in', 'ODE solution')
% legend('Proportion of heroin/fentanyl users simulated','Proportion of heroin/fentanyl users data' )
 



 % Note: if want to display Data# points explicitly in command window, 
 % write the following in the code, for example: Data1 
 
 