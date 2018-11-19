function value = HeroinModel_ODE45_Testing(z)

%%%GOING TO CHANGE FOR MY MODEL (right now it's a different model) 

% Final time 
N = 4;
T=N;
tspan=linspace(0,T,N+1);
global value 
 

% Estimated  values of parameters from "HeroinModel_MultiStart.m"
z0=[0.1  0.001  0.5   0.0001  0.8  0.001  0.001  0.002  0.0001  0.4  0.1  0.1  0.0001  0.7  0.05 0.1 0.01 0.01];

z=z0;

%Parameters
alpha=z(1); 
 
beta_A=z(2); 
 
beta_P=z(3);
 
theta_1=z(4);
 
epsilon=z(5);
 
mu=z(6);  
 
mu_A=z(7);   
 
mu_H=z(8);
 
gamma=z(9);   
 
theta_2=z(10); 
 
sigma_A=z(11);
 
zeta=z(12);
 
theta_3=z(13);
 
sigma_H=z(14);
 
nu=z(15);

A0=z(16);

H0=z(17);

R0=z(18);


%Initials
%MADE UP VALUES IN ORDER TO RUN CODE
S0=1-0.1-z(16)-z(17)-z(18); 
P0=0.1;
A0=z(16);
H0=z(17);
R0=z(18); 
X0=0;
Z0=0;
initials = [S0,P0,A0,H0,R0,X0,Z0];


[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  Z=y(:,7);
  
  
  % int of harvesting , h*A

 prescript=zeros(1,4);
 for i=1:4
 prescript(i) = y(i,2)+y(i+1,6)-y(i,6); 
 end
 
%yearly output from the model as a fraction
 Estim1=[prescript(1),prescript(2),prescript(3),prescript(4)];
       
% Proportions of population of prescription opioid users (MADE UP FRACTIONS
% FOR NOW)
 Data1=[.1 .2 .3 .25];
 

% the difference between estimated value and data 
 diff1= Estim1-Data1;
 
 
 
 %Using ODE y(7) in order to account for new admissions coming into recovery class
 %(NOT total in recovery class), going to run from 2013-2015, so want:
 admitted=zeros(1,2);
 for i=1:2
 admitted(i) = y(i+1,7)-y(i,7);
 end
 
 %yearly output from the model as a fraction
 Estim2=[admitted(1), admitted(2)];
  
 %Proportions of population being admitted into recovery (MADE UP FRACTIONS FOR NOW)
 Data2=[.001,.0015];

 % the difference between estimated value and data 
 diff2=Estim2-Data2;
 
 %Proportion of opioid addicts in 2015
 addicts(2)=y(2,3); 
 %yearly output from the model as a fraction
 Estim3=[addicts(2)];
%Made up for now
 Data3=[.4];
 %the difference between estimated value and data
 diff3=Estim3-Data3;

 
 %Proportion of heroin/fentanyl addicts in 2015
 heroin(2)=y(2,4);
 Estim4=[heroin(2)];
 %Made up for now
 Data4=[.2];
 %the difference between estimated value and data
 diff4=Estim4-Data4;
 
 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 %normalized by norm of the data
 value = norm(diff1,2)./norm(Data1)+ norm(diff2,2)./norm(Data2)+norm(diff3,2)./norm(Data3)+norm(diff4,2)./norm(Data4);
 
 
 figure(1)
 hold all
 scatter(1:1:15,Estim,'filled')
 scatter(1:1:15, Data, 'filled')
 set(gca, 'xtick', [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2002','2003','2004','2005','2006','2007','2008', '2009','2010', '2011', '2012','2013','2014','2015','2016'})
 ylabel('Landing')
 legend('Landing Simulated','Yearly Landing Data' )
 
 figure(2)

           subplot(2,2,1);plot(t,y(:,1),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Time')
           subplot(2,2,1);ylabel('Anchovy')
           xlim([0 , T])
           
           subplot(2,2,2);plot(t,y(:,2),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Time')
           subplot(2,2,2);ylabel('Jelly Fish')
           xlim([0 , T])
           
           subplot(2,2,3);plot(t,y(:,3) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Time')
           subplot(2,2,3);ylabel('Zoo-Plankton')
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
           subplot(2,2,4);
           plot(t,y(:,1),'b-','LineWidth',1);
           hold on
           plot(t,y(:,2),'r-','LineWidth',1);
           hold on
           plot(t,y(:,3),' g-','LineWidth',1); 
           xlabel('time')
           ylabel('Size of Populations');
           legend('A','P','Z')
           xlim([0 , T])
           ylim([0 , 2e+5])
           legend('A','P','Z')
 