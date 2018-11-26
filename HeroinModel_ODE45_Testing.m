function value = HeroinModel_ODE45_Testing(z)

%%%GOING TO CHANGE FOR MY MODEL (right now it's a different model) 

% Final time 
N = 4;
T=N;
tspan=linspace(0,T,N+1);
global value 
 

% Estimated  values of parameters from "HeroinModel_MultiStart.m"
z0=[0.1  0.001  0.5   0.0001  0.8  0.0001  0.4  0.1  0.1  0.0001  0.7  0.05 0.1 0.01 0.01 0.01];

z=z0;

%Parameters
alpha=z(1); 
 
beta_A=z(2); 
 
beta_P=z(3);
 
theta_1=z(4);
 
epsilon=z(5);
 
mu=0.00864;  
 
mu_A=0.00775;   
 
mu_H=0.0271;
 
gamma=z(6);   
 
theta_2=z(7); 
 
sigma_A=z(8);
 
zeta=z(9);
 
theta_3=z(10);
 
sigma_H=z(11);
 
nu=z(12);

P0=z(13);

A0=z(14);

H0=z(15);

R0=z(16);



%Initials
%MADE UP VALUES IN ORDER TO RUN CODE
S0=1-z(13)-z(14)-z(15)-z(16); 
P0=z(13);
A0=z(14);
H0=z(15);
R0=z(16); 
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
  
  
 %COMPARING MODEL ESTIMATES TO DATA 
 
 %%%%%
 % In order to count the total number of individuals in P at some point throughout a certain year
 % (the number who are in the class AT ALL during the year, it's okay if they leave), 
 % we want to add the total number of individuals using prescriptions at the beginning of the
 % year to those who come into P at some point during the year: 
 % prescript(i) is number of prescription users at beginning of the year (time step) +
 % the number of new cases that came in from X'=dy(6) ODE from that year until the beginning 
 % of the next year. Thus, we are adding in those who come into P at some point during 
 % the year by integrating y(6) ODE but just focusing in on the one year we care about 
 %(so have to subtract: integrating gives total number of new cases from t=0 to t=i, so have to 
 % subtract off the number from t=0 to t=i-1). 
 % Here, calculating for years 2014-2017: 
 
 total_prescription_users=zeros(1,4);
 for i=1:4
 total_prescription_users(i) = y(i,2)+y(i+1,6)-y(i,6); 
 end
 
% yearly output from the model as a proportion from 2014 to 2017, excludes
%2013
 Estim1=[total_prescription_users(1),total_prescription_users(2),total_prescription_users(3), total_prescription_users(4)];
       
% actual proportions of population that were prescription opioid users for 2014-2017 (MADE UP FRACTIONS
% FOR NOW)
 Data1=[.1 .2 .3 .25];
 
% the difference between estimated value and data 
 diff1= Estim1-Data1;
 
%%%%%
% We have total number of individuals who take prescription opioids for the
% year 2013; so must take IC (which is estimated because don't know number
% at beginning of the year) and add on the number of individuals that enter
% the P class at any point during the year 2013, which comes from
% integrating ODE X'=dy(6) from t=0 to t=1

% output from the model as a proportion in 2013, MUST PUT IN z(13) value
% here=0.1!!
 Estim2=[0.1+y(1,6)];
       
% actual proportion of population that were prescription opioid users in 2013 (MADE UP FRACTIONS
% FOR NOW)
 Data2=[.2];
 

% the difference between estimated value and data 
 diff2= Estim2-Data2;
 
 
%%%%%
 %To calculate number of new admissions coming into the recovery class
 %(NOT total in recovery class), we use ODE Z'=dy(7); going to run from
 %2013-2015 because those are the only years we have data for:
 new_admissions=zeros(1,2);
 for i=1:2
 new_admissions(i) = y(i+1,7)-y(i,7);
 end
 
 % yearly output from the model as a proportion in the recovery class
 Estim3=[new_admissions(1), new_admissions(2)];
  
 % actual proportions of population being admitted into recovery (MADE UP FRACTIONS FOR NOW)
 Data3=[.001 .0015];

 % the difference between estimated value and data 
 diff3=Estim3-Data3;
 
 
 %%%%%
 % output from the model of the proportion of opioid addicts in 2015
 Estim4=[y(2,3)];
 % made up for now
 Data4=[.4];
 % the difference between estimated value and data
 diff4=Estim4-Data4;
 

 %%%%%
 % output from the model of the proportion of heroin/fentanyl addicts in 2015
 Estim5=[y(2,4)];
 % Made up for now
 Data5=[.2];
 % the difference between estimated value and data
 diff5=Estim5-Data5;
 
 
 %%%%%
 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 %normalized by norm of the data
 value = norm(diff1,2)./norm(Data1)+ norm(diff2,2)./norm(Data2)+...
 norm(diff3,2)./norm(Data3)+norm(diff4,2)./norm(Data4)+norm(diff5,2)./norm(Data5);
 
 
 figure(1)
 hold all
 scatter(1:1:4, Estim1,'filled')
 scatter(1:1:4, Data1, 'filled')
 set(gca, 'xtick', [ 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014','2015','2016','2017'})
 xlabel('Year')
 ylabel('Proportion of prescription users')
 legend('Prescription users simulated','Prescription users data' )
 
 
 figure(2)
 hold all
 scatter(1:1:1, Estim2,'filled')
 scatter(1:1:1, Data2, 'filled')
 set(gca, 'xtick', [ 1 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013'})
 xlabel('Year')
 ylabel('Proportion of prescription users')
 legend('Prescription users simulated','Prescription users data' )
 
 
  
 figure(3)
 hold all
 scatter(1:1:1, Estim4,'filled')
 scatter(1:1:1, Data4, 'filled')
 set(gca, 'xtick', [ 1 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2015'})
 xlabel('Year')
 ylabel('Proportion of opioid addicts')
 legend('Opioid addicts simulated','Opioid addicts data' )
 
  
 figure(4)
 hold all
 scatter(1:1:1, Estim5,'filled')
 scatter(1:1:1, Data5, 'filled')
 set(gca, 'xtick', [ 1 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2015'})
 xlabel('Year')
 ylabel('Proportion of heroin/fentanyl users')
 legend('Heroin/fentanyl users simulated','Heroin/fentanyl users data' )
 
 
 
 figure(5)

           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           xlim([0 , T])
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           xlim([0 , T])
           
           subplot(2,2,3);plot(t,y(:,4) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Individuals in recovery')
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
           %All together, later make separate figure: 
           %subplot(2,2,4);
           %plot(t,y(:,1),'b-','LineWidth',1);
           %hold on
           %plot(t,y(:,2),'r-','LineWidth',1);
           %hold on
           %plot(t,y(:,3),' g-','LineWidth',1); 
           %xlabel('time')
           %ylabel('Size of Populations');
           %legend('A','P','Z')
           %xlim([0 , T])
           %ylim([0 , 2e+5])
           %legend('A','P','Z')
 