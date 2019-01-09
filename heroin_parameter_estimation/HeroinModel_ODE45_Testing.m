%File name: HeroinModel_ODE45_Testing.m
function value = HeroinModel_ODE45_Testing(z)


% Final time 
N = 4;
T = N;
tspan=linspace(0,T,N+1);
global value 
 

% Estimated  values of parameters from "HeroinModel_MultiStart.m"
% For testing, selected values based on opioid paper or something realistic (then
% simulated data, put in HeroinModel_ODE45.m file, then ran
% HeroinModel_MultiStart.m to see if got back these values for z
z0=[0.15 0.00094 0.00266 0.0001 3.25 0.00744 0.0002 0.5 0.05 0.0004 0.8 0.05 0.1 0.0057 0.0013 0.009];
% When tested got:
%z0=[[0.150161107255521,0.00380869386967680,0.00149298769055621,0.0872984749985357,3.25439679479616,0.00151684474187485,0.104431997965411,0.432381864908839,0.0504849321396042,0.304161057476428,0.460884286435396,0.0519159523800158,0.0999954889122865,0.00589489064530614,0.00189961986359587,0.00908242522815079]
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
S0=1-z(13)-z(14)-z(15)-z(16); 
P0=z(13);
A0=z(14);
H0=z(15);
R0=z(16); 
X0=0;
Z0=0;
K0=0;
L0=0;
M0=0;
initials = [S0,P0,A0,H0,R0,X0,Z0,K0,L0,M0];



[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  Z=y(:,7);
  K=y(:,8);
  L=y(:,9);
  M=y(:,10);
  
  
  
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
 
% We have total number of individuals who take prescription opioids for the
% year 2013; so must take IC (which is estimated because don't know number
% at beginning of the year) and add on the number of individuals that enter
% the P class at any point during the year 2013, which comes from
% integrating ODE X'=dy(6) from t=0 to t=1; this gives first value in
% Estim1

 
% yearly output from the model as a proportion from 2013 to 2017
 Estim1=[z(13)+y(1,6),total_prescription_users(1),total_prescription_users(2),total_prescription_users(3), total_prescription_users(4)];
       
% actual proportions of population that were prescription opioid users for
% 2013-2017 (total number of prescription opioid users in each year in TN that are 12 and older divided by
% total population in TN 12 and older) 
 Data1=[1845144./5517176 1824342./5559006 1819581./5602117 1761363./5651993 1636374./5708586];
%Data simulated when put in 0.1 for all parameters
%Data1=[0.1 0.16 0.18 0.2 0.21];
% the difference between estimated value and data 
 diff1= Estim1-Data1;
 
%print Estim1 values 
 Estim1
 
%%%%%
 %To calculate number of new admissions coming into the recovery class 
 %(NOT total in recovery class) from the opioid addict class, we use ODE Z'=dy(7); going to run from
 %2013-2015 because those are the only years among 2013-2017 we have data for; 
 %new_opioid_admissions is for years 2014-2015, and y(1,7) is the number that are
 %admitted in the first year 2013
 new_opioid_admissions=zeros(1,2);
 for i=1:2
 new_opioid_admissions(i) = y(i+1,7)-y(i,7);
 end
 
 % yearly output from the model as a proportion in the recovery class
 Estim2=[y(1,7), new_opioid_admissions(1), new_opioid_admissions(2)];
  
 % actual proportions of population each year being admitted into recovery from opioid addict class
 Data2=[4485./5517176 4530./5559006 4326./5602117];
 %%Data simulated when put in 0.1 for all parameters
 %Data2=[0 0.011 0.013];

 % the difference between estimated value and data 
 diff2=Estim2-Data2;
 
 %print Estim2 values 
 Estim2
 
 
 
 %%%%%
 %To calculate number of new admissions coming into the recovery class 
 %(NOT total in recovery class) from the heroin/fentanyl class, we use ODE K'=dy(8); going to run from
 %2013-2015 because those are the only years among 2013-2017 we have data for; 
 %new_heroin_admissions is for years 2014-2015, and y(1,8) is the number that are
 %admitted in the first year 2013
 new_heroin_admissions=zeros(1,2);
 for i=1:2
 new_heroin_admissions(i) = y(i+1,8)-y(i,8);
 end
 
 % yearly output from the model as a proportion that have been admitted into recovery class
 Estim3=[y(1,8), new_heroin_admissions(1), new_heroin_admissions(2)];
  
 % actual proportions of population each year being admitted into recovery from
 % heroin/fentanyl class
 Data3=[555./5517176 743./5559006 1083./5602117];
 %%Data simulated when put in 0.1 for all parameters
 %Data3=[0 0.01 0.011];
 % the difference between estimated value and data 
 diff3=Estim3-Data3;
 
 %print Estim3 values
 Estim3

 %%%%%
 % output from the model of the proportion of opioid addicts in 2015; we take
 % initial number of opioid addicts in 2015, y(2,3), and add the number of individuals that enter
% the A class at any point during the year 2015, which comes from
% integrating ODE L'=dy(9) but just focusing in on the one year, 2015, we care about 
 %(so have to subtract: integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2). 
 % MAY NEED TO FIX VALUE IF DECIDE THIS NUMBER COUNTS INDIVIDUALS IN RECOVERY, TOO 
 Estim4=[y(2,3)+y(3,9)-y(2,9)];
 % made up for now
 Data4=[48000./5602117];
 %Data simulated when put in 0.1 for all parameters
 %Data4=[0.12];
 % the difference between estimated value and data
 diff4=Estim4-Data4;
 
 %print Estim4 value 
 Estim4

 %%%%%
 % output from the model of the proportion of heroin/fentanyl addicts in
 % 2015; we take initial number of heroin/fentanyl addicts in 2015, y(2,4),
 % and add the number of individuals that enter
 % the H class at any point during the year 2015, which comes from
 % integrating ODE M'=dy(10) but just focusing in on the one year, 2015, we care about 
 %(so have to subtract: integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2). 
 Estim5=[y(2,4)+y(3,10)-y(2,10)];
 % Made up for now
 Data5=[14000./5602117];
 %Data simulated when put in 0.1 for all parameters
 %Data5=[0.11];
 % the difference between estimated value and data
 diff5=Estim5-Data5;
 
 
 %print Estim5 value 
 Estim5
 
 %%%%%
 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 %normalized by norm of the data
 value = norm(diff1,2)./norm(Data1)+ norm(diff2,2)./norm(Data2)+...
 norm(diff3,2)./norm(Data3)+norm(diff4,2)./norm(Data4)+norm(diff5,2)./norm(Data5);



 %{
 %Data points from Estim1, Data1
 figure(1)
 hold all
 scatter(1:1:5, Estim1,'filled')
 scatter(1:1:5, Data1, 'filled')
 set(gca, 'xtick', [1 2 3 4 5])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017'})
 xlabel('Year')
 ylabel('Proportion of prescription users')
 legend('Proportion entering P simulated','Proportion entering P data')

 
 
 Data points from Estim2, Data2
 figure(2)
 hold all
 scatter(1:1:3, Estim2,'filled')
 scatter(1:1:3, Data2, 'filled')
 set(gca, 'xtick', [1 2 3])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015'})
 xlabel('Year')
 ylabel('Proportion of new admissions into R from A')
 legend('Proportion entering R from A simulated','Proportion entering R from A data' )
 
 
 %Data points from Estim3, Data3
 figure(3)
 hold all
 scatter(1:1:3, Estim3,'filled')
 scatter(1:1:3, Data3, 'filled')
 set(gca, 'xtick', [1 2 3])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015'})
 xlabel('Year')
 ylabel('Proportion of new admissions into R from H')
 legend('Proportion entering R from H simulated','Proportion entering R from H data' )
  
 
 
 %Data points from Estim4, Data4
 figure(4)
 hold all
 scatter(1:1:1, Estim4,'filled')
 scatter(1:1:1, Data4, 'filled')
 set(gca, 'xtick', [ 1 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2015'})
 xlabel('Year')
 ylabel('Proportion of opioid addicts')
 legend('Proportion of pioid addicts simulated','Proportion of opioid addicts data' )
 
  
 %Data points from Estim5, Data5
 figure(5)
 hold all
 scatter(1:1:1, Estim5,'filled')
 scatter(1:1:1, Data5, 'filled')
 set(gca, 'xtick', [ 1 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2015'})
 xlabel('Year')
 ylabel('Proportion of heroin/fentanyl users')
 legend('Proportion of heroin/fentanyl users simulated','Proportion of heroin/fentanyl users data' )
 

 %Plots of ODES used in Estim1-Estim5
 figure(6)
          
           subplot(2,2,1);plot(t,y(:,6),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Proportion entering P')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           
           subplot(2,2,2);plot(t,y(:,7),'r-','LineWidth',1)
           hold on 
           subplot(2,2,2);plot(t,y(:,8) ,' g-','LineWidth',1)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Proportion entering R')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           legend('Proportion entering R from A','Proportion entering R from H' )
           xlim([0 , T])
           
           subplot(2,2,3);plot(t,y(:,9) ,' y-','LineWidth',1)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Proportion entering A')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
          
           subplot(2,2,4);plot(t,y(:,10) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Proportion entering H')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
           %}
 
 

 
           
         
 %ODE solutions plotted separately 
 figure(7)

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
           %ylim([0 , 1e+2])
           
          
           subplot(2,2,4);plot(t,y(:,5) ,' m-','LineWidth',1)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Individuals in recovery')
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
                
                 
 %ODE Solutions all plotted together
 figure(8)
           plot(t,y(:,2),'b-','LineWidth',1);
           hold all
           plot(t,y(:,3),'r-','LineWidth',1);
           hold all
           plot(t,y(:,4),' g-','LineWidth',1); 
           hold all
           plot(t,y(:,5),' m-','LineWidth',1); 
           xlabel('time')
           ylabel('Size of Populations');
           set(gca, 'xtick', [ 0 1 2 3 4 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
           legend('P','A','H','R')
           xlim([0 , T])
           ylim([0 , 0.2])
           legend('P','A','H','R')
           
           

 
