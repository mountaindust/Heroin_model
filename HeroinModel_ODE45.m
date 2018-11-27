%File name: HeroinModel_ODE45.m
function value = HeroinModel_ODE45(z)


% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-2017 where t=0 represents 2013 and t=4 represents 2017:
N = 4; 
T = N;

%Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
tspan=linspace(0,T,N+1);
global value
 
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
S0=1-z(13)-z(14)-z(15)-z(16); 
P0=z(13);
A0=z(14);
H0=z(15);
R0=z(16); 
X0=0;
Z0=0;
K0=0;
initials = [S0,P0,A0,H0,R0,X0,Z0,K0];



[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  Z=y(:,7);
  K=y(:,8);
  
  

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
 Estim1=[z(13)+y(1,6), total_prescription_users(1),total_prescription_users(2),total_prescription_users(3), total_prescription_users(4)];
       
% actual proportions of population that were prescription opioid users for
% 2013-2017 (total number of prescription opioid users in each year in TN that are 12 and older divided by
% total population in TN 12 and older) 
 Data1=[1845144./5517176 1824342./5559006 1819581./5602117 1761363./5651993 1636374./5708586];
 
% the difference between estimated value and data 
 diff1= Estim1-Data1;
 
 
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

 % the difference between estimated value and data 
 diff2=Estim2-Data2;
 
 
 
 %%%%%
 %To calculate number of new admissions coming into the recovery class 
 %(NOT total in recovery class) from the heroin/fentanyl class, we use ODE K'=dy(8); going to run from
 %2013-2015 because those are the only years among 2013-2017 we have data for; 
 %new_heroin_admissions is for years 2014-2015, and y(1,8) is the number that are
 %admitted in the first year 2013
 new_heroin_admissions=zeros(1,2);
 for i=1:2
 new_heroin_admissions(i) = y(i+1,7)-y(i,7);
 end
 
 % yearly output from the model as a proportion that have been admitted into recovery class
 Estim3=[y(1,8), new_heroin_admissions(1), new_heroin_admissions(2)];
  
 % actual proportions of population each year being admitted into recovery from
 % heroin/fentanyl class
 Data3=[555./5517176 743./5559006 1083./5602117];

 % the difference between estimated value and data 
 diff3=Estim3-Data3;
 

 %%%%%
 % output from the model of the proportion of opioid addicts in 2015
 % MAY NEED TO FIX VALUE IF DECIDE THIS NUMBER COUNTS INDIVIDUALS IN RECOVERY, TOO 
 Estim4=[y(2,3)];
 % made up for now
 Data4=[48000./5602117];
 % the difference between estimated value and data
 diff4=Estim4-Data4;
 

 %%%%%
 % output from the model of the proportion of heroin/fentanyl addicts in 2015
 Estim5=[y(2,4)];
 % Made up for now
 Data5=[14000./5602117];
 % the difference between estimated value and data
 diff5=Estim5-Data5;
 
 
 %%%%%
 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 %normalized by norm of the data
 value = norm(diff1,2)./norm(Data1)+ norm(diff2,2)./norm(Data2)+...
 norm(diff3,2)./norm(Data3)+norm(diff4,2)./norm(Data4)+norm(diff5,2)./norm(Data5);