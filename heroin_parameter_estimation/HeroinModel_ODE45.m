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
 
mu=0.00868;  
 
mu_A=0.00775;   
 
mu_H=0.0271;
 
gamma=z(6);   
 
theta_2=z(7); 
 
sigma=z(8);
 
zeta=z(9);
 
theta_3=z(10);
 
nu=z(11);

%Although we know total number of prescription users in 2013, we do not
%know the initial number right at the start of 2013, so must be estimated
P0=z(12);

%Do not know number of opioid addicts at start of 2013
A0=z(13);

%Do not know number of heroin users at start of 2013
H0=z(14);

%Do not know number of individuals in recovery at start of 2013
R0=z(15);


%Initials
S0=1-z(12)-z(13)-z(14)-z(15); 
P0=z(12);
A0=z(13);
H0=z(14);
R0=z(15); 
X0=0;
%Z0=0;
%K0=0;
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
 % Z=y(:,7);
 % K=y(:,8);
  L=y(:,7);
  M=y(:,8);
  
  

 %COMPARING MODEL ESTIMATES TO DATA 
 
 %%%%%
 % In order to count the total number of individuals in P at some point throughout a certain year
 % (the number who are in the class AT ALL during the year, it's okay if they leave), 
 % we want to add the total number of individuals using prescriptions at the beginning of the
 % year to those who come into P at some point during the year: 
 % prescript(i) is number of prescription users at beginning of the year (time step) +
 % the number of new cases that came in from X'=dy(6) ODE from that year until the beginning 
 % of the next year. Thus, we are adding in those who come into P at some point during 
 % the year by integrating X'=dy(6) ODE but just focusing in on the one year we care about 
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
% Data1=[1845144./5517176 1824342./5559006 1819581./5602117 1761363./5651993 1636374./5708586];
 %Data simulated when testing codes
 Data1=[0.1  0.2387  0.1861   0.1843  0.1842];
% the difference between estimated value and data 
 diff1= Estim1-Data1;
 
 
%%%%%
%CANNOT USE ANYMORE SINCE CHANGED DEFINITION OF "RECOVERY CLASS"
 %To calculate number of new admissions coming into the recovery class 
 %(NOT total in recovery class) from the opioid addict class, we use ODE Z'=dy(7); going to run from
 %2013-2015 because those are the only years among 2013-2017 we have data for; 
 %new_opioid_admissions is for years 2014-2015, and y(1,7) is the number that are
 %admitted in the first year 2013
 
 %new_opioid_admissions=zeros(1,2);
 %for i=1:2
 %new_opioid_admissions(i) = y(i+1,7)-y(i,7);
 %end
 
 % yearly output from the model as a proportion in the recovery class
 %Estim2=[y(1,7), new_opioid_admissions(1), new_opioid_admissions(2)];
  
 % actual proportions of population each year being admitted into recovery from opioid addict class
% Data2=[4485./5517176 4530./5559006 4326./5602117];
  %Data simulated when testing codes
 %Data2=[0 0.0003677 0.0004421];

 % the difference between estimated value and data 
% diff2=Estim2-Data2;
 
 
 
 %%%%%
 %CANNOT USE ANYMORE SINCE CHANGED DEFINITION OF "RECOVERY CLASS"
 %To calculate number of new admissions coming into the recovery class 
 %(NOT total in recovery class) from the heroin/fentanyl class, we use ODE K'=dy(8); going to run from
 %2013-2015 because those are the only years among 2013-2017 we have data for; 
 %new_heroin_admissions is for years 2014-2015, and y(1,8) is the number that are
 %admitted in the first year 2013
 
 %new_heroin_admissions=zeros(1,2);
 %for i=1:2
 %new_heroin_admissions(i) = y(i+1,8)-y(i,8);
 %end
 
 % yearly output from the model as a proportion that have been admitted into recovery class
 %Estim3=[y(1,8), new_heroin_admissions(1), new_heroin_admissions(2)];
  
 % actual proportions of population each year being admitted into recovery from
 % heroin/fentanyl class
 %Data3=[555./5517176 743./5559006 1083./5602117];
 %Data simulated when testing codes
% Data3=[0  0.0001825  0.0002859];
 % the difference between estimated value and data 
 %diff3=Estim3-Data3;
 

 %%%%%
 % output from the model of the proportion of opioid addicts in 2015; we take
 % initial number of opioid addicts in 2015, y(2,3), and add the number of individuals that enter
% the A class at any point during the year 2015, which comes from
% integrating ODE L'=dy(9) but just focusing in on the one year, 2015, we care about 
 %(so have to subtract: integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2). 
 Estim4=[y(2,3)+y(3,7)-y(2,7)];
 %actual proportion of opioid addicted individuals in the population in 2015
 %Data4=[48000./5602117];
 %Data simulated when testing codes
 Data4=[0.0097];
 % the difference between estimated value and data
 diff4=Estim4-Data4;
 

 %%%%%
 % output from the model of the proportion of heroin/fentanyl addicts in
 % 2015; we take initial number of heroin/fentanyl addicts in 2015, y(2,4),
 % and add the number of individuals that enter
 % the H class at any point during the year 2015, which comes from
 % integrating ODE M'=dy(10) but just focusing in on the one year, 2015, we care about 
 %(so have to subtract: integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2). 
 Estim5=[y(2,4)+y(3,8)-y(2,8)];
 % Made up for now
 %Data5=[14000./5602117];
 %Data simulated when testing codes
 Data5=[0.0066];
 % the difference between estimated value and data
 diff5=Estim5-Data5;
 
 
 %%%%%
 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 %normalized by norm of the data
 value = norm(diff1,2)./norm(Data1)+norm(diff4,2)./norm(Data4)+norm(diff5,2)./norm(Data5);

%Want value=f(x) to be small value when run MultiStart  