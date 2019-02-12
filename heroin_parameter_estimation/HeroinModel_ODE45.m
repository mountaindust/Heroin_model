% File name: HeroinModel_ODE45.m (used to be in Heroin_model folder)
function value = HeroinModel_ODE45(z)


% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-2017 where t=0 represents 2013 and t=4 represents 2017:
N = 4; 
T = N;

% Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
tspan=linspace(0,T,N+1);
global value
 
% Parameters
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

omega=0.00001;

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
  
 
 %COMPARING MODEL ESTIMATES TO DATA 
 
 %%%%%
 % In order to count the total number of individuals in P at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
  
 % To get the output from the model of the proportion of non-addicted
 % prescription opioid users in 2013 (first entry of Estim1), 
 % we take the total number of non-addiction prescription opioid users 
 % at the beginning of 2013 (P_0=IC, which is estimated because we do not
 % know the number at a particular instant at the beginning of the year)
 % and add on the number of individuals that enter the P class at any point during the year 2013, 
 % which comes fromintegrating ODE X'=dy(6) from t=0 to t=1; this gives
 % first value in Estim1. 
 
 % To get the output from the model of the proportion of non-addicted prescription opioid users in 2014, 2015
 % 2016, and 2017 (remaining entries of Estim1),
 % we take the initial number of non-addicted prescription opioid users 
 % in 2014: y(1,2), 2015: y(2,2), 2016: y(3,2), and 2017: y(4,2) and add the number of individuals that enter
 % the P class at any point during each of these years, which comes from
 % integrating ODE X'=dy(6) but just focusing in on these four specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 % for 2016, we have to subtract because integrating gives total number of new cases from t=0 to t=4, so have to 
 % subtract off the number from t=0 to t=3. 
 % for 2017, we have to subtract because integrating gives total number of new cases from t=0 to t=5, so have to 
 % subtract off the number from t=0 to t=4. 
 
 
 % Yearly output from the model as a proportion of P individuals for 2013-2017
 Estim1=[z(13)+y(1,6), y(1,2)+y(2,6)-y(1,6), y(2,2)+y(3,6)-y(2,6), y(3,2)+y(4,6)-y(3,6), y(4,2)+y(5,6)-y(4,6)];
 
 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 Data1=[1660630./5517176 1641908./5559006 1637623./5602117 1585227./5651993 1472737./5708586];
 
 % Data simulated when testing codes
 % Data1=[0.0057  0.2387  0.1860   0.1843  0.1842];
 
 % The difference between estimated value and data: 
 Diff1= Estim1-Data1;
 


 %%%%%
 % In order to count the total number of individuals in A at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of opioid addicts in 2014 and 2015 (Estim2),
 % we take the initial number of opioid addicts in 2014, y(1,3), and 2015, y(2,3),
 % and add the number of individuals that enter
 % the A class at any point during the year 2014 or 2015, which comes from
 % integrating ODE L'=dy(9) but just focusing in on these two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 
 % Yearly output from the model as a proportion of A individuals for 2014-2015 
 Estim2=[y(1,3)+y(2,7)-y(1,7), y(2,3)+y(3,7)-y(2,7)];
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 Data2=[42000./5651993 48000./5602117];
 
 % Data simulated when testing codes
 % Data2=[0.0113 0.0112];
 
 % The difference between estimated value and data
 Diff2=Estim2-Data2;
 

 %%%%%
 % In order to count the total number of individuals in H at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of heroin/fentanyl addicts in
 % 2014 and 2015 (Estim3), we take initial number of heroin/fentanyl addicts in 2014, y(2,3), 
 % and 2015, y(2,4),and add the number of individuals that enter the H class at any point
 % during the year 2014 or 2015, which comes from
 % integrating ODE M'=dy(10) but just focusing in on the two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1. 
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 
 % Yearly output from the model as a proportion of P individuals for 2015-2015 
 Estim3=[y(1,4)+y(2,8)-y(1,8), y(2,4)+y(3,8)-y(2,8), y(3,4)+y(4,8)-y(3,8)];
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 Data3=[14000./5559006 14000./5602117 19000./5651993];
 
 % Data simulated when testing codes
 % Data3=[0.0021 0.0022 0.0023];
 
 % The difference between estimated value and data
 Diff3=Estim3-Data3;
 
 
 %%%%%
 % STUDY MORE ABOUT 
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; 
 % gives least squares percentage error so each piece weighted evenly)
 value = norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);

 % Want value=f(x) to be small value when run MultiStart  