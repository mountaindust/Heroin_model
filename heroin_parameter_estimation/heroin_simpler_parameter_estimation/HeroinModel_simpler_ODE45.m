% File name: HeroinModel_ODE45.m (used to be in Heroin_model folder)
function value = HeroinModel_simpler_ODE45(z)


% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-2017 where t=0 represents 2013 and t=4 represents 2017:
N = 25; 
T = N;
% Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
tspan=linspace(0,T,N+1);
% global value
 
% Initial conditions

S0=1-0.05-0.0062-0.00062-0.00062; 
P0=0.05;
A0=0.0062;
H0=0.00062;
R0=0.00062;
X0=0;
L0=0;
M0=0;
initials = [S0,P0,A0,H0,R0,X0,L0,M0];


[t,y]=ode45(@(t,y) HeroinModel_simpler(t,y,z),tspan,initials);


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

 
 % Yearly output from the model as a proportion of P individuals for 2013-2017 (individually typing in equations, same as below with for loop):
 %Estim1=[z(12)+y(1,6), y(1,2)+y(2,6)-y(1,6), y(2,2)+y(3,6)-y(2,6), y(3,2)+y(4,6)-y(3,6), y(4,2)+y(5,6)-y(4,6)];
 
 
 % Yearly output from the model as a proportion of P individuals for
 % 2013-final year, Estim1 is a row vector
 Estim1=zeros(1,25);
 % For 2013:
 Estim1(1)=P0+y(2,6);  
 % For 2014-final year:
    for i=2:25
       Estim1(i)= y(i,2)+y(i+1,6)-y(i,6);
    end
 
 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Estimations to use to test code with simulated data with proportions in P:
 %Estim1=zeros(1,26);  
 % For 2013-final year:
    %for i=1:26
    %   Estim1(i)=y(i,2);
   % end
 
 % Data simulated when testing codes (Estim1 data)
 Data1=[0.188316779660538,0.217794221758152,0.223222856140525,0.224120312140478,0.224157326248006,0.224039108266267,0.223895051697782,0.223731799823345,0.223608204028267,0.223424310558097,0.223286018518309,0.223137118315429,0.222998912872893,0.222849914695050,0.222726481869693,0.222575998870438,0.222452913163314,0.222316736063000,0.222172768835784,0.222058084569280,0.221914172633211,0.221795716713520,0.221677976757133,0.221545160666303,0.221432788017319];
 
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

 
 % Yearly output from the model as a proportion of A individuals for 2014-2015 (individually typing in equations, same as below with for loop except that adds on beginning year): 
 %Estim2=[y(1,3)+y(2,7)-y(1,7), y(2,3)+y(3,7)-y(2,7)];
 
  
 % Yearly output from the model as a proportion of A individuals for
 % 2013-final year, Estim2 is a row vector
 Estim2=zeros(1,25);  
 % For 2013:
 Estim2(1)=A0+y(2,7);  
 % For 2014-final year:
    for i=2:25
       Estim2(i)= y(i,3)+y(i+1,7)-y(i,7);
    end
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data2=[42000/5651993 48000/5602117];
 
 
  
 % Estimations to use to test code with simulated data for proportions in A
 %Estim2=zeros(1,26);  
 % For 2013-final year:
    %for i=1:26
      % Estim2(i)=y(i,3);
   % end
 
 % Data simulated when testing codes (Estim2 data)
 Data2=[0.00718962158934690,0.00792087052033291,0.00866262015698704,0.00939190998540586,0.0101064909321721,0.0108072700283992,0.0114953978844816,0.0121712404780388,0.0128359918511854,0.0134893773557100,0.0141320564820926,0.0147641389964138,0.0153860047360764,0.0159975932818400,0.0165994545619014,0.0171914090539077,0.0177738717671610,0.0183470120726527,0.0189107133946403,0.0194655299149244,0.0200112679503111,0.0205482073688649,0.0210766253495816,0.0215964241679755,0.0221079310501131];
 
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
 
 % Yearly output from the model as a proportion of P individuals for 2014-2016 (individually typing in equations, same as below with for loop except that adds on beginning year): 
 %Estim3=[y(1,4)+y(2,8)-y(1,8), y(2,4)+y(3,8)-y(2,8), y(3,4)+y(4,8)-y(3,8)];
 
 
 % Yearly output from the model as a proportion of H individuals for
 % 2013-final year, Estim3 is a row vector 
 Estim3=zeros(1,25);  
 % For 2013:
 Estim3(1)=H0+y(2,8);  
 % For 2014-final year:
    for i=2:25
      Estim3(i)= y(i,4)+y(i+1,8)-y(i,8);
    end
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
  
 % Estimations to use to test code with simulated data for proportions in H
 %Estim3=zeros(1,26);  
 % For 2013-final year:
    %for i=1:26
      % Estim3(i)=y(i,4);
    %end
 
 % Data simulated when testing codes (Estim3 data)
 Data3=[0.000647349829272590,0.000620753228979919,0.000594848411244155,0.000569939553334447,0.000546148283453608,0.000523489554771551,0.000501927591302471,0.000481407407827558,0.000461867610495608,0.000443247409121684,0.000425490386163178,0.000408543482517175,0.000392358472359807,0.000376890921183491,0.000362100622905131,0.000347950149100495,0.000334405443854815,0.000321434932271192,0.000309009139487501,0.000297101133759905,0.000285685354189357,0.000274738219891521,0.000264237621250097,0.000254162656212776,0.000244493859803886];
 
 % The difference between estimated value and data
 Diff3=Estim3-Data3;
 
 
 %%%%%
 % STUDY MORE ABOUT 
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sqrt(sum from 1 to N of (diff#)^2)
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; 
 % gives least squares percentage error so each piece weighted evenly)
 value = norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);

 % Want value=f(x) to be small value when run MultiStart  