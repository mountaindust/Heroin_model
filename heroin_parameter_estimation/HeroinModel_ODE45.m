% File name: HeroinModel_ODE45.m (used to be in Heroin_model folder)
function value = HeroinModel_ODE45(z)


% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-2017 where t=0 represents 2013 and t=4 represents 2017:
N = 25; 
T = N;
% Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
tspan=linspace(0,T,N+1);
% global value
 
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

 
 % Yearly output from the model as a proportion of P individuals for 2013-2017 (individually typing in equations, same as below with for loop):
 %Estim1=[z(12)+y(1,6), y(1,2)+y(2,6)-y(1,6), y(2,2)+y(3,6)-y(2,6), y(3,2)+y(4,6)-y(3,6), y(4,2)+y(5,6)-y(4,6)];
 
 
 % Yearly output from the model as a proportion of P individuals for
 % 2013-final year, Estim1 is a row vector
 Estim1=zeros(1,25);
 % For 2013:
 Estim1(1)=z(12)+y(1,6);  
 % For 2014-final year:
    for i=2:25
       Estim1(i)= y(i-1,2)+y(i,6)-y(i-1,6);
    end
 
 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Data simulated when testing codes
 Data1=[0.0450594908386941 0.0432485942841840 0.0431809613437141 0.0431746328338971 0.0431509231830043 0.0431572067366684 0.0431577566860147 0.0431531937401041 0.0431347533745317 0.0431397645886815 0.0431373119796772 0.0431340779134710 0.0431302208131826 0.0431252747679749 0.0431148760419780 0.0431125795370239 0.0431216311030822 0.0431165112228240 0.0431098534831810 0.0430851299333931 0.0430725105777350 0.0430792376956774 0.0431113861643432 0.0430975797808445 0.0430623842165597];
 
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
 Estim2(1)=z(13)+y(1,7);  
 % For 2014-final year:
    for i=2:25
       Estim2(i)= y(i-1,3)+y(i,7)-y(i-1,7);
    end
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data2=[42000/5651993 48000/5602117];
 
 % Data simulated when testing codes
 Data2=[0.00877395531086108 0.0105557126141506 0.0116682100213010 0.0123948986983610 0.0128976146814453 0.0132694400674311 0.0135647071650127 0.0138144599264532 0.0140367046970982 0.0142417046066421 0.0144355690075931 0.0146217506851487 0.0148022742661666 0.0149783276408538 0.0151506296277036 0.0153195670387838 0.0154853970872745 0.0156483686527950 0.0158085651085436 0.0159661240743621 0.0161210279063084 0.0162733152031094 0.0164230213041524 0.0165704043568830 0.0167154384087158];
 
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
 Estim3(1)=z(14)+y(1,8);  
 % For 2014-final year:
    for i=2:25
       Estim3(i)= y(i-1,4)+y(i,8)-y(i-1,8);
    end
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes
 Data3=[0.00180537360151720 0.00203550456542268 0.00212183752598483 0.00213179027966765 0.00210158358946534 0.00205077268415127 0.00199000995916635 0.00192501892113224 0.00185888739670270 0.00179325702916532 0.00172898611797921 0.00166650161876701 0.00160599882845433 0.00154754727901121 0.00149114939911912 0.00143677190698099 0.00138436349911963 0.00133386429675168 0.00128521069777462 0.00123833840177292 0.00119318358922620 0.00114968394992442 0.00110777897037671 0.00106741033083450 0.00102852141783479];
 
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