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
 Estim1(1)=P0+y(1,6);  
 % For 2014-final year:
    for i=2:25
       Estim1(i)= y(i-1,2)+y(i,6)-y(i-1,6);
    end

 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Data simulated when testing codes (Estim1 data)
 Data1=[0.100000000000000,0.236920644285165,0.197496735756373,0.194699419098524,0.194452716509407,0.194388887657194,0.194329456360907,0.194290124565391,0.194246194271848,0.194200669271217,0.194152508688169,0.194099178081385,0.194049217468572,0.194033063704699,0.193968252418088,0.193925125309202,0.193874729431561,0.193837123814976,0.193813047326274,0.193752213403401,0.193711478104666,0.193663892319981,0.193629030434897,0.193608830035064,0.193548401618371];
 
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
 Estim2(1)=A0+y(1,7);  
 % For 2014-final year:
    for i=2:25
       Estim2(i)= y(i-1,3)+y(i,7)-y(i-1,7);
    end
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data2=[42000/5651993 48000/5602117];
 
 % Data simulated when testing codes (Estim2 data)
 Data2=[0.00570000000000000,0.00938467924914681,0.0114519564665786,0.0127778460727290,0.0136787440366524,0.0143323642818403,0.0148405968141525,0.0152632719362520,0.0156342161103079,0.0159734217737802,0.0162922444615463,0.0165972538437866,0.0168921743190029,0.0171795884089407,0.0174605638519617,0.0177360075537735,0.0180063391627154,0.0182718606710802,0.0185329955103435,0.0187896784622578,0.0190421337250643,0.0192903897693775,0.0195345115752703,0.0197748060458852,0.0200111092324205];
 
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
 Estim3(1)=H0+y(1,8);  
 % For 2014-final year:
    for i=2:25
       Estim3(i)= y(i-1,4)+y(i,8)-y(i-1,8);
    end
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes: NEED TO CHANGE TO BE ESTIM1 DATA, not
 % P data
 Data3=[0.00130000000000000,0.00193870237034378,0.00219291186092112,0.00228668767354006,0.00229636119393789,0.00226229757603101,0.00220603303143841,0.00213922769695117,0.00206805107285677,0.00199583766890079,0.00192432346358949,0.00185440202431088,0.00178651470839317,0.00172085606563564,0.00165748825531202,0.00159640340181048,0.00153755715171060,0.00148088734862659,0.00142632371160503,0.00137379225022955,0.00132321917484624,0.00127453163018483,0.00122765878302782,0.00118253237471051,0.00113908598281541];
 
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