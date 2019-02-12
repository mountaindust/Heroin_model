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
 Estim1(1)=0.05+y(1,6);  
 % For 2014-final year:
    for i=2:25
       Estim1(i)= y(i-1,2)+y(i,6)-y(i-1,6);
    end
 
 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Data simulated when testing codes
 Data1=[0.0500000000000000,0.0437858768547988,0.0435673566034397,0.0435442092824961,0.0435376327082335,0.0435349036497300,0.0435135588646178,0.0434837571559799,0.0434650157573139,0.0434886581256490,0.0434504184350985,0.0434234275728627,0.0434029084840386,0.0434299680294224,0.0433891738680858,0.0433644490386973,0.0433374701116465,0.0433584497124210,0.0433293034062243,0.0433080608746508,0.0432771233946776,0.0432706969265679,0.0432886433415970,0.0432634791181139,0.0432471763079038];
 
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
 Estim2(1)=0.0062+y(1,7);  
 % For 2014-final year:
    for i=2:25
       Estim2(i)= y(i-1,3)+y(i,7)-y(i-1,7);
    end
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data2=[42000/5651993 48000/5602117];
 
 % Data simulated when testing codes
 Data2=[0.00620000000000000,0.00405433554272812,0.00276124782316749,0.00198866243896308,0.00152722132461754,0.00125151036401257,0.00108685937938708,0.000988448901509814,0.000929575577325199,0.000894196448057172,0.000873142093678575,0.000860483069916289,0.000852842470289440,0.000848052333429853,0.000845291524608411,0.000843552288979749,0.000842464418826318,0.000841599235538164,0.000841134707640878,0.000840781862473436,0.000840545192994104,0.000840275736233286,0.000839954178216526,0.000839797467581590,0.000839622311748481];
 
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
 Estim3(1)=0.00062+y(1,8);  
 % For 2014-final year:
    for i=2:25
       Estim3(i)= y(i-1,4)+y(i,8)-y(i-1,8);
    end
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes
 Data3=[0.000620000000000000,0.000569135839301842,0.000522515834892917,0.000479797752581738,0.000440659630647755,0.000404798025839577,0.000371929386649117,0.000341792160947820,0.000314148183757158,0.000288781898393431,0.000265498523439792,0.000244121643950184,0.000224491134942967,0.000206461222625318,0.000189898989374732,0.000174682916587517,0.000160701884826775,0.000147854143541916,0.000136046542303742,0.000125193672740303,0.000115217276291321,0.000106045591505194,9.76128280871813e-05,8.98586722741504e-05,8.27277714811727e-05];
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