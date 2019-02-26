%File name: heroin_multistart_final.m 

clf;

% We wish to estimate the parameter vector 
% x =[alpha,beta_A,beta_P,theta_1,epsilon,gamma,theta_2,sigma,zeta,theta_3,nu]
%Give ranges on each of the parameters 
LowerBounds=[0.01  0.00001 0.0001 0.00001 0.8 0.001 0.00001 0.01 0.01 0.0001 0.01];
UpperBounds=[0.7    0.1     0.009   0.1    4   0.1    0.3     2   1     0.6   1  ];
%Initial starting points for parameters, starting in the middle of each of the ranges
xstart=0.5*(LowerBounds + UpperBounds); 

% Create MultiStart problem using optimization function fmincon;
% x0 is xstart, objective is what we are trying to minimize which comes from 
% value = HeroinModel_ODE45(z) = fval(x) as output
problem=createOptimProblem('fmincon','objective',@HeroinModel_ODE45,...
         'x0', xstart,...
         'lb',LowerBounds,...
         'ub',UpperBounds);

%problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

% Define a multistart problem; results are reported after each local solver run, in addition to the final summary
ms=MultiStart('Display', 'iter'); 

% Number of times I want to run optimization scheme
numstartpoints=10;

% Manymins is a vector of solutions containing the distinct local minima found during the run;
%  runs MultiStart on numstartpoints to find a solution or multiple local solutions to problem
%[x,fval, exitflag, output, solutions]=run(ms,problem,numstartpoints);
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

alpha=x(1);
beta_A=x(2); 
beta_P=x(3);
theta_1=x(4);
epsilon=x(5);
mu=0.00868;  
mu_A=0.00775;   
mu_H=0.0271;
gamma=x(6);   
theta_2=x(7); 
sigma=x(8);
zeta=x(9);
theta_3=x(10);
nu=x(11);
omega=0.0000000001;

pars=[alpha,beta_A,beta_P,theta_1,epsilon,0.00868,0.00775,0.0271,gamma,theta_2,sigma,zeta,theta_3,nu,0.0000000001];

x
fval

% Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
N = 25; 
tspan=linspace(0,N,N+1);


% Initial conditions: assume for now 
% Later: although we know total number of prescription users in 2013, we do not
% know the initial number right at the start of 2013, so must be estimated;
% same with opioid addicts, heroin users, and stably recovered individuals.
S0=1-0.13-0.01-0.001-0.0003;
P0=0.13;
A0=0.01;
H0= 0.001;
R0=0.0003;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];


[t,y]=ode45(@HeroinModel,tspan,initials,[],pars);

 %Gives solution vector for each compartment 
  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);


 State1=y(:,1);
 State_data_1=[0.858700000000000;0.859494463775423;0.858748922841086;0.857888060790142;0.857012726633369;0.856143335717068;0.855429573386591;0.854364904792347;0.853660407345004;0.852866926507885;0.852090164613941;0.851279869123105;0.850622823206032;0.849718517221746;0.849010929616645;0.848316326549004;0.847627808045856;0.846896704487492;0.846243960086674;0.845527127998511;0.844837821720664;0.844187161871940;0.843571198906193;0.842963047939877;0.842249711450292;0.841659152142498];
 
 %Simulated data points from Susceptibles and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1))
 plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 legend('Proportion of susceptibles simulated','Proportion of susceptibles (simulated) data' )

 State2=y(:,2);
 State_data_2=[0.130000000000000;0.128152585424329;0.127865587080704;0.127709875635552;0.127584653930935;0.127469497873961;0.127214002052962;0.127326587283523;0.127093057124644;0.126964173553611;0.126833570563762;0.126751474355174;0.126530109061266;0.126571392521579;0.126429308527153;0.126288119940791;0.126154536342878;0.126077210481862;0.125934492622786;0.125869201770668;0.125789144468306;0.125682952403122;0.125554391508952;0.125430259570332;0.125423820551677;0.125305992978402];
 
 
 %Simulated data points from Prescription Users and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2))
 plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 legend('Proportion of prescription users simulated','Proportion of prescription users (simulated) data' )


 State3=y(:,3);
 State_data_3=[0.0100000000000000;0.0104810615372110;0.0111854899699225;0.0119844214734340;0.0128181287745922;0.0136590321538875;0.0144958539257199;0.0153214335394139;0.0161362840340819;0.0169377687796236;0.0177261097843963;0.0185011024044664;0.0192636824884010;0.0200122968625343;0.0207490676501564;0.0214734032456553;0.0221854777372510;0.0228852813088750;0.0235735362798669;0.0242498181648154;0.0249147145455269;0.0255684591523685;0.0262112170574558;0.0268430503395490;0.0274636446369883;0.0280741593962154];
 
 
 %Simulated data points from Opioid Addicts and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3))
 plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend('Proportion of opioid addicts simulated','Proportion of opioid addicts (simulated) data' )

 
 State4=y(:,4);
 State_data_4=[0.00100000000000000;0.000955552913699605;0.000936827079652513;0.000929247591160181;0.000926810471213153;0.000926934121294542;0.000928524446324821;0.000931088681845971;0.000934371037338425;0.000938242761454411;0.000942625993198063;0.000947465790457975;0.000952723853644378;0.000958364078721766;0.000964365369072255;0.000970703789798523;0.000977360591214322;0.000984319193538382;0.000991566923336990;0.000999090162812635;0.00100687925692138;0.00101492500608098;0.00102321920356836;0.00103175441121314;0.00104052328020724;0.00104952239283311];
 
 %Simulated data points from Heroin/Fentanyl Addicts and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4))
 plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend('Proportion of heroin/fentanyl addicts simulated','Proportion of heroin/fentanyl addicts (simulated) data' )
 

 State5=y(:,5);
 State_data_5=[0.000300000000000000;0.000916336349337721;0.00126317302863535;0.00148839450971127;0.00165768018989119;0.00180120013378857;0.00193204618840229;0.00205598570287056;0.00217588045893161;0.00229288839742625;0.00240752904470281;0.00252008832679688;0.00263066139065648;0.00273942931541916;0.00284632883697332;0.00295144647475144;0.00305481728280103;0.00315648452823205;0.00325644408733570;0.00335476190319288;0.00345144000858277;0.00354650156648922;0.00363997332383043;0.00373188773902871;0.00382230008083507;0.00391117309005228];

  %Simulated data points from Stably Recovered Individuals and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5))
 plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 legend('Proportion of stably recovered addicts simulated','Proportion of stably recovered addicts (simulated) data' )

 %{
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 Data1=[0.387803312960829;0.385859884749359;0.385291909033797;0.384844182844686;0.384426363863043;0.384001142924548;0.383502489332366;0.383286474912885;0.382780015767370;0.382368169620462;0.381962853339432;0.381583803634593;0.381120645490644;0.380863382568114;0.380449265557906;0.380039206710655;0.379645083482821;0.379293608791355;0.378896962400306;0.378569085690886;0.378225964168557;0.377859610734567;0.377476726514315;0.377114708562721;0.376842862140627];
 
 %Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t(1:end-1),Estim1)
 plot(t(1:end-1), Data1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('Data points interested in', 'ODE solution')

 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 Data2=[0.0117455980073673;0.0126003611353261;0.0135702699068408;0.0145862134466473;0.0156163578770581;0.0166472360917028;0.0176707728614604;0.0186868024270982;0.0196925574245404;0.0206879590347452;0.0216726879014664;0.0226474706126398;0.0236108065439615;0.0245644710781424;0.0255079130690137;0.0264411627138692;0.0273641185151985;0.0282773526420237;0.0291804040459065;0.0300737126980973;0.0309574258792886;0.0318316178524296;0.0326962699927304;0.0335510122628854;0.0343968111832996];
 %Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t(1:end-1),Estim2)
 plot(t(1:end-1), Data2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('Data points interested in', 'ODE solution')

 Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);  
 Data3=[0.00103907582372071;0.00101755059193969;0.00100856050223205;0.00100532420926167;0.00100495607550998;0.00100622617515543;0.00100858046437510;0.00101173893007554;0.00101555379888286;0.00101993773451460;0.00102482901368315;0.00103018380897159;0.00103596238440371;0.00104213951408922;0.00104868879606276;0.00105558903032135;0.00106282166296150;0.00107037211192985;0.00107822549298304;0.00108637072936078;0.00109479753199721;0.00110349676404172;0.00111246022459328;0.00112168003192911;0.00113115181543487];
 
 %Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t(1:end-1),Estim3)
 plot(t(1:end-1), Data3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('Data points interested in', 'ODE solution')
 
 %}
 
function value = HeroinModel_ODE45(z)

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
omega=0.0000000001;


pars=[alpha,beta_A,beta_P,theta_1,epsilon,0.00868,0.00775,0.0271,gamma,theta_2,sigma,zeta,theta_3,nu,0.0000000001];

% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-? where t=0 represents 2013
% and t=N represents ?; with spacing (T-0)/((N+1)-1)=1 between the points:
N = 25; 
tspan=linspace(0,N,N+1);

% Initial conditions: assumed for now; although we know total number of prescription users in 2013, we do not
% know the initial number right at the start of 2013, so must be estimated;
% same with opioid addicts, heroin users, and stably recovered individuals.
S0=1-0.13-0.01-0.001-0.0003;
P0=0.13;
A0=0.01;
H0= 0.001;
R0=0.0003;
X0=0;
L0=0;
M0=0;
initials = [S0,P0,A0,H0,R0,X0,L0,M0];


[t,y]=ode45(@HeroinModel,tspan,initials,[],pars);
  
 
 
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

 value=0;
 % Yearly output from the model as a proportion of P individuals for
 % 2013-final year, Estim1 is a row vector
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  

 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Data simulated when testing codes (Data1)
 Data1=[0.387803312960829;0.385859884749359;0.385291909033797;0.384844182844686;0.384426363863043;0.384001142924548;0.383502489332366;0.383286474912885;0.382780015767370;0.382368169620462;0.381962853339432;0.381583803634593;0.381120645490644;0.380863382568114;0.380449265557906;0.380039206710655;0.379645083482821;0.379293608791355;0.378896962400306;0.378569085690886;0.378225964168557;0.377859610734567;0.377476726514315;0.377114708562721;0.376842862140627];
 
 Diff1=Estim1-Data1; 
 % The difference between estimated value and data: 
% State1=y(1:end,1);
 State1=y(:,1);
 State_data_1=[0.858700000000000;0.859494463775423;0.858748922841086;0.857888060790142;0.857012726633369;0.856143335717068;0.855429573386591;0.854364904792347;0.853660407345004;0.852866926507885;0.852090164613941;0.851279869123105;0.850622823206032;0.849718517221746;0.849010929616645;0.848316326549004;0.847627808045856;0.846896704487492;0.846243960086674;0.845527127998511;0.844837821720664;0.844187161871940;0.843571198906193;0.842963047939877;0.842249711450292;0.841659152142498];
 
 State_diff_1= State1-State_data_1;
 


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

 
 % Yearly output from the model as a proportion of A individuals for
 % 2013-final year, Estim2 is a row vector
 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
% Data2=[42000/5651993 48000/5602117];
 
 % Data simulated when testing codes (Data2)
 Data2=[0.0117455980073673;0.0126003611353261;0.0135702699068408;0.0145862134466473;0.0156163578770581;0.0166472360917028;0.0176707728614604;0.0186868024270982;0.0196925574245404;0.0206879590347452;0.0216726879014664;0.0226474706126398;0.0236108065439615;0.0245644710781424;0.0255079130690137;0.0264411627138692;0.0273641185151985;0.0282773526420237;0.0291804040459065;0.0300737126980973;0.0309574258792886;0.0318316178524296;0.0326962699927304;0.0335510122628854;0.0343968111832996];
 
 Diff2=Estim2-Data2; 
 % The difference between estimated value and data
 %State2=y(1:end,2);
 State2=y(:,2);
 State_data_2=[0.130000000000000;0.128152585424329;0.127865587080704;0.127709875635552;0.127584653930935;0.127469497873961;0.127214002052962;0.127326587283523;0.127093057124644;0.126964173553611;0.126833570563762;0.126751474355174;0.126530109061266;0.126571392521579;0.126429308527153;0.126288119940791;0.126154536342878;0.126077210481862;0.125934492622786;0.125869201770668;0.125789144468306;0.125682952403122;0.125554391508952;0.125430259570332;0.125423820551677;0.125305992978402];
 
 State_diff_2=State2-State_data_2;
 

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
 

 % Yearly output from the model as a proportion of H individuals for
 % 2013-final year, Data3 is a row vector 
 Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);  
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes (Data3)
 Data3=[0.00103907582372071;0.00101755059193969;0.00100856050223205;0.00100532420926167;0.00100495607550998;0.00100622617515543;0.00100858046437510;0.00101173893007554;0.00101555379888286;0.00101993773451460;0.00102482901368315;0.00103018380897159;0.00103596238440371;0.00104213951408922;0.00104868879606276;0.00105558903032135;0.00106282166296150;0.00107037211192985;0.00107822549298304;0.00108637072936078;0.00109479753199721;0.00110349676404172;0.00111246022459328;0.00112168003192911;0.00113115181543487];
 
 Diff3=Estim3-Data3;
 % The difference between estimated value and data
 %State3=y(1:end,3);
 State3=y(:,3);
 State_data_3=[0.0100000000000000;0.0104810615372110;0.0111854899699225;0.0119844214734340;0.0128181287745922;0.0136590321538875;0.0144958539257199;0.0153214335394139;0.0161362840340819;0.0169377687796236;0.0177261097843963;0.0185011024044664;0.0192636824884010;0.0200122968625343;0.0207490676501564;0.0214734032456553;0.0221854777372510;0.0228852813088750;0.0235735362798669;0.0242498181648154;0.0249147145455269;0.0255684591523685;0.0262112170574558;0.0268430503395490;0.0274636446369883;0.0280741593962154];
 
 State_diff_3=State3-State_data_3;
 
 %State4=y(1:end,4);
 State4=y(:,4);
 State_data_4=[0.00100000000000000;0.000955552913699605;0.000936827079652513;0.000929247591160181;0.000926810471213153;0.000926934121294542;0.000928524446324821;0.000931088681845971;0.000934371037338425;0.000938242761454411;0.000942625993198063;0.000947465790457975;0.000952723853644378;0.000958364078721766;0.000964365369072255;0.000970703789798523;0.000977360591214322;0.000984319193538382;0.000991566923336990;0.000999090162812635;0.00100687925692138;0.00101492500608098;0.00102321920356836;0.00103175441121314;0.00104052328020724;0.00104952239283311];
 
 State_diff_4=State4-State_data_4;
 
 
 %State5=y(1:end,5);
 State5=y(:,5);
 State_data_5=[0.000300000000000000;0.000916336349337721;0.00126317302863535;0.00148839450971127;0.00165768018989119;0.00180120013378857;0.00193204618840229;0.00205598570287056;0.00217588045893161;0.00229288839742625;0.00240752904470281;0.00252008832679688;0.00263066139065648;0.00273942931541916;0.00284632883697332;0.00295144647475144;0.00305481728280103;0.00315648452823205;0.00325644408733570;0.00335476190319288;0.00345144000858277;0.00354650156648922;0.00363997332383043;0.00373188773902871;0.00382230008083507;0.00391117309005228];
 
 State_diff_5=State5-State_data_5;
 
 %%%%%
 % STUDY MORE ABOUT 
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sqrt(sum from 1 to N of (diff#)^2)
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; 
 % gives least squares percentage error so each piece weighted evenly)
 
 %value=norm(State_diff_1,2)+norm(State_diff_2,2)+norm(State_diff_3,2)+norm(State_diff_4,2)+norm(State_diff_5,2);
 
 value = norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5);
 
 %value=sum(State_diff_1.^2);
 
 %value=norm(Diff1,2)+norm(Diff2,2)+norm(Diff3,2);
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5);
 
 % Want value=f(x) to be small value when run MultiStart  
 
end
 

function f = HeroinModel(t,y,pars)
f=zeros(5,1);
f(1)=-pars(1)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=pars(1)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
f(6) = pars(1)*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));
% Transpose column vector solution into row vector 
 % f=f';  
 
end








           
 

