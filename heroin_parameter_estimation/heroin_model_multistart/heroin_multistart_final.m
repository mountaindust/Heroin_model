%File name: heroin_multistart_final.m 

function heroin_multistart_final
clf

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
numstartpoints=100;

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
 State_data_1=[0.858700000000000;0.859426382221539;0.858612750962439;0.857680026489993;0.856728627418263;0.855778719000533;0.854981985805504;0.853823981529086;0.853027641683518;0.852136485264308;0.851258558244580;0.850343756535407;0.849577072330486;0.848561235741790;0.847738135336314;0.846924588052323;0.846114111479355;0.845257995125053;0.844477966333902;0.843629583190915;0.842807390205450;0.842020984909400;0.841266458741277;0.840517686175110;0.839661576789117;0.838925559955817];
 
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
 State_data_2=[0.130000000000000;0.128146752382792;0.127850272546011;0.127684209002999;0.127547822881584;0.127420861575431;0.127150984714692;0.127254640305519;0.127005918887905;0.126863057848528;0.126717499541811;0.126619468588702;0.126383153635609;0.126407329040188;0.126248261503800;0.126089783843546;0.125938322363150;0.125842698871172;0.125680620681888;0.125597021425835;0.125496899221430;0.125370537843949;0.125221780321309;0.125076774697224;0.125048976016212;0.124909978261869];
 
 
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
 State_data_3=[0.0100000000000000;0.0105521986905989;0.0113279484674021;0.0122010261109671;0.0131127408816685;0.0140358756701086;0.0149592278436868;0.0158755538671669;0.0167852882262808;0.0176856622676330;0.0185767718376313;0.0194582876246347;0.0203309911869715;0.0211932556342387;0.0220470263202622;0.0228916103977938;0.0237270657018265;0.0245532718684156;0.0253708274478466;0.0261792051297063;0.0269788861682556;0.0277699903595808;0.0285525798990412;0.0293266223811932;0.0300917205162188;0.0308489038337722];
 
 
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
 State_data_4=[0.00100000000000000;0.000955474051020676;0.000936517341675051;0.000928614467030446;0.000925813409980114;0.000925563989978451;0.000926786882742384;0.000928992670475523;0.000931928796634418;0.000935465244291553;0.000939522999751411;0.000944045900759262;0.000948994285772316;0.000954331395841673;0.000960034832376391;0.000966080092296467;0.000972447827651086;0.000979121005419110;0.000986086427265061;0.000993330263648804;0.00100084252305340;0.00100861375781570;0.00101663559133985;0.00102490048518095;0.00103340112048735;0.00104213374400898];
 
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
 State_data_5=[0.000300000000000000;0.000919192654049952;0.00127251068247289;0.00150612392901015;0.00168499540850452;0.00183897976394859;0.00198101475337481;0.00211683162775210;0.00224922240566137;0.00237932937523934;0.00250764737622560;0.00263444135049720;0.00275978856116059;0.00288384818794147;0.00300654200724719;0.00312793761404071;0.00324805262801645;0.00336691312994073;0.00348449910909832;0.00360085998989521;0.00371598188181106;0.00382987312925404;0.00394254544703346;0.00405401626129209;0.00416432555796432;0.00427342420453181];
 
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
 
end

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
 value=0;
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

 % Yearly output from the model as a proportion of the population in P at some point during the year for
 % 2013-final year, Estim1 is a column vector
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  

 % Actual proportions of population that were non-addicted prescription opioid users at some point
 % during the year for 2013-2017 
 % (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Data simulated when testing codes 
 Data1=[0.387803312960829;0.385859884749359;0.385291909033797;0.384844182844686;0.384426363863043;0.384001142924548;0.383502489332366;0.383286474912885;0.382780015767370;0.382368169620462;0.381962853339432;0.381583803634593;0.381120645490644;0.380863382568114;0.380449265557906;0.380039206710655;0.379645083482821;0.379293608791355;0.378896962400306;0.378569085690886;0.378225964168557;0.377859610734567;0.377476726514315;0.377114708562721;0.376842862140627];
 
 % The difference between estimated value and data
 Diff1=Estim1-Data1; 


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

 
 % Yearly output from the model as a proportion of population in A at some point during the year for
 % 2013-final year, Estim2 is a column vector
 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population at some point during the year in 2014 and 2015 
 % (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data2=[42000/5651993 48000/5602117];
 
 % Data simulated when testing codes 
 Data2=[0.0117455980073673;0.0126003611353261;0.0135702699068408;0.0145862134466473;0.0156163578770581;0.0166472360917028;0.0176707728614604;0.0186868024270982;0.0196925574245404;0.0206879590347452;0.0216726879014664;0.0226474706126398;0.0236108065439615;0.0245644710781424;0.0255079130690137;0.0264411627138692;0.0273641185151985;0.0282773526420237;0.0291804040459065;0.0300737126980973;0.0309574258792886;0.0318316178524296;0.0326962699927304;0.0335510122628854;0.0343968111832996];
 
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
 

 % Yearly output from the model as a proportion of population in H at some point during the year for
 % 2013-final year, Data3 is a column vector 
 Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);  
 
 % Actual proportion of heroin addicted individuals in the population at some point during the year 
 % in 2014 and 2015
 % (total number of heroin addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes 
 Data3=[0.00103907582372071;0.00101755059193969;0.00100856050223205;0.00100532420926167;0.00100495607550998;0.00100622617515543;0.00100858046437510;0.00101173893007554;0.00101555379888286;0.00101993773451460;0.00102482901368315;0.00103018380897159;0.00103596238440371;0.00104213951408922;0.00104868879606276;0.00105558903032135;0.00106282166296150;0.00107037211192985;0.00107822549298304;0.00108637072936078;0.00109479753199721;0.00110349676404172;0.00111246022459328;0.00112168003192911;0.00113115181543487];
 
 % The difference between estimated value and data
 Diff3=Estim3-Data3;

 
 %Comparing simulated data for susceptibles to output of model for
 %susceptibles 
 State1=y(:,1);
 State_data_1=[0.858700000000000;0.859426382221539;0.858612750962439;0.857680026489993;0.856728627418263;0.855778719000533;0.854981985805504;0.853823981529086;0.853027641683518;0.852136485264308;0.851258558244580;0.850343756535407;0.849577072330486;0.848561235741790;0.847738135336314;0.846924588052323;0.846114111479355;0.845257995125053;0.844477966333902;0.843629583190915;0.842807390205450;0.842020984909400;0.841266458741277;0.840517686175110;0.839661576789117;0.838925559955817];
 
 State_diff_1= State1-State_data_1;
 
 %Comparing simulated data for prescription users to output of model for
 %prescription users
 State2=y(:,2);
 State_data_2=[0.130000000000000;0.128146752382792;0.127850272546011;0.127684209002999;0.127547822881584;0.127420861575431;0.127150984714692;0.127254640305519;0.127005918887905;0.126863057848528;0.126717499541811;0.126619468588702;0.126383153635609;0.126407329040188;0.126248261503800;0.126089783843546;0.125938322363150;0.125842698871172;0.125680620681888;0.125597021425835;0.125496899221430;0.125370537843949;0.125221780321309;0.125076774697224;0.125048976016212;0.124909978261869];
 
 State_diff_2=State2-State_data_2;
 
 
 %Comparing simulated data for opioid addicts to output of model for
 %opioid addicts
 State3=y(:,3);
 State_data_3=[0.0100000000000000;0.0105521986905989;0.0113279484674021;0.0122010261109671;0.0131127408816685;0.0140358756701086;0.0149592278436868;0.0158755538671669;0.0167852882262808;0.0176856622676330;0.0185767718376313;0.0194582876246347;0.0203309911869715;0.0211932556342387;0.0220470263202622;0.0228916103977938;0.0237270657018265;0.0245532718684156;0.0253708274478466;0.0261792051297063;0.0269788861682556;0.0277699903595808;0.0285525798990412;0.0293266223811932;0.0300917205162188;0.0308489038337722];
 
 State_diff_3=State3-State_data_3;
 
 %Comparing simulated data for heroin addicts to output of model for
 %heroin addicts
 State4=y(:,4);
 State_data_4=[0.00100000000000000;0.000955474051020676;0.000936517341675051;0.000928614467030446;0.000925813409980114;0.000925563989978451;0.000926786882742384;0.000928992670475523;0.000931928796634418;0.000935465244291553;0.000939522999751411;0.000944045900759262;0.000948994285772316;0.000954331395841673;0.000960034832376391;0.000966080092296467;0.000972447827651086;0.000979121005419110;0.000986086427265061;0.000993330263648804;0.00100084252305340;0.00100861375781570;0.00101663559133985;0.00102490048518095;0.00103340112048735;0.00104213374400898];
 
 State_diff_4=State4-State_data_4;
 
 
 %Comparing simulated data for stably recovered individuals to output of model for
 %stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.000300000000000000;0.000919192654049952;0.00127251068247289;0.00150612392901015;0.00168499540850452;0.00183897976394859;0.00198101475337481;0.00211683162775210;0.00224922240566137;0.00237932937523934;0.00250764737622560;0.00263444135049720;0.00275978856116059;0.00288384818794147;0.00300654200724719;0.00312793761404071;0.00324805262801645;0.00336691312994073;0.00348449910909832;0.00360085998989521;0.00371598188181106;0.00382987312925404;0.00394254544703346;0.00405401626129209;0.00416432555796432;0.00427342420453181];
 
 State_diff_5=State5-State_data_5;
 
 %%%%%
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
f=zeros(8,1);
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








           
 

