%File name: heroin_multistart_final.m 

clf

% We wish to estimate the parameter vector 
% x =[alpha,beta_A,beta_P,theta_1,epsilon,gamma,theta_2,sigma,zeta,theta_3,nu]
%Give ranges on each of the parameters 
%LowerBounds=[0.01  0.00001 0.0001 0.00001 0.8 0.001 0.00001 0.01 0.01 0.00001 0.01];
%UpperBounds=[0.7    0.2     0.009   0.1    4   0.1    0.3     2   1     0.6   1  ];

LowerBounds=[0.00001  0.00001  0.00001 0.00001 0.00001 0.00001 0.00001 0.00001 0.00001 0.00001 0.00001 ];
UpperBounds=[2  2  2  2  4  2  2  2  2  2  2 ];

%Initial starting points for parameters, starting in the middle of each of the ranges
xstart=0.5*(LowerBounds + UpperBounds); 

% Create MultiStart problem using optimization function fmincon;
% x0 is xstart, objective is what we are trying to minimize which comes from 
% value = HeroinModel_ODE45(z) = fval(x) as output
problem=createOptimProblem('fmincon','objective',@HeroinModel_ODE45,...
         'x0', xstart,...
         'lb',LowerBounds,...
         'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

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


[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

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
 State_data_1=[0.858700000000000;0.859382451177455;0.858604607539590;0.857660390456534;0.856728664053505;0.855831462756203;0.854818448510738;0.853966835015160;0.853041562066129;0.852114347630275;0.851249398614380;0.850427766329290;0.849413233153236;0.848670385354288;0.847793611075849;0.846900699142143;0.846067150196388;0.845276414532071;0.844482605275038;0.843580460110938;0.842912609498823;0.842028155481626;0.841177707588313;0.840465739401827;0.839720285363812;0.838934387878802];
 
 %Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1))
 plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 legend('Proportion of susceptibles simulated','Proportion of susceptibles (simulated) data' )

 State2=y(:,2);
 State_data_2=[0.130000000000000;0.128190859549039;0.127858447724637;0.127703922719689;0.127547785245572;0.127367904999871;0.127315177685713;0.127111212503539;0.126991941779825;0.126885283206232;0.126726694822195;0.126535122369615;0.126547645752294;0.126297743058357;0.126192563903323;0.126113766710047;0.125985469095978;0.125824205168737;0.125675962416061;0.125646337897149;0.125391262375665;0.125363337961644;0.125310880773010;0.125128924989896;0.124990034793451;0.124901114240573];
 
 
 %Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2))
 plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 legend('Proportion of prescription users simulated','Proportion of prescription users (simulated) data' )


 State3=y(:,3);
 State_data_3=[0.0100000000000000;0.0105520810379251;0.0113279580377919;0.0122011851371189;0.0131126407121267;0.0140360059502953;0.0149586483662231;0.0158761481158622;0.0167853312635482;0.0176855631381937;0.0185767307602861;0.0194586430063721;0.0203302972695245;0.0211937186551860;0.0220472616701950;0.0228915099342311;0.0237268679526164;0.0245533504786638;0.0253708475997679;0.0261789989397853;0.0269793300987756;0.0277700211359088;0.0285522077321395;0.0293264053319200;0.0300919676699944;0.0308489418129934];
 
 
 %Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3))
 plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend('Proportion of opioid addicts simulated','Proportion of opioid addicts (simulated) data' )

 
 State4=y(:,4);
 State_data_4=[0.00100000000000000;0.000955479358429735;0.000936520540829022;0.000928630846892882;0.000925807292459248;0.000925559336853647;0.000926791423597724;0.000928993889769114;0.000931928300308530;0.000935464908638258;0.000939522970371164;0.000944046932237310;0.000948992764027780;0.000954332716854014;0.000960035584757600;0.000966080053462942;0.000972447580934136;0.000979121388931770;0.000986086672997716;0.000993330033986466;0.00100084366197055;0.00100861402760792;0.00101663504777262;0.00102490027248497;0.00103340182097105;0.00104213403539852];
 
 %Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4))
 plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend('Proportion of heroin/fentanyl addicts simulated','Proportion of heroin/fentanyl addicts (simulated) data' )
 

 State5=y(:,5);
 State_data_5=[0.000300000000000000;0.000919128877151769;0.00127246615715200;0.00150587083976547;0.00168510269633713;0.00183906695677702;0.00198093401372812;0.00211681047566992;0.00224923659018886;0.00237934111666128;0.00250765283276720;0.00263442136248605;0.00275983106091751;0.00288382021531443;0.00300652776587515;0.00312794416011594;0.00324806517408285;0.00336690843159594;0.00348449803613539;0.00360087301814071;0.00371595436476510;0.00382987139321337;0.00394256885876428;0.00405403000387223;0.00416431035177222;0.00427342203223374];
 
  %Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5))
 plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 legend('Proportion of stably recovered addicts simulated','Proportion of stably recovered addicts (simulated) data' )


 
 State6=y(:,6);
 State_data_6=[0;0.257803312960829;0.515516445327395;0.772958081815182;1.03011805565687;1.28699659663833;1.54357687798745;1.79992838260512;2.05596021721248;2.31173431409195;2.56723942586388;2.82248477966150;3.07744911470740;3.33218660656243;3.58664266009036;3.84084366414446;4.09479308701157;4.34849984813124;4.60195075805143;4.85516709976985;5.10813916403490;5.36086822898202;5.61335730187264;5.86561224806565;6.11765018193115;6.36944406805556];
 
 %Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6))
 plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend('X ODE solution','data simulated' )


 State7=y(:,7);
 State_data_7=[0;0.00174559800736727;0.00379376045209450;0.00603608189153314;0.00842126922721341;0.0109248862226030;0.0135362466441972;0.0162477916619708;0.0190590402219021;0.0219663094201617;0.0249686061872739;0.0280645222511090;0.0312537052391140;0.0345335205961041;0.0379047360400078;0.0413656227887593;0.0449151751048347;0.0485522279182066;0.0522763086918147;0.0560858852898746;0.0599803928582655;0.0639589325692985;0.0680205600621473;0.0721642501558365;0.0763886400375288;0.0806937307046096];
 
 %Simulated data for X and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7))
 plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend('L ODE solution','data simulated' )
  
 State8=y(:,8);
 State_data_8=[0;3.90758237207137e-05;0.000101152364639723;0.000173195525196723;0.000249905267427941;0.000329047932957803;0.000409710118134776;0.000491503699767489;0.000574249959367504;0.000657874961615947;0.000742347451838990;0.000827653465770727;0.000913791373983054;0.00100075947261445;0.00108856759086200;0.00117722155454837;0.00126673049257325;0.00135710432788366;0.00144835543439440;0.00154049450011238;0.00163353496582435;0.00172748997476816;0.00182237298099418;0.00191819761424761;0.00201497716099577;0.00211272785594329];
  
 %Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8))
 plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend('M ODE solution','data simulated' )
 
 
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 Data1=[0.387803312960829;0.385859884749359;0.385291909033797;0.384844182844686;0.384426363863043;0.384001142924548;0.383502489332366;0.383286474912885;0.382780015767370;0.382368169620462;0.381962853339432;0.381583803634593;0.381120645490644;0.380863382568114;0.380449265557906;0.380039206710655;0.379645083482821;0.379293608791355;0.378896962400306;0.378569085690886;0.378225964168557;0.377859610734567;0.377476726514315;0.377114708562721;0.376842862140627];
 
 %Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
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
 figure(10)
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
 figure(11)
 hold all
 plot(t(1:end-1),Estim3)
 plot(t(1:end-1), Data3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('Data points interested in', 'ODE solution')
 
 

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
initials = [S0;P0;A0;H0;R0;X0;L0;M0];


[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);
  
 
 
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
 % in 2014, 2015, and 2016
 % (total number of heroin addicted individuals in 2014, 2015, and 2016 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes 

 Data3=[0.00103907582372071;0.00101755059193969;0.00100856050223205;0.00100532420926167;0.00100495607550998;0.00100622617515543;0.00100858046437510;0.00101173893007554;0.00101555379888286;0.00101993773451460;0.00102482901368315;0.00103018380897159;0.00103596238440371;0.00104213951408922;0.00104868879606276;0.00105558903032135;0.00106282166296150;0.00107037211192985;0.00107822549298304;0.00108637072936078;0.00109479753199721;0.00110349676404172;0.00111246022459328;0.00112168003192911;0.00113115181543487];
 
 % The difference between estimated value and data
 
 Diff3=Estim3-Data3;

 
 %Comparing simulated data for susceptibles to output of model for
 %susceptibles 
 State1=y(:,1);
 State_data_1=[0.858700000000000;0.859382451177455;0.858604607539590;0.857660390456534;0.856728664053505;0.855831462756203;0.854818448510738;0.853966835015160;0.853041562066129;0.852114347630275;0.851249398614380;0.850427766329290;0.849413233153236;0.848670385354288;0.847793611075849;0.846900699142143;0.846067150196388;0.845276414532071;0.844482605275038;0.843580460110938;0.842912609498823;0.842028155481626;0.841177707588313;0.840465739401827;0.839720285363812;0.838934387878802];
 
 State_diff_1= State1-State_data_1;
 
 %Comparing simulated data for prescription users to output of model for
 %prescription users
 State2=y(:,2);
 State_data_2=[0.130000000000000;0.128190859549039;0.127858447724637;0.127703922719689;0.127547785245572;0.127367904999871;0.127315177685713;0.127111212503539;0.126991941779825;0.126885283206232;0.126726694822195;0.126535122369615;0.126547645752294;0.126297743058357;0.126192563903323;0.126113766710047;0.125985469095978;0.125824205168737;0.125675962416061;0.125646337897149;0.125391262375665;0.125363337961644;0.125310880773010;0.125128924989896;0.124990034793451;0.124901114240573];
 
 State_diff_2=State2-State_data_2;
 
 
 %Comparing simulated data for opioid addicts to output of model for
 %opioid addicts
 State3=y(:,3);
 State_data_3=[0.0100000000000000;0.0105520810379251;0.0113279580377919;0.0122011851371189;0.0131126407121267;0.0140360059502953;0.0149586483662231;0.0158761481158622;0.0167853312635482;0.0176855631381937;0.0185767307602861;0.0194586430063721;0.0203302972695245;0.0211937186551860;0.0220472616701950;0.0228915099342311;0.0237268679526164;0.0245533504786638;0.0253708475997679;0.0261789989397853;0.0269793300987756;0.0277700211359088;0.0285522077321395;0.0293264053319200;0.0300919676699944;0.0308489418129934];
 
 State_diff_3=State3-State_data_3;
 
 %Comparing simulated data for heroin addicts to output of model for
 %heroin addicts
 State4=y(:,4);
 State_data_4=[0.00100000000000000;0.000955479358429735;0.000936520540829022;0.000928630846892882;0.000925807292459248;0.000925559336853647;0.000926791423597724;0.000928993889769114;0.000931928300308530;0.000935464908638258;0.000939522970371164;0.000944046932237310;0.000948992764027780;0.000954332716854014;0.000960035584757600;0.000966080053462942;0.000972447580934136;0.000979121388931770;0.000986086672997716;0.000993330033986466;0.00100084366197055;0.00100861402760792;0.00101663504777262;0.00102490027248497;0.00103340182097105;0.00104213403539852];
 
 State_diff_4=State4-State_data_4;
 
 
 %Comparing simulated data for stably recovered individuals to output of model for
 %stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.000300000000000000;0.000919128877151769;0.00127246615715200;0.00150587083976547;0.00168510269633713;0.00183906695677702;0.00198093401372812;0.00211681047566992;0.00224923659018886;0.00237934111666128;0.00250765283276720;0.00263442136248605;0.00275983106091751;0.00288382021531443;0.00300652776587515;0.00312794416011594;0.00324806517408285;0.00336690843159594;0.00348449803613539;0.00360087301814071;0.00371595436476510;0.00382987139321337;0.00394256885876428;0.00405403000387223;0.00416431035177222;0.00427342203223374];
 
 State_diff_5=State5-State_data_5;
 
 %Comparing simulated data for 6th ODE and the model output 
 State6=y(:,6);
 State_data_6=[0;0.257803312960829;0.515516445327395;0.772958081815182;1.03011805565687;1.28699659663833;1.54357687798745;1.79992838260512;2.05596021721248;2.31173431409195;2.56723942586388;2.82248477966150;3.07744911470740;3.33218660656243;3.58664266009036;3.84084366414446;4.09479308701157;4.34849984813124;4.60195075805143;4.85516709976985;5.10813916403490;5.36086822898202;5.61335730187264;5.86561224806565;6.11765018193115;6.36944406805556];
 
 State_diff_6=State6-State_data_6;
 
 %Comparing simulated data for 7th ODE and the model output 
 State7=y(:,7);
 State_data_7=[0;0.00174559800736727;0.00379376045209450;0.00603608189153314;0.00842126922721341;0.0109248862226030;0.0135362466441972;0.0162477916619708;0.0190590402219021;0.0219663094201617;0.0249686061872739;0.0280645222511090;0.0312537052391140;0.0345335205961041;0.0379047360400078;0.0413656227887593;0.0449151751048347;0.0485522279182066;0.0522763086918147;0.0560858852898746;0.0599803928582655;0.0639589325692985;0.0680205600621473;0.0721642501558365;0.0763886400375288;0.0806937307046096];
 
 State_diff_7=State7-State_data_7;
 
 %Comparing simulated data for 8th ODE and the model output 
 State8=y(:,8);
 State_data_8=[0;3.90758237207137e-05;0.000101152364639723;0.000173195525196723;0.000249905267427941;0.000329047932957803;0.000409710118134776;0.000491503699767489;0.000574249959367504;0.000657874961615947;0.000742347451838990;0.000827653465770727;0.000913791373983054;0.00100075947261445;0.00108856759086200;0.00117722155454837;0.00126673049257325;0.00135710432788366;0.00144835543439440;0.00154049450011238;0.00163353496582435;0.00172748997476816;0.00182237298099418;0.00191819761424761;0.00201497716099577;0.00211272785594329];
 
 State_diff_8=State8-State_data_8;
 
 %%%%%
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sqrt(sum from 1 to N of (diff#)^2)
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; 
 % gives least squares percentage error so each piece weighted evenly)
 
 %value=norm(State_diff_1,2)+norm(State_diff_2,2)+norm(State_diff_3,2)+norm(State_diff_4,2)+norm(State_diff_5,2);
 
 %Run1
 %value = norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5);
 
 %Run2 and Run3
 %value = norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8);
 
 %value=norm(Diff1,2)+norm(Diff2,2)+norm(Diff3,2);
 
 %Run4
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);
 
 %value=sum(Diff1.^2)./sum(Data1.^2)+sum(Diff2.^2)./sum(Data2.^2)+sum(Diff3.^2)./sum(Data3.^2);
 
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








           
 

