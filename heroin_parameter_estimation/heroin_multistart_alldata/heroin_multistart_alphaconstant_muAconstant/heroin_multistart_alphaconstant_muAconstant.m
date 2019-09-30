%File name: heroin_multistart_alphalinear_muAconstant.m 

clf;
clear all;

% Realistic parameter bounds
%           [alpha  betaA     betaP  theta1  epsilon  gamma   theta2   sigma    zeta   theta3       nu        P0        A0        H0       R0   ]
LowerBounds=[0.1  0.00001   0.00001   0.1      1       0.005    0.1     0.1    0.0001     10      0.0001     0.001   0.0001   0.00001  0.00001  ];
UpperBounds=[0.5  0.01     0.01       4.5      5       0.1       9       2       0.2      19       0.2       0.35      0.01    0.002     0.1    ];
 
%LowerBounds=[ 0.1      0.00001  0.000001  0.00001   0.8    0.001    0.0001  0.0001    0.0001  0.001   0.0001   0.0001   0.00001  0.00001  0.00001  ];
%UpperBounds=[ 0.8        0.01     0.01    0.001      8       0.1       2       1        0.5     4      0.1       0.5      0.1       0.1      0.1   ];
%LowerBounds=[ 0.1     0.0000001 0.00000001  0.00001    0.8    0.001    0.0001  0.0001   0.00001  0.001   0.000001   0.0001   0.00001  0.00001  0.00001  ];
%UpperBounds=[ 0.8        0.001     0.001     0.1        8      0.1       5       5       0.5     10      0.1       0.5      0.1       0.1      0.1   ];
  


% Initial starting points for parameters, starting in the middle of each of the ranges
xstart=0.5*(LowerBounds + UpperBounds); 

% Create MultiStart problem using optimization function fmincon;
% x0 is xstart, objective is what we are trying to minimize which comes from 
% value = HeroinModel_ODE45(z) = fval(x) as output
problem=createOptimProblem('fmincon','objective',@HeroinModel_ODE15s,...
         'x0', xstart,...
         'lb',LowerBounds,...
         'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

% Define a multistart problem; results are reported after each local solver run, in addition to the final summary
ms=MultiStart('Display', 'iter'); 

% Number of times I want to run optimization scheme
numstartpoints=50;

% Runs MultiStart with numstartpoints to find a solution or multiple local solutions to problem; 
% solutions contains the distinct local minima found during the run
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

% x vector to estimate 
alpha=x(1);
beta_A=x(2);
beta_P=x(3);
theta_1=x(4);
epsilon=x(5);
mu=0.00710;  
mu_A=0.00884;
mu_H=0.0466;
gamma=x(6);   
theta_2=x(7);
sigma=x(8);
zeta=x(9);
theta_3=x(10);
nu=x(11);
omega=0.0000000001;

pars=[alpha,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega];



%syms alpha(t)
%alpha(t)=piecewise(t<3,pars(1)*t+pars(16),t>=3,pars(17)*t+pars(18));


% Print optimal parameter solution and objective function value in command
% window when completed 
format short 
x
fval

% Final time N; will run from beginning of 2013 to beginning of 2019 where t=0 represents 2013
% and t=6 represents 2019, with spacing (N-0)/(25-1)=0.25 between the points to represent quarters of a year:
N = 6; 
tspan=linspace(0,N,25);


% Initial conditions
P0=x(12);
A0=x(13);
H0=x(14);%0.0001;%0.00136;%x(7);
R0=x(15);
S0=1-P0-A0-H0-R0;
X0=0;
L0=0;
M0=0;
J0=0;
K0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0;J0;K0];

% Run stiff ODE solver 
[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

% Gives solution vector for each state 
  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);
  J=y(:,9);
  K=y(:,10);
 
  
 % Making sure S+P+A+H+R=1
 total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);
 
 %%% For testing purposes: states and corresponding simulated data 
 %%% For testing purposes: states and corresponding simulated data 
 State1=y(:,1);
 State_data_1=[0.893950000000000;0.906260600617363;0.912545816733428;0.915724783803624;0.917366389706360;0.918220782107385;0.918645570140307;0.918863357868658;0.918976778036070;0.919031780140936;0.919051256511692;0.919053758411835;0.919045648338514;0.919035975765080;0.919022474851329;0.919006134856080;0.918986663881597;0.918964126869529;0.918938870967058;0.918911056747263;0.918880655644427;0.918847520725198;0.918811740095758;0.918772864809202;0.918730759236667];
 
 
 % Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1),'k-','LineWidth',3)
 %plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State2=y(:,2);
 State_data_2=[0.100000000000000;0.0876765941095081;0.0813811988978871;0.0781929848513978;0.0765421167897582;0.0756779442808696;0.0752425685752547;0.0750131739731454;0.0748869651307160;0.0748178482978492;0.0747827914269986;0.0747631025595947;0.0747522722173598;0.0747411023148207;0.0747316973714157;0.0747228980213334;0.0747148219315794;0.0747072109451840;0.0746995179317213;0.0746913866260949;0.0746826467515788;0.0746732190086565;0.0746628203581209;0.0746516662815590;0.0746396783408485];
 
 % Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2),'b-','LineWidth',3)
 %plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State3=y(:,3);
 State_data_3=[0.00550000000000000;0.00519732432494098;0.00497205474401894;0.00480305005343622;0.00467390050465430;0.00457282943867184;0.00449157121905707;0.00442371082304013;0.00436469418294693;0.00431125947866790;0.00426112678315408;0.00421233723960143;0.00416366432587812;0.00411375446836516;0.00406205869108483;0.00400797607409889;0.00395103492610914;0.00389101538172291;0.00382772907533032;0.00376095542909859;0.00369050839394453;0.00361650485921811;0.00353871571667588;0.00345740892422991;0.00337258987022563];
 
 % Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3),'r-','LineWidth',3)
 %plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State4=y(:,4);
 State_data_4=[0.000300000000000000;0.000317084315026625;0.000337755929226056;0.000362399371169064;0.000391204766321878;0.000424297063757189;0.000461800642558676;0.000503874648307842;0.000550701168117925;0.000602491688627867;0.000659493010754619;0.000721963838700273;0.000790205905141539;0.000864511574176533;0.000945211364422046;0.00103264517106098;0.00112714434320039;0.00122899389540499;0.00133844916020935;0.00145578172364143;0.00158124318987831;0.00171498167083110;0.00185726813417492;0.00200811617166065;0.00216768734810337];
 
 % Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
 %plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 

 State5=y(:,5);
 State_data_5=[0.000250000000000000;0.000548396633370740;0.000763173696281766;0.000916781921102061;0.00102638823303348;0.00110414710905707;0.00115848942244782;0.00119588268627553;0.00122086148130251;0.00123662039286863;0.00124533226624262;0.00124883794906389;0.00124820921187674;0.00124465587630515;0.00123855772046318;0.00123034587610976;0.00122033491616810;0.00120865290675648;0.00119543286421359;0.00118081947236558;0.00116494601855533;0.00114777373439557;0.00112945569347732;0.00110994381147619;0.00108928520222351];
 
 % Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5),'Color', [0.7,0,0.7],'LineWidth',3)
 %plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State6=y(:,6);
 State_data_6=[0;0.0450317135860872;0.0905164160227824;0.136234471548033;0.182067821556238;0.227959997542694;0.273884056841525;0.319823209084668;0.365769684194975;0.411719970817370;0.457672285000904;0.503625165901536;0.549578059722995;0.595530210887578;0.641481701304694;0.687432362060003;0.733382116891694;0.779330847516952;0.825278401202380;0.871224645631923;0.917169455557892;0.963112674676869;1.00905417612518;1.05499380026506;1.10093140335687];
 
 % Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6),'LineWidth',3)
 %plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend({'Proportion that enter P at some point during the year'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State7=y(:,7);
 State_data_7=[0;0.000120503932005548;0.000297760231461671;0.000515744624047440;0.000761939032019806;0.00102728245935279;0.00130534222165395;0.00159106452679277;0.00188079015357472;0.00217182520541625;0.00246220326860219;0.00275019438692853;0.00303460140330786;0.00331413374617644;0.00358808557697787;0.00385570253450357;0.00411632693699396;0.00436948996442154;0.00461474986305352;0.00485163919115730;0.00507972112605774;0.00529879892658896;0.00550839979514543;0.00570848351989035;0.00589877603177667];
 
 % Simulated data for L and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7),'LineWidth',3)
 %plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend({'Proportion that enter A at some point during the year'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State8=y(:,8);
 State_data_8=[0;2.12324412384357e-05;4.63046825944511e-05;7.56523790654478e-05;0.000109521007318534;0.000148092719669448;0.000191550462444680;0.000240113655186178;0.000294027003690504;0.000353567314609661;0.000419050535307124;0.000490806494931895;0.000569213960526294;0.000654642161223710;0.000747506737025057;0.000848236332382209;0.000957253086573106;0.00107494106892904;0.00120165942649867;0.00133778101957834;0.00148366077531215;0.00163956540883621;0.00180586686433111;0.00198270204214003;0.00217034446337051];
 
 % Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8),'LineWidth',3)
 %plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend({'Proportion that enter H at some point during the year'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) %correspond to the actual t values from t vector that I want to label
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})

 
 State9=y(:,9);
 State_data_9=[0;1.18060035477295e-05;2.30338050685866e-05;3.38259234060662e-05;4.42906459950823e-05;5.45034224142530e-05;6.45150583111077e-05;7.43631913230769e-05;8.40724322193780e-05;9.36582602200049e-05;0.000103129050289736;0.000112490814690644;0.000121743999768370;0.000130889991731897;0.000139924046362428;0.000148841530933083;0.000157636908096722;0.000166302907737293;0.000174832189647631;0.000183217573098077;0.000191451810372029;0.000199526063783094;0.000207433246671944;0.000215164837490734;0.000222712912088113];
 
 % Simulated data for J and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 plot(t,y(:,9),'LineWidth',3)
 %plot(t(1:end), State_data_9, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('J(t)')
 legend({'Proportion that overdose from A'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) %correspond to the actual t values from t vector that I want to label
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})


 State10=y(:,10);
 State_data_10=[0;3.59298664397514e-06;7.40468228362448e-06;1.14793711499641e-05;1.58649968476876e-05;2.06111071658900e-05;2.57684313517117e-05;3.13891769631262e-05;3.75275824877416e-05;4.42402262262404e-05;5.15870008254115e-05;5.96295126518397e-05;6.84344866395628e-05;7.80685010534169e-05;8.86052855682655e-05;0.000100121712227842;0.000112696420845752;0.000126415209788056;0.000141368000176557;0.000157642512845790;0.000175328243153766;0.000194527921775159;0.000215329011630933;0.000237838319048378;0.000262152817321765];
 
 % Simulated data for K and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 plot(t,y(:,10),'LineWidth',3)
 %plot(t(1:end), State_data_10, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('K(t)')
 legend({'Proportion that overdose from H'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) %correspond to the actual t values from t vector that I want to label
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})

 
 %%% Data points we are interested in 


 % Yearly simulation of individuals in P class at all during the year for years
 % 2013-2017
 
 Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6); y(21,2)+y(25,6)-y(21,6)];
     
 % Actual Data for years 2013-2018
 Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
 
 %Testing Data 
 %Data1=[0.282067821556238;0.260243979428496;0.258695340658736;0.258556329386059;0.258502160597777;0.258444594550555];
 
 % Data points from proportion that is in P at some point in the year and corresponding ODE solution points 
 figure(11)
 hold all
 z1 = linspace(0,5,6); %defines mesh where going to plot Estim1, Data1 values 
 scatter(z1, Estim1, 100,'o');
 scatter(z1, Data1, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 5 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 

 % Yearly simulation of individuals in A class at all during the year for years
 % 2013-2017
 Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7); y(21,3)+y(25,7)-y(21,7)];
 
 
 % Actual Data for years 2013-2018
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
 
 %Testing Data
 %Data2=[0.00626193903201981;0.00579275162620922;0.00551850543268008;0.00524538985956421;0.00491442911517291;0.00450956329966346];
 
 
  % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
 figure(12)
 hold all
 z2 = linspace(0,5,6);
 scatter(z2, Estim2, 100,'o');
 scatter(z2, Data2, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14) 
 set(gca, 'xtick', [ 0 1 2 3 4 5 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017', '2018'})



 % Yearly simulation of individuals in H class at all during the year for years
 % 2014-2016
 
 Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8); y(13,4)+y(17,8)-y(13,8)];
 
 % Actual Data for years 2014-2016
 Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 
 %Testing Data
 %Data3=[0.000575710762693847;0.000825888124953715;0.00117824503118835];
 
 
 % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
 figure(13)
 hold all
 z3 = linspace(0,2,3);
 scatter(z3, Estim3, 100,'o');
 scatter(z3, Data3, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})

 
 
 % Yearly simulation of individuals in P class at all during the quarters in years
 % 2013-2018
 Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);

        
 %Actual Data for years 2013-2018 
 Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
        833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
        827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
        832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
        775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
        688451./5754509; 683498./5754509; 641894./5754509; 625054./5754509];
 
 %Testing Data
 %Data4=[0.145031713586087;0.133161296546203;0.127099254423138;0.124026334859603;0.122434292776214;0.121602003579701;0.121181720818398;0.120959649083453;0.120837251753111;0.120770162481383;0.120735672327631;0.120715996381054;0.120704423381943;0.120692592731936;0.120682358126725;0.120672652853024;0.120663552556837;0.120654764630612;0.120645762361264;0.120636196552063;0.120625865870556;0.120614720456963;0.120602444498004;0.120589269373368];
 
 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(14)
 hold all
 z4 = linspace(0,23,24);
 scatter(z4, Estim4, 100,'o');
 scatter(z4, Data4, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in P') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 23])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
                   
                   
                   
                   
% Yearly simulation of individuals overdosing from A class during the year for years
% 2013-2017 (anyone in A class throughout the year times mu_A)
 
 Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];
     
 % Actual Data for years 2013-2017
 Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
 
 %Testing Data
 %Data5=[4.42906459950823e-05;3.97817862242957e-05;3.76715675489923e-05;3.58929083283522e-05];
 
 
  % Data points from proportion that is in A at some point and overdoses in the year and corresponding ODE solution points 
 figure(15)
 hold all
 z5 = linspace(0,3,4); %defines mesh where going to plot Estim5, Data5 values 
 scatter(z5, Estim5, 100,'o');
 scatter(z5, Data5, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})                   
                   
 
 % Yearly simulation of individuals overdosing from H class during the year for years
% 2013-2017 (anyone in H class throughout the year times mu_H)
 
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
        y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual Data for years 2013-2017
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 
 %Testing Data
 %Data6=[1.58649968476876e-05;2.16625856400539e-05;3.09069041518212e-05;4.42619342061892e-05;6.26318223080137e-05];
 
  % Data points from proportion that is in H at some point and overdoses in the year and corresponding ODE solution points 
 figure(16)
 hold all
 z6 = linspace(0,4,5); %defines mesh where going to plot Estim6, Data6 values 
 scatter(z6, Estim6, 100,'o');
 scatter(z6, Data6, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from H') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})                   
    

function value = HeroinModel_ODE15s(z)

% Parameters
alpha=z(1);
beta_A=z(2);
beta_P=z(3);
theta_1=z(4);
epsilon=z(5);
mu=0.00710;  
mu_A=0.00884;
mu_H=0.0466;
gamma=z(6);   
theta_2=z(7);
sigma=z(8);
zeta=z(9);
theta_3=z(10);
nu=z(11);
omega=0.0000000001;

% Parameter vector
pars=[alpha,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega];

% Final time N; will run from beginning of 2013 to beginning of 2019 where t=0 represents 2013
% and t=6 represents 2019, with spacing (N-0)/(25-1)=0.25 between the points to represent quarters of a year:
N = 6; 
tspan=linspace(0,N,25);

% Initial conditions
P0=z(12);
A0=z(13);
H0=z(14);
R0=z(15);
S0=1-P0-A0-H0-R0;
X0=0;
L0=0;
M0=0;
J0=0;
K0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0;J0;K0];


% Run stiff ODE solver 
[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);
  
 
 % COMPARING MODEL ESTIMATES TO DATA 
 value=0;
 %%%%%
 % In order to count the total number of individuals in P at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
  
 % To get the output from the model of the proportion of non-addicted
 % prescription opioid users in 2013 (first entry of Estim1), 
 % we take the total number of non-addiction prescription opioid users 
 % at the beginning of 2013 (P_0=IC)
 % and add on the number of individuals that enter the P class at any point during the year 2013, 
 % which comes from integrating ODE X'=dy(6) from t=0 to t=1; this gives
 % first value in Estim1. 
 
 % To get the output from the model of the proportion of non-addicted prescription opioid users in 2014, 2015
 % 2016, and 2017 (remaining entries of Estim1),
 % we take the initial number of non-addicted prescription opioid users at
 % the beginning of the year
 % in 2014: y(5,2), 2015: y(9,2), 2016: y(13,2), and 2017: y(17,2) and add the number of individuals that enter
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
 % 2013-2017, Estim1 is a column vector

 Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6); y(21,2)+y(25,6)-y(21,6)];
 


 % Actual proportions of population that were non-addicted prescription opioid users at some point
 % during the year for 2013-2017 
 % (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 
 Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];

 %Testing Data
 %Data1=[0.282067821556238;0.260243979428496;0.258695340658736;0.258556329386059;0.258502160597777;0.258444594550555];
 
 % The difference between estimated values and data
 
 Diff1=Estim1-Data1; 


 %%%%%
 % In order to count the total number of individuals in A at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of opioid
 % addicts in 2013 (first entry of Estim2), 
 % we take the total number of opioid addicts
 % at the beginning of 2013 (A_0=IC)
 % and add on the number of individuals that enter the P class at any point during the year 2013, 
 % which comes from integrating ODE L'=dy(7) from t=0 to t=1; this gives
 % first value in Estim2. 
 
 % To get the output from the model of the proportion of opioid addicts in 2014, 2015, 2016, and 2017 
 % (remaining entries of Estim2),
 % we take the initial number of opioid addicts at the beginning of the year
 % in 2014: y(5,3), 2015: y(9,3), 2016: y(13,3), and 2017: y(17,3)
 % and add the number of individuals that enter
 % the A class at any point during the year, which comes from
 % integrating ODE L'=dy(9) but just focusing in on these specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 % for 2016, we have to subtract because integrating gives total number of new cases from t=0 to t=4, so have to 
 % subtract off the number from t=0 to t=3. 
 % for 2017, we have to subtract because integrating gives total number of new cases from t=0 to t=5, so have to 
 % subtract off the number from t=0 to t=4. 

 
 % Yearly output from the model as a proportion of population in A at some point during the year for
 % 2013-2017, Estim2 is a column vector
 
 Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7); y(21,3)+y(25,7)-y(21,7)];
    
 
 % Actual proportions of population that were opioid addicted individuals in
 % the population at some point during the year 
 % (total number of opioid addicted individuals in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
 
 %Testing Data
 %Data2=[0.00626193903201981;0.00579275162620922;0.00551850543268008;0.00524538985956421;0.00491442911517291;0.00450956329966346];
 
 % The difference between estimated value and data
 
 Diff2=Estim2-Data2; 


 
 %%%%%
 % In order to count the total number of individuals in H at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of heroin/fentanyl addicts in
 % 2014, 2015, and 2016 (Estim3), we take initial number of heroin/fentanyl addicts in 2014: y(5,4), 
 % 2015: y(9,4), and 2016: y(13,4), and add the number of individuals that enter the H class at any point
 % during the year 2014, 2015, or 2016, which comes from
 % integrating ODE M'=dy(8) but just focusing in on the three specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1. 
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 % for 2016, we have to subtract because integrating gives total number of new cases from t=0 to t=4, so have to 
 % subtract off the number from t=0 to t=3. 
 

 % Yearly output from the model as a proportion of population in H at some point during the year for
 % 2014-2016, Estim3 is a column vector 

 Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8); y(13,4)+y(17,8)-y(13,8)];
 
 
 % Actual proportion of heroin addicted individuals in the population at some point during the year 
 % in 2014, 2015, and 2016
 % (total number of heroin addicted individuals in 2014, 2015, and 2016 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 
 Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 
 %Testing Data
 %Data3=[0.000575710762693847;0.000825888124953715;0.00117824503118835];
 

 % The difference between estimated value and data
 
 Diff3=Estim3-Data3;

 
  %%%%%
 % In order to count the total number of individuals in P at some point throughout a certain quarter of a year, 
 % we need to count the number who are in the class AT ALL during the quarter,
 % even if they leave or come back at some point. 
  
 % To get the output from the model of the proportion of non-addicted
 % prescription opioid users in the first quarter of 2013 (first entry of Estim4), 
 % we take the total number of non-addiction prescription opioid users 
 % at the beginning of 2013 (P_0=IC)
 % and add on the number of individuals that enter the P class at any point during the first quarter of 2013, 
 % which comes from integrating ODE X'=dy(6) from t=0 to t=0.25; this gives
 % first value in Estim1. 
 
 % To get the output from the model of the proportion of non-addicted prescription opioid users in the 
 % rest of the quarters of years 2013-2018 (remaining entries of Estim4),
 % we take the initial number of non-addicted prescription opioid users at
 % the beginning of each of the quarters and add the number of individuals that enter
 % the P class at any point during each of the quarters, which comes from
 % integrating ODE X'=dy(6) but just focusing in oon the specific quarter
 % for 2013 Q2, we have to subtract because integrating gives total number of new cases from t=0 to t=.5, so have to 
 % subtract off the number from t=0 to t=0.25;
 % for 2013 Q3, we have to subtract because integrating gives total number of new cases from t=0 to t=0.75, so have to 
 % subtract off the number from t=0 to t=0.5. 
 % for 2013 Q4, we have to subtract because integrating gives total number of new cases from t=0 to t=1, so have to 
 % subtract off the number from t=0 to t=0.75. 
 % for 2014 Q1: gives total number of new cases from t=0 to t=1.25, so subtract the number from t=0 to t=1. 
 % for 2014 Q2: gives total number of new cases from t=0 to t=1.5, so subtract the number from t=0 to t=1.25. 
 % for 2014 Q3: gives total number of new cases from t=0 to t=1.75, so subtract the number from t=0 to t=1.5. 
 % for 2014 Q4: gives total number of new cases from t=0 to t=2, so subtract the number from t=0 to t=1.75. 
 % for 2015 Q1: gives total number of new cases from t=0 to t=2.25, so subtract the number from t=0 to t=2. 
 % for 2015 Q2: gives total number of new cases from t=0 to t=2.5, so subtract the number from t=0 to t=2.25. 
 % for 2015 Q3: gives total number of new cases from t=0 to t=2.75, so subtract the number from t=0 to t=2.5. 
 % for 2015 Q4: gives total number of new cases from t=0 to t=3, so subtract the number from t=0 to t=2.75. 
 % for 2016 Q1: gives total number of new cases from t=0 to t=3.25, so subtract the number from t=0 to t=3. 
 % for 2016 Q2: gives total number of new cases from t=0 to t=3.5, so subtract the number from t=0 to t=3.25. 
 % for 2016 Q3: gives total number of new cases from t=0 to t=3.75, so subtract the number from t=0 to t=3.5. 
 % for 2016 Q4: gives total number of new cases from t=0 to t=4, so subtract the number from t=0 to t=3.75. 
 % for 2017 Q1: gives total number of new cases from t=0 to t=4.25, so subtract the number from t=0 to t=4. 
 % for 2017 Q2: gives total number of new cases from t=0 to t=4.5, so subtract the number from t=0 to t=4.25. 
 % for 2017 Q3: gives total number of new cases from t=0 to t=4.75, so subtract the number from t=0 to t=4.5. 
 % for 2017 Q4: gives total number of new cases from t=0 to t=5, so subtract the number from t=0 to t=4.75. 
 % for 2018 Q1: gives total number of new cases from t=0 to t=5.25, so subtract the number from t=0 to t=5. 
 % for 2018 Q2: gives total number of new cases from t=0 to t=5.5, so subtract the number from t=0 to t=5.25. 
 % for 2018 Q3: gives total number of new cases from t=0 to t=5.75, so subtract the number from t=0 to t=5.5. 
 % for 2018 Q4: gives total number of new cases from t=0 to t=6, so subtract the number from t=0 to t=5.75. 
        
 
 % Yearly output from the model as a proportion of the population in P at some point during a quarter of year for
 % years 2013-2018, Estim4 is a column vector

 Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);



 %Actual proportion of non-addicted prescription opioid users for every quarter in years 2013-2018 
 Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
        833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
        827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
        832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
        775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
        688451./5754509; 683498./5754509; 641894./5754509; 625054./5754509];
      
 %Testing Data
 %Data4=[0.145031713586087;0.133161296546203;0.127099254423138;0.124026334859603;0.122434292776214;0.121602003579701;0.121181720818398;0.120959649083453;0.120837251753111;0.120770162481383;0.120735672327631;0.120715996381054;0.120704423381943;0.120692592731936;0.120682358126725;0.120672652853024;0.120663552556837;0.120654764630612;0.120645762361264;0.120636196552063;0.120625865870556;0.120614720456963;0.120602444498004;0.120589269373368];
 

 %The difference between estimated value and data 
 Diff4=Estim4-Data4;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from A class at some point during the year for
 % years 2013-2017, Estim5 is a column vector
 Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];
    
 % Actual proportion of addicted prescription opioid users that overdose
 % each year 2013-2017
 
 Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
 
 %Testing Data
 %Data5=[4.42906459950823e-05;3.97817862242957e-05;3.76715675489923e-05;3.58929083283522e-05];
 
 
 %The difference between estimated value and data
 Diff5=Estim5-Data5;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from H class at some point during the year for
 % years 2013-2017, Estim6 is a column vector
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10); y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual proportion of heroin addicts that overdose each year 2013-2017
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 
 %Testing Data
 %Data6=[1.58649968476876e-05;2.16625856400539e-05;3.09069041518212e-05;4.42619342061892e-05;6.26318223080137e-05];
 
 

 %The difference between estimated value and data 
 Diff6=Estim6-Data6;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.893950000000000;0.906260600617363;0.912545816733428;0.915724783803624;0.917366389706360;0.918220782107385;0.918645570140307;0.918863357868658;0.918976778036070;0.919031780140936;0.919051256511692;0.919053758411835;0.919045648338514;0.919035975765080;0.919022474851329;0.919006134856080;0.918986663881597;0.918964126869529;0.918938870967058;0.918911056747263;0.918880655644427;0.918847520725198;0.918811740095758;0.918772864809202;0.918730759236667];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.100000000000000;0.0876765941095081;0.0813811988978871;0.0781929848513978;0.0765421167897582;0.0756779442808696;0.0752425685752547;0.0750131739731454;0.0748869651307160;0.0748178482978492;0.0747827914269986;0.0747631025595947;0.0747522722173598;0.0747411023148207;0.0747316973714157;0.0747228980213334;0.0747148219315794;0.0747072109451840;0.0746995179317213;0.0746913866260949;0.0746826467515788;0.0746732190086565;0.0746628203581209;0.0746516662815590;0.0746396783408485];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00550000000000000;0.00519732432494098;0.00497205474401894;0.00480305005343622;0.00467390050465430;0.00457282943867184;0.00449157121905707;0.00442371082304013;0.00436469418294693;0.00431125947866790;0.00426112678315408;0.00421233723960143;0.00416366432587812;0.00411375446836516;0.00406205869108483;0.00400797607409889;0.00395103492610914;0.00389101538172291;0.00382772907533032;0.00376095542909859;0.00369050839394453;0.00361650485921811;0.00353871571667588;0.00345740892422991;0.00337258987022563];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.000300000000000000;0.000317084315026625;0.000337755929226056;0.000362399371169064;0.000391204766321878;0.000424297063757189;0.000461800642558676;0.000503874648307842;0.000550701168117925;0.000602491688627867;0.000659493010754619;0.000721963838700273;0.000790205905141539;0.000864511574176533;0.000945211364422046;0.00103264517106098;0.00112714434320039;0.00122899389540499;0.00133844916020935;0.00145578172364143;0.00158124318987831;0.00171498167083110;0.00185726813417492;0.00200811617166065;0.00216768734810337];
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.000250000000000000;0.000548396633370740;0.000763173696281766;0.000916781921102061;0.00102638823303348;0.00110414710905707;0.00115848942244782;0.00119588268627553;0.00122086148130251;0.00123662039286863;0.00124533226624262;0.00124883794906389;0.00124820921187674;0.00124465587630515;0.00123855772046318;0.00123034587610976;0.00122033491616810;0.00120865290675648;0.00119543286421359;0.00118081947236558;0.00116494601855533;0.00114777373439557;0.00112945569347732;0.00110994381147619;0.00108928520222351];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0450317135860872;0.0905164160227824;0.136234471548033;0.182067821556238;0.227959997542694;0.273884056841525;0.319823209084668;0.365769684194975;0.411719970817370;0.457672285000904;0.503625165901536;0.549578059722995;0.595530210887578;0.641481701304694;0.687432362060003;0.733382116891694;0.779330847516952;0.825278401202380;0.871224645631923;0.917169455557892;0.963112674676869;1.00905417612518;1.05499380026506;1.10093140335687];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000120503932005548;0.000297760231461671;0.000515744624047440;0.000761939032019806;0.00102728245935279;0.00130534222165395;0.00159106452679277;0.00188079015357472;0.00217182520541625;0.00246220326860219;0.00275019438692853;0.00303460140330786;0.00331413374617644;0.00358808557697787;0.00385570253450357;0.00411632693699396;0.00436948996442154;0.00461474986305352;0.00485163919115730;0.00507972112605774;0.00529879892658896;0.00550839979514543;0.00570848351989035;0.00589877603177667];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;2.12324412384357e-05;4.63046825944511e-05;7.56523790654478e-05;0.000109521007318534;0.000148092719669448;0.000191550462444680;0.000240113655186178;0.000294027003690504;0.000353567314609661;0.000419050535307124;0.000490806494931895;0.000569213960526294;0.000654642161223710;0.000747506737025057;0.000848236332382209;0.000957253086573106;0.00107494106892904;0.00120165942649867;0.00133778101957834;0.00148366077531215;0.00163956540883621;0.00180586686433111;0.00198270204214003;0.00217034446337051];
 
 State_diff_8=State8-State_data_8;
 
 % Comparing simulated data for proportion of individuals that overdose out of A throughout the year and the model output 
 State9=y(:,9);
 State_data_9=[0;1.18060035477295e-05;2.30338050685866e-05;3.38259234060662e-05;4.42906459950823e-05;5.45034224142530e-05;6.45150583111077e-05;7.43631913230769e-05;8.40724322193780e-05;9.36582602200049e-05;0.000103129050289736;0.000112490814690644;0.000121743999768370;0.000130889991731897;0.000139924046362428;0.000148841530933083;0.000157636908096722;0.000166302907737293;0.000174832189647631;0.000183217573098077;0.000191451810372029;0.000199526063783094;0.000207433246671944;0.000215164837490734;0.000222712912088113];
 
 State_diff_9=State9-State_data_9;
 
 % Comparing simulated data for proportion of individuals that overdose out of H the year and the model output 
 State10=y(:,10);
 State_data_10=[0;3.59298664397514e-06;7.40468228362448e-06;1.14793711499641e-05;1.58649968476876e-05;2.06111071658900e-05;2.57684313517117e-05;3.13891769631262e-05;3.75275824877416e-05;4.42402262262404e-05;5.15870008254115e-05;5.96295126518397e-05;6.84344866395628e-05;7.80685010534169e-05;8.86052855682655e-05;0.000100121712227842;0.000112696420845752;0.000126415209788056;0.000141368000176557;0.000157642512845790;0.000175328243153766;0.000194527921775159;0.000215329011630933;0.000237838319048378;0.000262152817321765];
 
 State_diff_10=State10-State_data_10;
 
 %%%%%
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors [norm gives sqrt(sum from 1 to N
 % of (Diff#)^2)]
 % normalized by norm of the data [because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; gives least squares percentage
 % error so each piece weighted evenly].
 
 
 % Objective function value we wish to minimize; want value=fval(x) to be small  when run MultiStart
 
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 
 %%%%
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4);
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5);

 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff6,2)./norm(Data6);
 
 %value=norm(Diff5,2)./norm(Data5);
 
 % For testing purposes with states and data sets
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8)+norm(State_diff_9,2)./norm(State_data_9)+norm(State_diff_10,2)./norm(State_data_10);
  
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8);
 
 
end

          

function f = HeroinModel(t,y,pars)
f=zeros(10,1);
f(1)=-pars(1)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=pars(1)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time;
% i.e. individuals who enter the P class at any time from S (used in Estim1, Estim4)
f(6) = pars(1)*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
% i.e. individuals who enter the A class at any time (used in Estim2)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction over time; 
% i.e. individuals who enter the H class at any time (used in Estim3)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));

%J' ODE to calculate number of prescription opioid addict overdoses over
%time; i.e. individuals who overdose at any time (used in Estim5)
f(9) = pars(7)*y(3);

%K' ODE to calculate number of heroin/fentanayl addict overdoses over
%time; i.e. individuals who overdose at any time (used in Estim6)
f(10) = pars(8)*y(4);
end

