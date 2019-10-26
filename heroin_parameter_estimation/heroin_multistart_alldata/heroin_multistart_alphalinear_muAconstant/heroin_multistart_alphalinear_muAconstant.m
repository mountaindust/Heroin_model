%File name: heroin_multistart_alphalinear_muAconstant.m

clf;
clear all;

% Realistic parameter bounds
%           [m      betaA     betaP   theta1   epsilon  gamma   theta2   sigma    zeta   theta3    nu        b     P0        A0       H0       R0   ]
LowerBounds=[-0.1  0.00001   0.00001   0.05      1       0.005    0.1     0.1    0.0001     10      0.0001   0.1   0.001   0.0001   0.00001  0.00001  ];
UpperBounds=[-0.001  0.01     0.01     0.3       5        0.1     0.4       2       0.2     20       0.2     0.5   0.35      0.01    0.002     0.1    ];



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
numstartpoints=2000;

% Runs MultiStart with numstartpoints to find a solution or multiple local solutions to problem; 
% solutions contains the distinct local minima found during the run
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

% x vector to estimate, alpha=m*t+b 
m=x(1);
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
b=x(12);


pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];



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
P0=x(13);
A0=x(14);
H0=x(15);
R0=x(16);
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
 State_data_1=[0.878100000000000;0.876666546924509;0.876299025258515;0.876456501469636;0.876817216893106;0.877238447639042;0.877672079502912;0.878100829731545;0.878514572010537;0.878906382579047;0.879274368965907;0.879620361360020;0.879947129249510;0.880256909496128;0.880551306364888;0.880831851912734;0.881099618342216;0.881355610096019;0.881600759905256;0.881835897612911;0.882061691029446;0.882278621065846;0.882487397282168;0.882688214984969;0.882881651859619];
 
 
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
 State_data_2=[0.0950000000000000;0.0973363109753449;0.0984508468393681;0.0989117003789873;0.0990651765958730;0.0990738048231760;0.0990011909105158;0.0988767452130333;0.0987199880499770;0.0985459706612804;0.0983627847256185;0.0981733785267092;0.0979787464548160;0.0977796504809884;0.0975768337511818;0.0973709115912034;0.0971623854968742;0.0969515327123125;0.0967385559902748;0.0965236140570888;0.0963068188614997;0.0960882897844022;0.0958681612021066;0.0956465042318401;0.0954234383019021];
 
 
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
 State_data_3=[0.00650000000000000;0.00534025859240171;0.00440211211768753;0.00364564788724283;0.00303462637577537;0.00253958472919088;0.00213726010824393;0.00180924919326778;0.00154114999921931;0.00132041330244590;0.00113735604148011;0.000984487069446751;0.000855957809625031;0.000747209597832733;0.000654816451374683;0.000575698149156979;0.000507764981778466;0.000449379608462006;0.000399178079206152;0.000356010579962034;0.000319003973060903;0.000287446940300732;0.000260365466678502;0.000237366764969711;0.000217632351408772];
 
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
 State_data_4=[0.000400000000000000;0.000422833634091739;0.000448376447745582;0.000477058398459706;0.000509428415831872;0.000546138746324545;0.000587907130170575;0.000635476470834669;0.000689553786618447;0.000750819547836161;0.000819874546724446;0.000897179741241305;0.000983019848273445;0.00107749649311482;0.00118052418109363;0.00129196383627505;0.00141152150616845;0.00153887815595476;0.00167366714685023;0.00181552586057717;0.00196412888244117;0.00211925785800203;0.00228059147864219;0.00244808892886509;0.00262149302544592];
 
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
 State_data_5=[0.0200000000000000;0.0202340498736817;0.0203996393368601;0.0205090918659392;0.0205735517195107;0.0206020240620355;0.0206015623476686;0.0205776993907347;0.0205347361530397;0.0204764139087785;0.0204056157196356;0.0203245933019090;0.0202351466370574;0.0201387339311777;0.0200365192506654;0.0199295745097974;0.0198187096720914;0.0197045994263161;0.0195878388774017;0.0194689518883880;0.0193483572524650;0.0192263843503360;0.0191034845692763;0.0189798250882066;0.0188557844604449];
 
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
 State_data_6=[0;0.0615574214314985;0.122914437620572;0.184120062732634;0.245202161528258;0.306171176146082;0.367028604678507;0.427773407097734;0.488404348282969;0.548920099350697;0.609319166869824;0.669599913391109;0.729760773546408;0.789800370203391;0.849717514326499;0.909511066848130;0.969180052506786;1.02872357159342;1.08814078083398;1.14743088690703;1.20659317654633;1.26562701214414;1.32453165671849;1.38330659537532;1.44195114879440];
 
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
 State_data_7=[0;0.000123541896271180;0.000246555120174957;0.000368108758626559;0.000487337115487758;0.000603446270111838;0.000715696021563742;0.000823396191506881;0.000925951967346226;0.00102282410830137;0.00111359587670225;0.00119802073755247;0.00127604699074029;0.00134781550560593;0.00141365579324808;0.00147394803658958;0.00152922207962071;0.00158002049234712;0.00162692373101272;0.00167049969947400;0.00171126215772119;0.00174959948086405;0.00178602600330057;0.00182072743112496;0.00185414118946404];
 
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
 State_data_8=[0;2.89763175546919e-05;6.10197750857941e-05;9.66061671888524e-05;0.000136335533293563;0.000180919893312983;0.000231147142704379;0.000287841274067736;0.000351802683616083;0.000423812463864439;0.000504582718911253;0.000594694354380050;0.000694557784432277;0.000804403206485998;0.000924273493879356;0.00105415638166838;0.00119388178908844;0.00134324926762895;0.00150200559373835;0.00166989770059057;0.00184670212778268;0.00203229510042487;0.00222645950677153;0.00242923867217349;0.00264047425002330];
 
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
 State_data_9=[0;0.000886299052649950;0.00162004074458737;0.00222786211471056;0.00273396572744100;0.00315799180359563;0.00351539079580603;0.00381832887222964;0.00407630225218005;0.00429757252785865;0.00448851973750013;0.00465413654092766;0.00479842163255471;0.00492461948558844;0.00503532536246295;0.00513293140202255;0.00521925090236949;0.00529579289466275;0.00536390716569781;0.00542478710270918;0.00547939814183631;0.00552850715846617;0.00557314926333864;0.00561372919366153;0.00565112132754567];
 
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
 State_data_10=[0;5.21560241779809e-06;1.07351402757885e-05;1.65975327332267e-05;2.28461758994058e-05;2.95318219505202e-05;3.67140379732337e-05;4.44616747707215e-05;5.28540172456642e-05;6.19765012497958e-05;7.19236389648114e-05;8.27972750223064e-05;9.47041359065530e-05;0.000107753389066177;0.000122054028641461;0.000137713726745808;0.000154837659834688;0.000173526516644733;0.000193876595278011;0.000215980913194933;0.000239926037292871;0.000265792272756848;0.000293668090604526;0.000323624932773257;0.000355746719019541];
 
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
 %Data1=[0.340202161528258;0.342267363350584;0.340076413313416;0.337398025415194;0.334575509536422;0.331664791109567];
 
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
 %Data2=[0.00698733711548776;0.00347324122763384;0.00189124502261338;0.00110913289850545;0.000689805059878946;0.000461883004803746];
 
 
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
 %Data3=[0.000724895566154391;0.00103230888743464;0.00148234385292961];
 
 
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
%Data4=[0.156557421431498;0.158693327164419;0.159656471951430;0.159993799174612;0.160034191213696;0.159931233355602;0.159745993329743;0.159507686398268;0.159235739117706;0.158945038180407;0.158643531246904;0.158334238682008;0.158018343111799;0.157696794604096;0.157370386272813;0.157039897249859;0.156705904583505;0.156368741952872;0.156028662063326;0.155685903696394;0.155340654459305;0.154992934358758;0.154643099858937;0.154291057650916];

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
% 2013-2016 (anyone in A class throughout the year times mu_A)
 
 Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];
     
 % Actual Data for years 2013-2017
 Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
 
 %Testing Data
 %Data5=[0.00273396572744100;0.00134233652473906;0.000722119380374654;0.000420829269814784];
 
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
 set(gca, 'xtick', [ 0 1 2 3])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})                   
                   
 
 % Yearly simulation of individuals overdosing from H class during the year for years
% 2013-2017 (anyone in H class throughout the year times mu_H)
 
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
        y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual Data for years 2013-2017
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 
 %Testing Data
 %Data6=[2.28461758994058e-05;3.00078413462584e-05;4.18501186608887e-05;6.01335239281354e-05;8.50883774581830e-05];
 
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
                   
fprintf('alpha values')

disp(a(0,pars))
disp(a(1,pars))
disp(a(2,pars))
disp(a(3,pars))
disp(a(4,pars))
disp(a(5,pars))
disp(a(6,pars))



function value = HeroinModel_ODE15s(z)

% Parameters
m=z(1);
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
b=z(12);


% Parameter vector
pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Final time N; will run from beginning of 2013 to beginning of 2019 where t=0 represents 2013
% and t=6 represents 2019, with spacing (N-0)/(25-1)=0.25 between the points to represent quarters of a year:
N = 6; 
tspan=linspace(0,N,25);

% Initial conditions
P0=z(13);
A0=z(14);
H0=z(15);
R0=z(16);
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
 %Data1=[0.340202161528258;0.342267363350584;0.340076413313416;0.337398025415194;0.334575509536422;0.331664791109567];
 
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
 %Data2=[0.00698733711548776;0.00347324122763384;0.00189124502261338;0.00110913289850545;0.000689805059878946;0.000461883004803746];
 
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
 %Data3=[0.000724895566154391;0.00103230888743464;0.00148234385292961];
 

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
 %Data4=[0.156557421431498;0.158693327164419;0.159656471951430;0.159993799174612;0.160034191213696;0.159931233355602;0.159745993329743;0.159507686398268;0.159235739117706;0.158945038180407;0.158643531246904;0.158334238682008;0.158018343111799;0.157696794604096;0.157370386272813;0.157039897249859;0.156705904583505;0.156368741952872;0.156028662063326;0.155685903696394;0.155340654459305;0.154992934358758;0.154643099858937;0.154291057650916];
 

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
 %Data5=[0.00273396572744100;0.00134233652473906;0.000722119380374654;0.000420829269814784];
 
 
 %The difference between estimated value and data
 Diff5=Estim5-Data5;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from H class at some point during the year for
 % years 2013-2017, Estim6 is a column vector
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10); y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual proportion of heroin addicts that overdose each year 2013-2017
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 
 %Testing Data
 %Data6=[2.28461758994058e-05;3.00078413462584e-05;4.18501186608887e-05;6.01335239281354e-05;8.50883774581830e-05];
 
 

 %The difference between estimated value and data 
 Diff6=Estim6-Data6;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.878100000000000;0.876666546924509;0.876299025258515;0.876456501469636;0.876817216893106;0.877238447639042;0.877672079502912;0.878100829731545;0.878514572010537;0.878906382579047;0.879274368965907;0.879620361360020;0.879947129249510;0.880256909496128;0.880551306364888;0.880831851912734;0.881099618342216;0.881355610096019;0.881600759905256;0.881835897612911;0.882061691029446;0.882278621065846;0.882487397282168;0.882688214984969;0.882881651859619];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.0950000000000000;0.0973363109753449;0.0984508468393681;0.0989117003789873;0.0990651765958730;0.0990738048231760;0.0990011909105158;0.0988767452130333;0.0987199880499770;0.0985459706612804;0.0983627847256185;0.0981733785267092;0.0979787464548160;0.0977796504809884;0.0975768337511818;0.0973709115912034;0.0971623854968742;0.0969515327123125;0.0967385559902748;0.0965236140570888;0.0963068188614997;0.0960882897844022;0.0958681612021066;0.0956465042318401;0.0954234383019021];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00650000000000000;0.00534025859240171;0.00440211211768753;0.00364564788724283;0.00303462637577537;0.00253958472919088;0.00213726010824393;0.00180924919326778;0.00154114999921931;0.00132041330244590;0.00113735604148011;0.000984487069446751;0.000855957809625031;0.000747209597832733;0.000654816451374683;0.000575698149156979;0.000507764981778466;0.000449379608462006;0.000399178079206152;0.000356010579962034;0.000319003973060903;0.000287446940300732;0.000260365466678502;0.000237366764969711;0.000217632351408772];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.000400000000000000;0.000422833634091739;0.000448376447745582;0.000477058398459706;0.000509428415831872;0.000546138746324545;0.000587907130170575;0.000635476470834669;0.000689553786618447;0.000750819547836161;0.000819874546724446;0.000897179741241305;0.000983019848273445;0.00107749649311482;0.00118052418109363;0.00129196383627505;0.00141152150616845;0.00153887815595476;0.00167366714685023;0.00181552586057717;0.00196412888244117;0.00211925785800203;0.00228059147864219;0.00244808892886509;0.00262149302544592];
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.0200000000000000;0.0202340498736817;0.0203996393368601;0.0205090918659392;0.0205735517195107;0.0206020240620355;0.0206015623476686;0.0205776993907347;0.0205347361530397;0.0204764139087785;0.0204056157196356;0.0203245933019090;0.0202351466370574;0.0201387339311777;0.0200365192506654;0.0199295745097974;0.0198187096720914;0.0197045994263161;0.0195878388774017;0.0194689518883880;0.0193483572524650;0.0192263843503360;0.0191034845692763;0.0189798250882066;0.0188557844604449];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0615574214314985;0.122914437620572;0.184120062732634;0.245202161528258;0.306171176146082;0.367028604678507;0.427773407097734;0.488404348282969;0.548920099350697;0.609319166869824;0.669599913391109;0.729760773546408;0.789800370203391;0.849717514326499;0.909511066848130;0.969180052506786;1.02872357159342;1.08814078083398;1.14743088690703;1.20659317654633;1.26562701214414;1.32453165671849;1.38330659537532;1.44195114879440];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000123541896271180;0.000246555120174957;0.000368108758626559;0.000487337115487758;0.000603446270111838;0.000715696021563742;0.000823396191506881;0.000925951967346226;0.00102282410830137;0.00111359587670225;0.00119802073755247;0.00127604699074029;0.00134781550560593;0.00141365579324808;0.00147394803658958;0.00152922207962071;0.00158002049234712;0.00162692373101272;0.00167049969947400;0.00171126215772119;0.00174959948086405;0.00178602600330057;0.00182072743112496;0.00185414118946404];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;2.89763175546919e-05;6.10197750857941e-05;9.66061671888524e-05;0.000136335533293563;0.000180919893312983;0.000231147142704379;0.000287841274067736;0.000351802683616083;0.000423812463864439;0.000504582718911253;0.000594694354380050;0.000694557784432277;0.000804403206485998;0.000924273493879356;0.00105415638166838;0.00119388178908844;0.00134324926762895;0.00150200559373835;0.00166989770059057;0.00184670212778268;0.00203229510042487;0.00222645950677153;0.00242923867217349;0.00264047425002330];
 
 State_diff_8=State8-State_data_8;
 
 % Comparing simulated data for proportion of individuals that overdose out of A throughout the year and the model output 
 State9=y(:,9);
 State_data_9=[0;0.000886299052649950;0.00162004074458737;0.00222786211471056;0.00273396572744100;0.00315799180359563;0.00351539079580603;0.00381832887222964;0.00407630225218005;0.00429757252785865;0.00448851973750013;0.00465413654092766;0.00479842163255471;0.00492461948558844;0.00503532536246295;0.00513293140202255;0.00521925090236949;0.00529579289466275;0.00536390716569781;0.00542478710270918;0.00547939814183631;0.00552850715846617;0.00557314926333864;0.00561372919366153;0.00565112132754567];
 
 State_diff_9=State9-State_data_9;
 
 % Comparing simulated data for proportion of individuals that overdose out of H the year and the model output 
 State10=y(:,10);
 State_data_10=[0;5.21560241779809e-06;1.07351402757885e-05;1.65975327332267e-05;2.28461758994058e-05;2.95318219505202e-05;3.67140379732337e-05;4.44616747707215e-05;5.28540172456642e-05;6.19765012497958e-05;7.19236389648114e-05;8.27972750223064e-05;9.47041359065530e-05;0.000107753389066177;0.000122054028641461;0.000137713726745808;0.000154837659834688;0.000173526516644733;0.000193876595278011;0.000215980913194933;0.000239926037292871;0.000265792272756848;0.000293668090604526;0.000323624932773257;0.000355746719019541];
 
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

           
function alpha = a(t,pars)
    alpha = pars(1)*t+pars(16);
end


function f = HeroinModel(t,y,pars)
f=zeros(10,1);
f(1)=-a(t,pars)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=a(t,pars)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time;
% i.e. individuals who enter the P class at any time from S (used in
                                                             % Estim1, Estim4)
f(6) = a(t,pars)*y(1);

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

  
