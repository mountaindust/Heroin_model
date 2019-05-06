%File name: heroin_multistart_quarters_no2016heroin.m 

clf;
clear all;


%Run1
%LowerBounds=[-0.1   0.000001     0.8      0.0001      0.0001    0.1   0.00001  0.00001  0.00001 0.00001  0.00001 0.001];
%UpperBounds=[0.1      0.001       8         0.1         0.1      0.8     0.5       0.1     0.01   0.01      0.9   0.9 ];
 
%Run2, Run3
%LowerBounds=[-0.1    0.00001     0.8      0.00001      0.0001    0.1    0.0001  0.00001  0.00001 0.00001   0.001 0.001];
%UpperBounds=[0.1      0.001       8        0.001         0.1      0.8     0.5       0.1     0.01   0.01       3     3 ];
 
%Run4, Run5
%LowerBounds=[-0.1    0.00001     0.8      0.001      0.0001    0.1    0.0001  0.00001  0.00001 0.00001   0.001  0.001];
%UpperBounds=[0.1      0.001       8        0.1         0.1     0.8     0.5       0.1     0.01   0.01       3      3 ];
 
%Run6
%LowerBounds=[-0.1    0.00001     0.8      0.001      0.0001    0.1    0.0001  0.00001  0.00001 0.00001   0.001 0.001  0.00001  0.00001];
%UpperBounds=[0.1      0.001       8        0.1         0.1      0.8     0.5       0.1     0.01   0.01       3     3       .1     .1];

%Run7--compare results with Run7c bounds if using this code 
LowerBounds=[-0.1    0.00001     0.8      0.0001      0.0001    0.1    0.0001  0.00001  0.00001  0.00001   0.001  0.001  0.00001  0.00001];
UpperBounds=[0.1      0.001       8        0.1           1        0.8     0.5       0.1     0.01     0.01       3     5    .1      .1];
  
%Run7b
%LowerBounds=[-0.1    0.00001     0.8      0.0001      0.00001    0.1    0.0001  0.00001  0.00001 0.00001   0.001 0.001  0.0001  0.0001];
%UpperBounds=[0.1      0.001       8        0.1         0.1       0.8     0.5       0.1     0.01   0.01       3     3      .1      .1];
 
%Run7c CONVERGES EXCEPT FIRST RUN--compare to other bounds if go with this
%code result 
%LowerBounds=[-0.1     0.0001      1      0.0001      0.0001    0.1    0.0001   0.0001   0.0001 0.0001     0.1  0.1   0.001  0.001];
%UpperBounds=[0.1      0.001       8       0.01        0.01     0.8      0.1       0.1     0.01   0.01       3     3      .1      .1];
 
%Run8
%LowerBounds=[-0.1    0.00001     0.8      0.0001      0.0001    0.1    0.0001  0.00001  0.00001 0.00001   0.001 0.001  0.0001  0.0001  0.00001  0.00001];
%UpperBounds=[0.1      0.001       8        0.1         0.1      0.8     0.5       0.1     0.01   0.01       3     3      .1      .1       0.01     0.01];
 

%Bounds for testing code
%LowerBounds=[-0.1     0.00001     0.8      0.0001      0.0001    0.1     0.00001  0.00001  0.00001 0.00001  0.0001   0.001];
%UpperBounds=[0.1       0.001      8         0.01        0.01      0.8     0.5       0.5     0.5    0.5       0.01     0.1 ];
 

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
numstartpoints=200;

% Runs MultiStart with numstartpoints to find a solution or multiple local solutions to problem; 
% solutions contains the distinct local minima found during the run
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

% x vector to estimate, alpha=m*t+b 
m=x(1);
beta_A=0.000273; 
beta_P=0.000777; 
theta_1=x(2);
epsilon=x(3);
mu=0.00868;  
mu_A=0.00870;   
mu_H=0.0507;
gamma=x(4);   
theta_2=x(11);%3*x(2);
sigma=x(5);
zeta=x(13);%0.0214;
theta_3=x(12);%16*x(2);
nu=x(14);%0.0155;
omega=0.0000000001;
b=x(6);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Print optimal parameter solution and objective function value in command
% window when completed 
format short 
x
fval

% Final time N; will run 2013-2018 where t=0 represents 2013
% and t=5 represents 2018, with spacing (N-0)/(25-1)=0.25 between the points:
N = 6; 
tspan=linspace(0,N,25);


% Initial conditions
P0=x(7);
A0=x(8);
H0=x(9);%0.0001;%0.00136;%x(7);
R0=x(10);
S0=1-P0-A0-H0-R0;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];

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
  %alpha=-pars(1)*t+pars(16);
  alpha=x(1)*t+x(6);
  
 % Making sure S+P+A+H+R=1
 total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);
 
 %%% For testing purposes: states and corresponding simulated data 
 %%% For testing purposes: states and corresponding simulated data 
 State1=y(:,1);
 State_data_1=[0.786760000000000;0.849096559775662;0.876641405482671;0.889186631768797;0.895079633937603;0.898073725651136;0.899792310676189;0.900979175329111;0.901922738227621;0.902744120907157;0.903525703988465;0.904298545363606;0.905066920146102;0.905831287155018;0.906593936240071;0.907358806007248;0.908125439496451;0.908894275003860;0.909666942933548;0.910442540242388;0.911221812863727;0.912004452126591;0.912790400006395;0.913579629672869;0.914372166786769];
 
 
 % Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1))
 %plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 legend('Proportion of susceptibles')%,'Proportion of susceptibles (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State2=y(:,2);
 State_data_2=[0.200000000000000;0.137391758114905;0.109661536807078;0.0969695223523026;0.0909479684349754;0.0878343855433950;0.0860013822580893;0.0847033481545710;0.0836511448173080;0.0827233486064567;0.0818373752443332;0.0809620909849510;0.0800932051887634;0.0792302525269369;0.0783709340706020;0.0775112965923273;0.0766517993152429;0.0757920001393280;0.0749302593788612;0.0740674851667191;0.0732029195065194;0.0723368760609821;0.0714694064578173;0.0706005345694132;0.0697302375639648];
 
 
 % Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2))
 %plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 legend('Proportion of prescription users')%,'Proportion of prescription users (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State3=y(:,3);
 State_data_3=[0.00760000000000000;0.00786444322401689;0.00804130390249546;0.00817861219365153;0.00829686308111957;0.00840530823284513;0.00850797931347344;0.00860673107272667;0.00870230694617434;0.00879502985091389;0.00888512351792990;0.00897268287277491;0.00905774190425214;0.00914032530040650;0.00922046009603785;0.00929819883398219;0.00937353299267190;0.00944648268218157;0.00951709193878415;0.00958533675200990;0.00965129538262477;0.00971493544592971;0.00977630240460125;0.00983541728034773;0.00989225996778463];
 
 % Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3))
%plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend('Proportion of opioid addicts')%,'Proportion of opioid addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State4=y(:,4);
 State_data_4=[0.00121000000000000;0.00118907553511728;0.00116847176811839;0.00114820623707514;0.00112828274821752;0.00110869987774823;0.00108945382934144;0.00107053991421861;0.00105195298489348;0.00103368780826489;0.00101573916468520;0.000998101856828860;0.000980770718142271;0.000963740647779967;0.000947006627497517;0.000930558102595570;0.000914396821188789;0.000898518254046546;0.000882912748057731;0.000867585559043041;0.000852518389733417;0.000837718811630932;0.000823176100578988;0.000808885297782603;0.000794851145885819];
 
 % Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4))
 %plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend('Proportion of heroin/fentanyl addicts')%,'Proportion of heroin/fentanyl addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 

 State5=y(:,5);
 State_data_5=[0.00443000000000000;0.00445816335040032;0.00448728203826519;0.00451702744562279;0.00454725179510159;0.00457788069153352;0.00460887391937061;0.00464020552540860;0.00467185701975483;0.00470381282284363;0.00473605808015047;0.00476857891718795;0.00480136203783277;0.00483439436476723;0.00486766296061053;0.00490114045861107;0.00493483136918280;0.00496872391532794;0.00500279299551623;0.00503705227465459;0.00507145385224523;0.00510601754975535;0.00514071502552731;0.00517553317452408;0.00521048453054189];
 
 % Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5))
 %plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 legend('Proportion of stably recovered addicts')%,'Proportion of stably recovered addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State6=y(:,6);
 State_data_6=[0;0.0613467709794854;0.125216002053204;0.189797148258713;0.254333735443525;0.318487743791413;0.382114062932921;0.445145758933416;0.507557333109014;0.569337519637458;0.630478088807434;0.690975090007418;0.750827452758148;0.810034490460835;0.868595139726926;0.926506845120101;0.983769962822696;1.04038353249203;1.09634539445447;1.15165672459813;1.20631343122552;1.26031720750166;1.31366562477631;1.36635756005405;1.41839408737959];
 
 % Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6))
 %plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend('Proportion that enter P at some point during the year')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 State7=y(:,7);
 State_data_7=[0;0.000339557118480829;0.000593618523869692;0.000809621916849858;0.00100779586570128;0.00119725794083937;0.00138196729145419;0.00156373182144740;0.00174326152162490;0.00192085104237566;0.00209669711603550;0.00227086937494662;0.00244337749460406;0.00261422218516008;0.00278340669390784;0.00295093301885284;0.00311680094430453;0.00328100871194106;0.00344355378182944;0.00360443749585012;0.00376364876754363;0.00392119218454272;0.00407705999500611;0.00423124861452216;0.00438376147126013];
 
 % Simulated data for L and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7))
 %plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend('Proportion that enter A at some point during the year')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State8=y(:,8);
 State_data_8=[0;1.53029752325150e-06;2.99264583798051e-06;4.41068275227863e-06;5.79459867040075e-06;7.14935514801869e-06;8.47750036062795e-06;9.78060923714028e-06;1.10597045031806e-05;1.23156068646828e-05;1.35490849447435e-05;1.47608267204845e-05;1.59514455887548e-05;1.71215220664733e-05;1.82716165258114e-05;1.94028686952344e-05;2.05150858790120e-05;2.16087210822438e-05;2.26847531079636e-05;2.37426485366423e-05;2.47841875401369e-05;2.58086334743688e-05;2.68169802646106e-05;2.78096874968083e-05;2.87863156180955e-05];
 
 % Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8))
 %plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend('Proportion that enter H at some point during the year')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) %correspond to the actual t values from t vector that I want to label
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 %%%
 
 
 %%% Data points we are interested in 


 % Yearly simulation of individuals in P class throughout year for years
 % 2013-2017
 
 Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6)];
     
 % Actual Data for years 2013-2017
 Data1=[1825910./5519417; 1805325./5559702; 1800613./5602187; 1744766./5648259; 1620951./5702475];
 
 %Data1=[0.454333735443525;0.344171566100465;0.313035715253311;0.299195267718064;0.285283575660588];
 
 % Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 z1 = linspace(0,4,5); %defines mesh where going to plot Estim1, Data1 values 
 scatter(z1, Estim1,'o');
 scatter(z1, Data1,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
 

 % Yearly simulation of individuals in A class throughout year for years
 % 2013-2017
 Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7)];
 
 
 % Actual Data for years 2013-2017 
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34816./5702475];
 
 %Data2=[0.00860779586570128;0.00903232873704320;0.00940242291915349;0.00973116535395262;0.0100203808159110];
 
 % Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 z2 = linspace(0,4,5);
 scatter(z2, Estim2,'o');
 scatter(z2, Data2,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017'})



 % Yearly simulation of individuals in H class throughout year for years
 % 2014-2016
 
 Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8)];
 
 % Actual Data for years 2014-2015
 Data3=[7560./5559702; 7560./5602187];
 %Data3=[14000./5559006; 14000./5602117];
 %Data3=[0.00113354785405030;0.00105684472597906];
 
  
 % Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(11)
 hold all
 z3 = linspace(0,1,2);
 scatter(z3, Estim3,'o');
 scatter(z3, Data3,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015'})

 
 tspan=linspace(0,N,25);
 Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);

 
 % Actual Data for years 2013-2018 INCLUDES addicts 
 %Data4=[856000./5517176; 870000./5517176; 874000./5517176; 856000./5517176;...
        %842000./5559006; 860000./5559006; 871000./5559006; 850000./5559006;...
        %836000./5602117; 861000./5602117; 865000./5602117; 854000./5602117;...
        %840000./5651993; 829000./5651993; 801000./5651993; 783000./5651993;...
        %783000./5708586; 772000./5708586; 747000./5708586; 713000./5708586;...
        %695000./5779971; 690000./5779971; 648000./5779971; 631000./5779971];
        
 %Actual Data for years 2013-2018 EXCLUDES addicts (take ratio of Rx users
 %(includes those addicted)
 %in quarter out of entire year and assume *same* ratio for opioid addicts in
 %the quarter to year so take out that number, i.e. .46*OA of year)
 Data4=[(856000-20142)./5519417; (870000-20472)./5519417; (874000-20566)./5519417; (856000-20142)./5519417;...
        (842000-19813)./5559702; (860000-20236)./5559702; (871000-20495)./5559702; (850000-20001)./5559702;...
        (836000-19672)./5602187; (861000-20260)./5602187; (865000-20354)./5602187; (854000-20095)./5602187;...
        (840000-17867)./5648259; (829000-17633)./5648259; (801000-17037)./5648259; (783000-16654)./5648259;...
        (783000-16659)./5702475; (772000-16425)./5702475; (747000-15893)./5702475; (713000-15170)./5702475;...
        (695000-15000)./5754509; (690000-15000)./5754509; (648000-15000)./5754509; (631000-15000)./5754509];
 
    %Data4=[0.261346770979485;0.201260989188623;0.174242683012588;0.161506109537114;0.155101976782864;0.151460704684903;0.149033078258584;0.147114922330169;0.145431331345751;0.143863917776433;0.142334376444317;0.140814453735682;0.139300242891450;0.137790901793028;0.136282639463777;0.134774414294923;0.133265368984576;0.131753862101771;0.130241589522517;0.128724191794107;0.127206695782661;0.125685293335629;0.124161341735557;0.122637061894953];
 
 
 figure(12)
 hold all
 z4 = linspace(0,23,24);
 scatter(z4, Estim4,'o');
 scatter(z4, Data4,'x');
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in P each quarter')
 legend('ODE solution', 'Data')
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
 
 
 %{                  
 Estim5=y(1:20,3)+y(2:21,7)-y(1:20,7);                 
                      
 %Actual Data for years 2013-2017
 Data5=[836027./(42.*5517176); 843027./(42.*5517176); 850027./(42.*5517176); 830027./(42.*5517176);...
        810038./(42.*5559006); 840038./(42.*5559006); 851038./(42.*5559006); 820038./(42.*5559006);...
        795305./(42.*5602117); 835305./(42.*5602117); 836305./(42.*5602117); 830305./(42.*5602117);...
        802767./(47.*5651993); 792767./(47.*5651993); 787767./(47.*5651993); 782767./(47.*5651993);...
        769633./(47.*5708586); 759633./(47.*5708586); 729633./(47.*5708586); 692633./(47.*5708586)];

 figure(13)
 hold all
 plot(t(1:20),Estim5, 'o')
 plot(t(1:20), Data5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in A each quarter')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 23])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017'})
                   
                   
                   
Estim6=y(5:16,4)+y(6:17,8)-y(5:16,8);                 
                      
 %Actual Data for years 2014-2016
Data6=[810038./(239.*5559006); 840038./(239.*5559006); 851038./(239.*5559006); 820038./(239.*5559006);...
        795305./(238.*5602117); 835305./(238.*5602117); 836305./(238.*5602117); 830305./(238.*5602117);...
        802767./(170.*5651993); 792767./(170.*5651993); 787767./(170.*5651993); 782767./(170.*5651993)];

figure(14)
 hold all
 plot(t(5:16),Estim6, 'o')
 plot(t(5:16), Data6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in H each quarter')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 23])
 set(gca,'xticklabel',{'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016'})
 
 
 %}
function value = HeroinModel_ODE15s(z)

%Parameters
m=z(1);
beta_A=0.000273; 
beta_P=0.000777; 
theta_1=z(2);
epsilon=z(3);
mu=0.00868;  
mu_A=0.00870;   
mu_H=0.0507;
gamma=z(4);   
theta_2=z(11);%3*z(2);
sigma=z(5);
zeta=z(13);%0.0214;
theta_3=z(12);%16*z(2);
nu=z(14);%0.0155;
omega=0.0000000001;
b=z(6);


pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Final time N; will run 2013-2018 where t=0 represents 2013
% and t=5 represents 2018, with spacing (T-0)/((N+1)-1)=1 between the points
N = 6; 
tspan=linspace(0,N,25);



% Initial conditions
P0=z(7);
A0=z(8);
H0=z(9);
R0=z(10);
S0=1-P0-A0-H0-R0;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];


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
 % 2013-2017, Estim1 is a column vector

Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6)];
 


 % Actual proportions of population (updated 3/12/19) that were non-addicted prescription opioid users at some point
 % during the year for 2013-2017 
 % (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 
 Data1=[1825910./5519417; 1805325./5559702; 1800613./5602187; 1744766./5648259; 1620951./5702475];
 
 %Data1=[0.454333735443525;0.344171566100465;0.313035715253311;0.299195267718064;0.285283575660588];
 
 % Data simulated when testing codes 
 %Data1=[0.399384466780766;0.476721593432771;0.469765087997695;0.449866817308604;0.427737148099884];
 
 % The difference between estimated values and data
 
 Diff1=Estim1-Data1; 


 %%%%%
 % In order to count the total number of individuals in A at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of opioid addicts in 2014 and 2015 (Estim2),
 % we take the initial number of opioid addicts in 2014, y(2,3), and 2015, y(3,3),
 % and add the number of individuals that enter
 % the A class at any point during the year 2014 or 2015, which comes from
 % integrating ODE L'=dy(9) but just focusing in on these two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 

 
 % Yearly output from the model as a proportion of population in A at some point during the year for
 % 2013-2017, Estim2 is a column vector
 
 Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7)];
 % When testing all points with simulated data
 %Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 % Actual proportions of population (updated 3/12/19) that were opioid addicted individuals in
 % the population at some point during the year in 2014 and 2015 
 % (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34816./5702475];
 
 %Data2=[0.00860779586570128;0.00903232873704320;0.00940242291915349;0.00973116535395262;0.0100203808159110];
 
 % Data simulated when testing codes  
 %Data2=[0.00709261474856600;0.0106675766257930;0.0130723402928730;0.0148284456654410;0.0162165492691030];
 
 % The difference between estimated value and data
 
 Diff2=Estim2-Data2; 


 
 %%%%%
 % In order to count the total number of individuals in H at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of heroin/fentanyl addicts in
 % 2014, 2015, and 2016 (Estim3), we take initial number of heroin/fentanyl addicts in 2014, y(2,4), 
 % and 2015, y(3,4), and 2016, y(4,4), and add the number of individuals that enter the H class at any point
 % during the year 2014, 2015, or 2016, which comes from
 % integrating ODE M'=dy(10) but just focusing in on the three specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1. 
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 % for 2016, we have to subtract because integrating gives total number of new cases from t=0 to t=4, so have to 
 % subtract off the number from t=0 to t=3. 
 

 % Yearly output from the model as a proportion of population in H at some point during the year for
 % 2014-2016, Data3 is a column vector 
 % Only for 2014 and 2015 
 %Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8);  

  
 Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8)];
 

 
 % Actual proportion (updated 3/12/19) of heroin addicted individuals in the population at some point during the year 
 % in 2014, 2015, and 2016
 % (total number of heroin addicted individuals in 2014, 2015, and 2016 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 
 Data3=[7560./5559702; 7560./5602187];
 %Data3=[14000./5559006; 14000./5602117];
 %Data3=[0.00113354785405030;0.00105684472597906];
 
 % Data simulated when testing codes 
 %Data3=[0.00116527288223448;0.00120952017524577;0.00118883157707289];
 
 % The difference between estimated value and data
 
 Diff3=Estim3-Data3;

 
  
 Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);
 %Data1=[0.399384466780766;0.476721593432771;0.469765087997695;0.449866817308604;0.427737148099884];
 %FILL IN 2018 total pop/data points and COMPARE TO YEARS that we know and
 %use 2nd graph for 2014 on!!
 
 %Actual Data for years 2013-2018 INCLUDES addicts 
 %Data4=[856000./5517176; 870000./5517176; 874000./5517176; 856000./5517176;...
        %842000./5559006; 860000./5559006; 871000./5559006; 850000./5559006;...
        %836000./5602117; 861000./5602117; 865000./5602117; 854000./5602117;...
        %840000./5651993; 829000./5651993; 801000./5651993; 783000./5651993;...
        %783000./5708586; 772000./5708586; 747000./5708586; 713000./5708586;...
        %695000./5779971; 690000./5779971; 648000./5779971; 631000./5779971];
        
 %Actual Data for years 2013-2018 EXCLUDES addicts 
 Data4=[(856000-20142)./5519417; (870000-20472)./5519417; (874000-20566)./5519417; (856000-20142)./5519417;...
        (842000-19813)./5559702; (860000-20236)./5559702; (871000-20495)./5559702; (850000-20001)./5559702;...
        (836000-19672)./5602187; (861000-20260)./5602187; (865000-20354)./5602187; (854000-20095)./5602187;...
        (840000-17867)./5648259; (829000-17633)./5648259; (801000-17037)./5648259; (783000-16654)./5648259;...
        (783000-16659)./5702475; (772000-16425)./5702475; (747000-15893)./5702475; (713000-15170)./5702475;...
        (695000-15000)./5754509; (690000-15000)./5754509; (648000-15000)./5754509; (631000-15000)./5754509];
 
 %Data4=[0.261346770979485;0.201260989188623;0.174242683012588;0.161506109537114;0.155101976782864;0.151460704684903;0.149033078258584;0.147114922330169;0.145431331345751;0.143863917776433;0.142334376444317;0.140814453735682;0.139300242891450;0.137790901793028;0.136282639463777;0.134774414294923;0.133265368984576;0.131753862101771;0.130241589522517;0.128724191794107;0.127206695782661;0.125685293335629;0.124161341735557;0.122637061894953];
 
 Diff4=Estim4-Data4;
 
%{ 
                  
Estim5=y(1:20,3)+y(2:21,7)-y(1:20,7);                 
                      
 %Actual Data for years 2013-2017
Data5=[836027./(42.*5517176); 843027./(42.*5517176); 850027./(42.*5517176); 830027./(42.*5517176);...
        810038./(42.*5559006); 840038./(42.*5559006); 851038./(42.*5559006); 820038./(42.*5559006);...
        795305./(42.*5602117); 835305./(42.*5602117); 836305./(42.*5602117); 830305./(42.*5602117);...
        802767./(47.*5651993); 792767./(47.*5651993); 787767./(47.*5651993); 782767./(47.*5651993);...
        769633./(47.*5708586); 759633./(47.*5708586); 729633./(47.*5708586); 692633./(47.*5708586)];
                   
Diff5=Estim5-Data5;
                   
 
Estim6=y(5:16,4)+y(6:17,8)-y(5:16,8);                 
                      
 %Actual Data for years 2014-2016
Data6=[810038./(239.*5559006); 840038./(239.*5559006); 851038./(239.*5559006); 820038./(239.*5559006);...
        795305./(238.*5602117); 835305./(238.*5602117); 836305./(238.*5602117); 830305./(238.*5602117);...
        802767./(170.*5651993); 792767./(170.*5651993); 787767./(170.*5651993); 782767./(170.*5651993)];
 
 
Diff6=Estim6-Data6;
 
 
%}
 

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.786760000000000;0.849096559775662;0.876641405482671;0.889186631768797;0.895079633937603;0.898073725651136;0.899792310676189;0.900979175329111;0.901922738227621;0.902744120907157;0.903525703988465;0.904298545363606;0.905066920146102;0.905831287155018;0.906593936240071;0.907358806007248;0.908125439496451;0.908894275003860;0.909666942933548;0.910442540242388;0.911221812863727;0.912004452126591;0.912790400006395;0.913579629672869;0.914372166786769];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.200000000000000;0.137391758114905;0.109661536807078;0.0969695223523026;0.0909479684349754;0.0878343855433950;0.0860013822580893;0.0847033481545710;0.0836511448173080;0.0827233486064567;0.0818373752443332;0.0809620909849510;0.0800932051887634;0.0792302525269369;0.0783709340706020;0.0775112965923273;0.0766517993152429;0.0757920001393280;0.0749302593788612;0.0740674851667191;0.0732029195065194;0.0723368760609821;0.0714694064578173;0.0706005345694132;0.0697302375639648];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00760000000000000;0.00786444322401689;0.00804130390249546;0.00817861219365153;0.00829686308111957;0.00840530823284513;0.00850797931347344;0.00860673107272667;0.00870230694617434;0.00879502985091389;0.00888512351792990;0.00897268287277491;0.00905774190425214;0.00914032530040650;0.00922046009603785;0.00929819883398219;0.00937353299267190;0.00944648268218157;0.00951709193878415;0.00958533675200990;0.00965129538262477;0.00971493544592971;0.00977630240460125;0.00983541728034773;0.00989225996778463];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.00121000000000000;0.00118907553511728;0.00116847176811839;0.00114820623707514;0.00112828274821752;0.00110869987774823;0.00108945382934144;0.00107053991421861;0.00105195298489348;0.00103368780826489;0.00101573916468520;0.000998101856828860;0.000980770718142271;0.000963740647779967;0.000947006627497517;0.000930558102595570;0.000914396821188789;0.000898518254046546;0.000882912748057731;0.000867585559043041;0.000852518389733417;0.000837718811630932;0.000823176100578988;0.000808885297782603;0.000794851145885819];
 
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.00443000000000000;0.00445816335040032;0.00448728203826519;0.00451702744562279;0.00454725179510159;0.00457788069153352;0.00460887391937061;0.00464020552540860;0.00467185701975483;0.00470381282284363;0.00473605808015047;0.00476857891718795;0.00480136203783277;0.00483439436476723;0.00486766296061053;0.00490114045861107;0.00493483136918280;0.00496872391532794;0.00500279299551623;0.00503705227465459;0.00507145385224523;0.00510601754975535;0.00514071502552731;0.00517553317452408;0.00521048453054189];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0613467709794854;0.125216002053204;0.189797148258713;0.254333735443525;0.318487743791413;0.382114062932921;0.445145758933416;0.507557333109014;0.569337519637458;0.630478088807434;0.690975090007418;0.750827452758148;0.810034490460835;0.868595139726926;0.926506845120101;0.983769962822696;1.04038353249203;1.09634539445447;1.15165672459813;1.20631343122552;1.26031720750166;1.31366562477631;1.36635756005405;1.41839408737959];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000339557118480829;0.000593618523869692;0.000809621916849858;0.00100779586570128;0.00119725794083937;0.00138196729145419;0.00156373182144740;0.00174326152162490;0.00192085104237566;0.00209669711603550;0.00227086937494662;0.00244337749460406;0.00261422218516008;0.00278340669390784;0.00295093301885284;0.00311680094430453;0.00328100871194106;0.00344355378182944;0.00360443749585012;0.00376364876754363;0.00392119218454272;0.00407705999500611;0.00423124861452216;0.00438376147126013];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;1.53029752325150e-06;2.99264583798051e-06;4.41068275227863e-06;5.79459867040075e-06;7.14935514801869e-06;8.47750036062795e-06;9.78060923714028e-06;1.10597045031806e-05;1.23156068646828e-05;1.35490849447435e-05;1.47608267204845e-05;1.59514455887548e-05;1.71215220664733e-05;1.82716165258114e-05;1.94028686952344e-05;2.05150858790120e-05;2.16087210822438e-05;2.26847531079636e-05;2.37426485366423e-05;2.47841875401369e-05;2.58086334743688e-05;2.68169802646106e-05;2.78096874968083e-05;2.87863156180955e-05];
 
 State_diff_8=State8-State_data_8;
 
 %%%%%
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sqrt(sum from 1 to N of (diff#)^2)
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; gives least squares percentage error so each piece weighted evenly)
 
 % For testing purposes with states 
 %value = norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8);
 
 % For testing purposes with states and data sets
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8);
 
 % Objective function value we wish to minimize; want value=fval(x) to be small  when run MultiStart
 
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4);
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 %value=norm(Diff4,2)./norm(Data4);
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff4,2)./norm(Data4);
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff4,2)./norm(Data4);


 
end
 

function f = HeroinModel(t,y,pars)
f=zeros(8,1);
f(1)=-(pars(1)*t+pars(16))*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=(pars(1)*t+pars(16))*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time;
% i.e. individuals who enter the P class at any time from S (used in Estim1)
f(6) =(pars(1)*t+pars(16))*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
% i.e. individuals who enter the A class at any time (used in Estim2)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction over time; 
% i.e. individuals who enter the H class at any time (used in Estim3)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));


end

