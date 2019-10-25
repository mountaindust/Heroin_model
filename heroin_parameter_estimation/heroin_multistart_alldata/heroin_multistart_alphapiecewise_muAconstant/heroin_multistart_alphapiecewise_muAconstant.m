%File name: heroin_multistart_quarters.m 

clf;
clear all;

% Realistic parameter bounds
%           [m      betaA     betaP   theta1   epsilon  gamma   theta2   sigma    zeta   theta3       nu     b    P0        A0       H0       R0       c  ]
LowerBounds=[-0.1  0.00001   0.00001   0.05      1       0.005    0.1     0.1    0.0001    10      0.0001   0.1   0.001   0.0001   0.00001  0.00001  -0.1  ];
UpperBounds=[-0.001  0.01     0.01     0.3       5        0.1     0.4       2    0.2       20       0.2     0.5   0.35      0.01    0.002     0.1    -0.001 ];
%LowerBounds=[-0.1  0.00001   0.00001   0.0001    1      0.005   0.0001     0.1    0.0001   0.01   0.0001   0.1   0.001   0.0001   0.00001  0.00001   -0.1  ];
%UpperBounds=[-0.001  0.01     0.01     0.01      5       0.1     0.05      2       0.2      0.5    0.2     0.5   0.35      0.01    0.002     0.1    -0.001 ];
%LowerBounds=[-0.1  0.00001   0.00001   0.001      1      0.005   0.02     0.1    0.0001     0.5    0.0001   0.1   0.001   0.0001   0.00001  0.00001   -0.1  ];
%UpperBounds=[-0.001  0.01     0.01     0.075      5       0.1     0.5      2       0.2     0.93    0.2     0.5   0.35      0.01    0.002     0.1    -0.001  ];
%LowerBounds=[-0.1  0.00001   0.00001   0.03      1      0.005    0.05     0.1    0.0001   0.5    0.0001   0.1   0.001   0.0001   0.00001  0.00001  -0.1 ];
%UpperBounds=[-0.001  0.01     0.01     0.15      5       0.1       1       2       0.2    1.85    0.2     0.5   0.35      0.01    0.002     0.1   -0.001];



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
c=x(17);
%d=x(18);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c];%,d];



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
P0=x(13);
A0=x(14);
H0=x(15);%0.0001;%0.00136;%x(7);
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
  %alpha=-pars(1)*t+pars(16);
  %alpha=x(1)*t+x(12);
  
 % Making sure S+P+A+H+R=1
 total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);
 
 %%% For testing purposes: states and corresponding simulated data 
 %%% For testing purposes: states and corresponding simulated data 
 State1=y(:,1);
 State_data_1=[0.911750000000000;0.905680727955066;0.902813042597302;0.901523613264462;0.900935031988554;0.900694391797014;0.900638160052312;0.900664213839019;0.900722410519065;0.900795143540474;0.900876217661415;0.900960863386460;0.901107231319879;0.901374182510484;0.901795540151398;0.902520354268059;0.903586175135996;0.904815887911269;0.906212853297238;0.907623832901443;0.909086827173730;0.910520668888514;0.911951178310007;0.913365165701522;0.914767207914087];
 
 
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
 State_data_2=[0.0800000000000000;0.0859490173625646;0.0886881050936439;0.0898438931926504;0.0902954905320082;0.0903963737672140;0.0903102233191685;0.0901391636572149;0.0899332071697220;0.0897097408417024;0.0894746922616239;0.0892325401477153;0.0889248988442954;0.0884926833271118;0.0879018628929939;0.0870035274028070;0.0857598900161739;0.0843480090254663;0.0827645199164663;0.0811618907325020;0.0795022320784812;0.0778654260050542;0.0762258403775289;0.0745952029674074;0.0729693034481355];
 
 
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
 State_data_3=[0.00760000000000000;0.00752730594326424;0.00746565175889461;0.00741106670530141;0.00736117071937759;0.00731472925780572;0.00727093143058155;0.00722910399462378;0.00718869049329482;0.00714920734521127;0.00711019732469982;0.00707120688884708;0.00703168281456537;0.00699106593675066;0.00694878857417397;0.00690384020241415;0.00685568392484352;0.00680315094228621;0.00674622496486570;0.00668357463384753;0.00661548199069936;0.00654104124259353;0.00646044961163867;0.00637302865690354;0.00627896276577220];
 
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
 State_data_4=[0.000600000000000000;0.000646196689087267;0.000696271162975088;0.000750522163407866;0.000809282140589654;0.000872921957772847;0.000941847647557150;0.00101649193832210;0.00109732316199227;0.00118484563706514;0.00127960345007685;0.00138217274896256;0.00149320565660478;0.00161334562894590;0.00174322627429242;0.00188366375459332;0.00203528595122269;0.00219927529055221;0.00237564815551886;0.00256641228003148;0.00277119357512908;0.00299215071686700;0.00322890834036054;0.00348359885316700;0.00375567395622300];
 
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
 State_data_5=[5.00000000000000e-05;0.000196752049973787;0.000336929386979425;0.000470904673961105;0.000599024619420713;0.000721583220285283;0.000838837550536363;0.000951026570998501;0.00105836865617952;0.00116106263590325;0.00125928930263327;0.00135321682852846;0.00144298136521223;0.00152872259729155;0.00161058210777824;0.00168861437312572;0.00176296497344789;0.00183367683254998;0.00190075366823626;0.00196428945422765;0.00202426518336119;0.00208071314777032;0.00213362336060356;0.00218300382060177;0.00222885191498770];
 
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
 State_data_6=[0;0.0605520659020888;0.120647236512263;0.180435233392232;0.239997109144841;0.299365189001551;0.358554288931915;0.417573609637546;0.476427627160312;0.535117701845272;0.593644111627245;0.652007012135567;0.710098584449818;0.767807821071126;0.825043709217991;0.881396018556176;0.936784639099194;0.991110719956659;1.04436879639724;1.09652311996434;1.14758256761387;1.19757508752387;1.24649343003076;1.29434688262690;1.34113249487163];
 
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
 State_data_7=[0;0.000129895380459722;0.000270664499136368;0.000418675192241953;0.000571832365078841;0.000729142984153997;0.000890019510841202;0.00105400682698996;0.00122076275148860;0.00139001667674692;0.00156152503011263;0.00173504917951496;0.00191025030869162;0.00208678327432817;0.00226429100573467;0.00244195449342770;0.00261945009372684;0.00279582272470218;0.00297106183526447;0.00314421760862893;0.00331549468802574;0.00348424248649901;0.00365060099728037;0.00381407205541882;0.00397477834811951];
 
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
 State_data_8=[0;5.46395000795889e-05;0.000113806105405925;0.000177854799328613;0.000247177282136865;0.000322208507167082;0.000403424056474009;0.000491330699730511;0.000586477020227823;0.000689454507286418;0.000800902185199254;0.000921497548563281;0.00105201101334084;0.00119320266132300;0.00134582140607685;0.00151084495484937;0.00168901448253222;0.00188176544946250;0.00208911711932428;0.00231350207207313;0.00255446790990392;0.00281465434495605;0.00309360775873427;0.00339396808851574;0.00371505972804691];
 
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
 State_data_9=[0;1.67131617260577e-05;3.32792520365464e-05;4.97167808008344e-05;6.60392251906314e-05;8.22554892795048e-05;9.83719812248542e-05;0.000114393905710695;0.000130325114180618;0.000146168195401521;0.000161924584475748;0.000177594891757280;0.000193177428704234;0.000208670729419151;0.000224073345975198;0.000239376857807449;0.000254580200988510;0.000269667575703308;0.000284638959485607;0.000299469095708998;0.000314162627316002;0.000328688620269372;0.000343052161331274;0.000357221050493541;0.000371203390093400];
 
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
 State_data_10=[0;7.25894819613566e-06;1.50761682145009e-05;2.35000153847639e-05;3.25814316622060e-05;4.23755203283973e-05;5.29420785160385e-05;6.43447653446587e-05;7.66525792221662e-05;8.99404677564156e-05;0.000104290056398674;0.000119788480839637;0.000136537446937669;0.000154637227101435;0.000174186958326299;0.000195325533802829;0.000218150360872733;0.000242878992722493;0.000269513906156840;0.000298420374713506;0.000329531070124557;0.000363259576921112;0.000399538614323286;0.000438804545522494;0.000480947914521110];
 
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
 %Data1=[0.319997109144841;0.326726008547479;0.323604164459228;0.315610953493671;0.296557818530854;0.273052159336238];
 
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
 %Data2=[0.00817183236507884;0.00801010110578735;0.00787817805049785;0.00774088259960058;0.00755172851914243;0.00727476565079313];
 
 
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
 %Data3=[0.00114858187868061;0.00156285715510528;0.00213020912579616];
 
 
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
 %Data4=[0.140552065902089;0.146044187972739;0.148476101973613;0.149405768945259;0.149663570388719;0.149585473697578;0.149329544024800;0.148993181179980;0.148623281854683;0.148236150623675;0.147837592769946;0.147324112461966;0.146634135465604;0.145728571473976;0.144254172231180;0.142392147945824;0.140085970873639;0.137606085466047;0.134918843483572;0.132221338382031;0.129494751988474;0.126783768511947;0.124079292973674;0.121380815212135];
 
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
 %Data5=[6.60392251906314e-05;6.42858889899863e-05;6.28523145236163e-05;6.14027722842754e-05];
 
 
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
 %Data6=[3.25814316622060e-05;4.40711475599602e-05;5.98848677155027e-05;8.16129139350641e-05;0.000111380709251824];
 
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
c=z(17);
%d=z(18);

% Parameter vector
pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c];%,d];

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
 %Data1=[0.319997109144841;0.326726008547479;0.323604164459228;0.315610953493671;0.296557818530854;0.273052159336238];
 
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
 %Data2=[0.00817183236507884;0.00801010110578735;0.00787817805049785;0.00774088259960058;0.00755172851914243;0.00727476565079313];
 
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
 %Data3=[0.00114858187868061;0.00156285715510528;0.00213020912579616];
 

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
 %Data4=[0.140552065902089;0.146044187972739;0.148476101973613;0.149405768945259;0.149663570388719;0.149585473697578;0.149329544024800;0.148993181179980;0.148623281854683;0.148236150623675;0.147837592769946;0.147324112461966;0.146634135465604;0.145728571473976;0.144254172231180;0.142392147945824;0.140085970873639;0.137606085466047;0.134918843483572;0.132221338382031;0.129494751988474;0.126783768511947;0.124079292973674;0.121380815212135];
 

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
 %Data5=[6.60392251906314e-05;6.42858889899863e-05;6.28523145236163e-05;6.14027722842754e-05];
 
 
 %The difference between estimated value and data
 Diff5=Estim5-Data5;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from H class at some point during the year for
 % years 2013-2017, Estim6 is a column vector
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10); y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual proportion of heroin addicts that overdose each year 2013-2017
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 
 %Testing Data
 %Data6=[3.25814316622060e-05;4.40711475599602e-05;5.98848677155027e-05;8.16129139350641e-05;0.000111380709251824];
 
 

 %The difference between estimated value and data 
 Diff6=Estim6-Data6;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.911750000000000;0.905680727955066;0.902813042597302;0.901523613264462;0.900935031988554;0.900694391797014;0.900638160052312;0.900664213839019;0.900722410519065;0.900795143540474;0.900876217661415;0.900960863386460;0.901107231319879;0.901374182510484;0.901795540151398;0.902520354268059;0.903586175135996;0.904815887911269;0.906212853297238;0.907623832901443;0.909086827173730;0.910520668888514;0.911951178310007;0.913365165701522;0.914767207914087];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.0800000000000000;0.0859490173625646;0.0886881050936439;0.0898438931926504;0.0902954905320082;0.0903963737672140;0.0903102233191685;0.0901391636572149;0.0899332071697220;0.0897097408417024;0.0894746922616239;0.0892325401477153;0.0889248988442954;0.0884926833271118;0.0879018628929939;0.0870035274028070;0.0857598900161739;0.0843480090254663;0.0827645199164663;0.0811618907325020;0.0795022320784812;0.0778654260050542;0.0762258403775289;0.0745952029674074;0.0729693034481355];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00760000000000000;0.00752730594326424;0.00746565175889461;0.00741106670530141;0.00736117071937759;0.00731472925780572;0.00727093143058155;0.00722910399462378;0.00718869049329482;0.00714920734521127;0.00711019732469982;0.00707120688884708;0.00703168281456537;0.00699106593675066;0.00694878857417397;0.00690384020241415;0.00685568392484352;0.00680315094228621;0.00674622496486570;0.00668357463384753;0.00661548199069936;0.00654104124259353;0.00646044961163867;0.00637302865690354;0.00627896276577220];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.000600000000000000;0.000646196689087267;0.000696271162975088;0.000750522163407866;0.000809282140589654;0.000872921957772847;0.000941847647557150;0.00101649193832210;0.00109732316199227;0.00118484563706514;0.00127960345007685;0.00138217274896256;0.00149320565660478;0.00161334562894590;0.00174322627429242;0.00188366375459332;0.00203528595122269;0.00219927529055221;0.00237564815551886;0.00256641228003148;0.00277119357512908;0.00299215071686700;0.00322890834036054;0.00348359885316700;0.00375567395622300];
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[5.00000000000000e-05;0.000196752049973787;0.000336929386979425;0.000470904673961105;0.000599024619420713;0.000721583220285283;0.000838837550536363;0.000951026570998501;0.00105836865617952;0.00116106263590325;0.00125928930263327;0.00135321682852846;0.00144298136521223;0.00152872259729155;0.00161058210777824;0.00168861437312572;0.00176296497344789;0.00183367683254998;0.00190075366823626;0.00196428945422765;0.00202426518336119;0.00208071314777032;0.00213362336060356;0.00218300382060177;0.00222885191498770];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0605520659020888;0.120647236512263;0.180435233392232;0.239997109144841;0.299365189001551;0.358554288931915;0.417573609637546;0.476427627160312;0.535117701845272;0.593644111627245;0.652007012135567;0.710098584449818;0.767807821071126;0.825043709217991;0.881396018556176;0.936784639099194;0.991110719956659;1.04436879639724;1.09652311996434;1.14758256761387;1.19757508752387;1.24649343003076;1.29434688262690;1.34113249487163];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000129895380459722;0.000270664499136368;0.000418675192241953;0.000571832365078841;0.000729142984153997;0.000890019510841202;0.00105400682698996;0.00122076275148860;0.00139001667674692;0.00156152503011263;0.00173504917951496;0.00191025030869162;0.00208678327432817;0.00226429100573467;0.00244195449342770;0.00261945009372684;0.00279582272470218;0.00297106183526447;0.00314421760862893;0.00331549468802574;0.00348424248649901;0.00365060099728037;0.00381407205541882;0.00397477834811951];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;5.46395000795889e-05;0.000113806105405925;0.000177854799328613;0.000247177282136865;0.000322208507167082;0.000403424056474009;0.000491330699730511;0.000586477020227823;0.000689454507286418;0.000800902185199254;0.000921497548563281;0.00105201101334084;0.00119320266132300;0.00134582140607685;0.00151084495484937;0.00168901448253222;0.00188176544946250;0.00208911711932428;0.00231350207207313;0.00255446790990392;0.00281465434495605;0.00309360775873427;0.00339396808851574;0.00371505972804691];
 
 State_diff_8=State8-State_data_8;
 
 % Comparing simulated data for proportion of individuals that overdose out of A throughout the year and the model output 
 State9=y(:,9);
 State_data_9=[0;1.67131617260577e-05;3.32792520365464e-05;4.97167808008344e-05;6.60392251906314e-05;8.22554892795048e-05;9.83719812248542e-05;0.000114393905710695;0.000130325114180618;0.000146168195401521;0.000161924584475748;0.000177594891757280;0.000193177428704234;0.000208670729419151;0.000224073345975198;0.000239376857807449;0.000254580200988510;0.000269667575703308;0.000284638959485607;0.000299469095708998;0.000314162627316002;0.000328688620269372;0.000343052161331274;0.000357221050493541;0.000371203390093400];
 
 State_diff_9=State9-State_data_9;
 
 % Comparing simulated data for proportion of individuals that overdose out of H the year and the model output 
 State10=y(:,10);
 State_data_10=[0;7.25894819613566e-06;1.50761682145009e-05;2.35000153847639e-05;3.25814316622060e-05;4.23755203283973e-05;5.29420785160385e-05;6.43447653446587e-05;7.66525792221662e-05;8.99404677564156e-05;0.000104290056398674;0.000119788480839637;0.000136537446937669;0.000154637227101435;0.000174186958326299;0.000195325533802829;0.000218150360872733;0.000242878992722493;0.000269513906156840;0.000298420374713506;0.000329531070124557;0.000363259576921112;0.000399538614323286;0.000438804545522494;0.000480947914521110];
 
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
 %TRYING WEIGHTING SCHEME
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+0.5*norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+0.1*norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 %value=0.5*norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+0.1*norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 
 
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
if  t<=3.25 
    alpha = pars(1)*t+pars(16);
else 
    alpha = pars(1)*3.25+pars(16)-pars(17)*3.25+pars(17)*t;
    %alpha = pars(17)*t+pars(18);
    %alpha = pars(1)*t+pars(16);
end
end


function f = HeroinModel(t,y,pars)
f=zeros(10,1);
f(1)=-a(t,pars)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=a(t,pars)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time;
% i.e. individuals who enter the P class at any time from S (used in Estim1, Estim4)
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

