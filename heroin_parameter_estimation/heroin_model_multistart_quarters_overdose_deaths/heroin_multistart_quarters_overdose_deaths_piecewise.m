%File name: heroin_multistart_quarters.m 

clf;
clear all;

% Realistic parameter bounds
%           [m      betaA     betaP   theta1   epsilon  gamma   theta2   sigma    zeta   theta3    nu     b     P0        A0       H0       R0     c      ]
LowerBounds=[-0.1  0.00001  0.000001  0.00001   0.8    0.001  0.0001  0.0001    0.0001  0.001   0.0001  0.1   0.0001   0.00001  0.00001  0.00001  -0.1   ];
UpperBounds=[ 0.1    0.01     0.01    0.001      8       0.1       2       1        0.5     4      0.1    0.8    0.5       0.1     0.1       0.1   0.1   ];
  


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
beta_A=x(2);
beta_P=x(3);
theta_1=x(4);
epsilon=x(5);
mu=0.00868;  
mu_A=0.00870;   
mu_H=0.0507;
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
  %alpha=x(1)*t+x(12);
  
 % Making sure S+P+A+H+R=1
 total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);
 
 %%% For testing purposes: states and corresponding simulated data 
 %%% For testing purposes: states and corresponding simulated data 
 State1=y(:,1);
 State_data_1=[0.896760000000000;0.899377751002331;0.900816047207897;0.901757581088491;0.902534299291740;0.903248119572510;0.903927682709151;0.904600206538894;0.905257559478322;0.905919953833542;0.906580191645060;0.907245523915862;0.907914330143132;0.908586956055880;0.909263370744073;0.909942988896194;0.910625760025024;0.911311708454185;0.912000612504646;0.912692586845141;0.913387641327484;0.914085772842294;0.914787022454977;0.915491422798859;0.916198949847726];
 
 
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
 State_data_2=[0.0900000000000000;0.0872297772640643;0.0856441516945581;0.0845587127481872;0.0836406688471712;0.0827878363379727;0.0819714945933662;0.0811643482258755;0.0803745523152129;0.0795818297333039;0.0787934015238919;0.0780019859663431;0.0772092104284578;0.0764147180614382;0.0756185403008752;0.0748212547431420;0.0740229108891868;0.0732234825287167;0.0724231813434837;0.0716218978360846;0.0708196076065692;0.0700163183559268;0.0692119821800804;0.0684065605448062;0.0676000818208384];
 
 
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
 State_data_3=[0.00760000000000000;0.00778116171807633;0.00795598793155036;0.00812622122260124;0.00829273593999195;0.00845582620765938;0.00861560173764675;0.00877218753335612;0.00892554390370710;0.00907580790418765;0.00922293072691581;0.00936701715144157;0.00950804376092786;0.00964608331193128;0.00978113258836163;0.00991325412183929;0.0100424555473468;0.0101687503565543;0.0102922132692219;0.0104128056347826;0.0105306324718292;0.0106456606349225;0.0107579396824775;0.0108675123622870;0.0109743472013863];
 
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
 State_data_4=[0.00121000000000000;0.00120033211418663;0.00119068154652081;0.00118105991023408;0.00117147653975267;0.00116193877131192;0.00115245321158380;0.00114302886231896;0.00113366266418249;0.00112436579125213;0.00111513426424128;0.00110597730572365;0.00109689285457656;0.00108788656978219;0.00107895819871130;0.00107011113452868;0.00106134581839050;0.00105266231229878;0.00104406321198758;0.00103554717475193;0.00102711626212135;0.00101876982344630;0.00101050810058698;0.00100233135814113;0.000994239401299281];
 
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
 State_data_5=[0.00443000000000000;0.00441097790135651;0.00439313161955009;0.00437642503057575;0.00436081938134173;0.00434627911044557;0.00433276774810818;0.00432022883937859;0.00430868163839562;0.00429804273755067;0.00428834183976756;0.00427949566053897;0.00427152281284557;0.00426435600093192;0.00425799816795863;0.00425239110428559;0.00424752772004798;0.00424339634824492;0.00423992967066933;0.00423716250925937;0.00423500233202908;0.00423347834345757;0.00423254758193795;0.00423217293597428;0.00423238172882072];
 
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
 State_data_6=[0;0.0697097381718135;0.138886288563879;0.207464449644868;0.275413967917287;0.342725600697384;0.409396295114333;0.475422109817353;0.540804313890159;0.605538973923890;0.669627486685807;0.733066687941862;0.795857232919680;0.857997028013618;0.919486107878724;0.980322718625911;1.04050656148992;1.10003726978466;1.15891249762302;1.21713344276747;1.27469678448976;1.33160349101113;1.38785195447904;1.44344062501239;1.49837059971882];
 
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
 State_data_7=[0;0.000260744110198658;0.000516938667824006;0.000770283115366340;0.00102160685326603;0.00127116469332545;0.00151902671867861;0.00176525606188035;0.00200983412779941;0.00225280578852487;0.00249415506431260;0.00273389360591882;0.00297201880955431;0.00320852691012397;0.00344341803454862;0.00367668704681515;0.00390833336701094;0.00413835632846670;0.00436675111502447;0.00459352021795610;0.00481865702955093;0.00504216363530130;0.00526403643387893;0.00548427235586919;0.00570287366112160];
 
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
 State_data_8=[0;1.28915665512514e-05;2.56208759021105e-05;3.81987487725662e-05;5.06350513046481e-05;6.29379634243770e-05;7.51151444917463e-05;8.71775334200859e-05;9.91214055831169e-05;0.000110961172973100;0.000122691703224832;0.000134325859421633;0.000145860766835309;0.000157305407370460;0.000168659384162855;0.000179929321219478;0.000191116046634964;0.000202220388803647;0.000213248928508481;0.000224198261683510;0.000235076319246200;0.000245880598080347;0.000256614256140880;0.000267280062167769;0.000277875977682533];
 
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

 
 %%% Data points we are interested in 


 % Yearly simulation of individuals in P class at all during the year for years
 % 2013-2017
 
 Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6); y(21,2)+y(25,6)-y(21,6)];
     
 % Actual Data for years 2013-2017
 Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
 
 
  % Data points from proportion that is in P at some point in the year and corresponding ODE solution points 
 figure(9)
 hold all
 z1 = linspace(0,5,6); %defines mesh where going to plot Estim1, Data1 values 
 scatter(z1, Estim1, 50,'o');
 scatter(z1, Data1, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P') % at some point during the year
 legend({'ODE solution', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 5 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 

 % Yearly simulation of individuals in A class at all during the year for years
 % 2013-2017
 Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7); y(21,3)+y(25,7)-y(21,7)];
 
 
 % Actual Data for years 2013-2017 
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];

 
  % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
 figure(10)
 hold all
 z2 = linspace(0,5,6);
 scatter(z2, Estim2, 50,'o');
 scatter(z2, Data2, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A') % at some point during the year
 legend({'ODE solution', 'Data'},'FontSize', 14) 
 set(gca, 'xtick', [ 0 1 2 3 4 5 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017', '2018'})



 % Yearly simulation of individuals in H class at all during the year for years
 % 2014-2016
 
 Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8); y(13,4)+y(17,8)-y(13,8)];
 
 % Actual Data for years 2014-2016
 Data3=[7560./5559702; 7560./5602187; 10260./5648259];

 
 
 % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
 figure(11)
 hold all
 z3 = linspace(0,2,3);
 scatter(z3, Estim3, 50,'o');
 scatter(z3, Data3, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H') % at some point during the year
 legend({'ODE solution', 'Data'},'FontSize', 14)
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
        688502./5754509; 683722./5754509; 641942./5754509; 625162./5754509];
 

 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(12)
 hold all
 z4 = linspace(0,23,24);
 scatter(z4, Estim4, 50,'o');
 scatter(z4, Data4, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in P') % at some point during the year
 legend({'ODE solution', 'Data'},'FontSize', 14)
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
 
 Estim5=[mu_A*(y(1,3)+y(5,7)-y(1,7)); mu_A*(y(5,3)+y(9,7)-y(5,7)); mu_A*(y(9,3)+y(13,7)-y(9,7));...
        mu_A*(y(13,3)+y(17,7)-y(13,7)); mu_A*(y(17,3)+y(21,7)-y(17,7))];
     
 % Actual Data for years 2013-2017
 Data5=[348./5519417; 381./5559702; 463./5602187; 551./5648259; 381./5702475];
 
 
  % Data points from proportion that is in A at some point and overdoses in the year and corresponding ODE solution points 
 figure(13)
 hold all
 z5 = linspace(0,4,5); %defines mesh where going to plot Estim5, Data5 values 
 scatter(z5, Estim5, 50,'o');
 scatter(z5, Data5, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from A') % at some point during the year
 legend({'ODE solution', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})                   
                   
 
 
 
 % Yearly simulation of individuals overdosing from H class during the year for years
% 2013-2017 (anyone in H class throughout the year times mu_H)
 
 Estim6=[mu_H*(y(1,4)+y(5,8)-y(1,8)); mu_H*y(5,4)+y(9,8)-y(5,8); mu_H*(y(9,4)+y(13,8)-y(9,8));...
        mu_H*(y(13,4)+y(17,8)-y(13,8)); mu_H*(y(17,4)+y(21,8)-y(17,8))];
     
    
 % Actual Data for years 2013-2017
 Data6=[116./5519417; 216./5559702; 374./5602187; 554./5648259; 811./5702475];
 
 
  % Data points from proportion that is in H at some point and overdoses in the year and corresponding ODE solution points 
 figure(14)
 hold all
 z6 = linspace(0,4,5); %defines mesh where going to plot Estim6, Data6 values 
 scatter(z6, Estim6, 50,'o');
 scatter(z6, Data6, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from H') % at some point during the year
 legend({'ODE solution', 'Data'},'FontSize', 14)
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
mu=0.00868;  
mu_A=0.00870;   
mu_H=0.0507;
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
        688502./5754509; 683722./5754509; 641942./5754509; 625162./5754509];

 %The difference between estimated value and data 
 Diff4=Estim4-Data4;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from A class at some point during the year for
 % years 2013-2017, Estim5 is a column vector
 Estim5=[mu_A*(y(1,3)+y(5,7)-y(1,7)); mu_A*(y(5,3)+y(9,7)-y(5,7)); mu_A*(y(9,3)+y(13,7)-y(9,7));...
        mu_A*(y(13,3)+y(17,7)-y(13,7)); mu_A*(y(17,3)+y(21,7)-y(17,7))];
     
 % Actual proportion of addicted prescription opioid users that overdose
 % each year 2013-2017
 Data5=[348./5519417; 381./5559702; 463./5602187; 551./5648259; 381./5702475];
 
 %The difference between estimated value and data
 Diff5=Estim5-Data5;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from H class at some point during the year for
 % years 2013-2017, Estim6 is a column vector
 Estim6=[mu_H*(y(1,4)+y(5,8)-y(1,8)); mu_H*y(5,4)+y(9,8)-y(5,8); mu_H*(y(9,4)+y(13,8)-y(9,8));...
        mu_H*(y(13,4)+y(17,8)-y(13,8)); mu_H*(y(17,4)+y(21,8)-y(17,8))];
     
    
 % Actual proportion of heroin addicts that overdose each year 2013-2017
 Data6=[116./5519417; 216./5559702; 374./5602187; 554./5648259; 811./5702475];
 
 %The difference between estimated value and data 
 Diff6=Estim6-Data6;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.896760000000000;0.899377751002331;0.900816047207897;0.901757581088491;0.902534299291740;0.903248119572510;0.903927682709151;0.904600206538894;0.905257559478322;0.905919953833542;0.906580191645060;0.907245523915862;0.907914330143132;0.908586956055880;0.909263370744073;0.909942988896194;0.910625760025024;0.911311708454185;0.912000612504646;0.912692586845141;0.913387641327484;0.914085772842294;0.914787022454977;0.915491422798859;0.916198949847726];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.0900000000000000;0.0872297772640643;0.0856441516945581;0.0845587127481872;0.0836406688471712;0.0827878363379727;0.0819714945933662;0.0811643482258755;0.0803745523152129;0.0795818297333039;0.0787934015238919;0.0780019859663431;0.0772092104284578;0.0764147180614382;0.0756185403008752;0.0748212547431420;0.0740229108891868;0.0732234825287167;0.0724231813434837;0.0716218978360846;0.0708196076065692;0.0700163183559268;0.0692119821800804;0.0684065605448062;0.0676000818208384];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00760000000000000;0.00778116171807633;0.00795598793155036;0.00812622122260124;0.00829273593999195;0.00845582620765938;0.00861560173764675;0.00877218753335612;0.00892554390370710;0.00907580790418765;0.00922293072691581;0.00936701715144157;0.00950804376092786;0.00964608331193128;0.00978113258836163;0.00991325412183929;0.0100424555473468;0.0101687503565543;0.0102922132692219;0.0104128056347826;0.0105306324718292;0.0106456606349225;0.0107579396824775;0.0108675123622870;0.0109743472013863];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.00121000000000000;0.00120033211418663;0.00119068154652081;0.00118105991023408;0.00117147653975267;0.00116193877131192;0.00115245321158380;0.00114302886231896;0.00113366266418249;0.00112436579125213;0.00111513426424128;0.00110597730572365;0.00109689285457656;0.00108788656978219;0.00107895819871130;0.00107011113452868;0.00106134581839050;0.00105266231229878;0.00104406321198758;0.00103554717475193;0.00102711626212135;0.00101876982344630;0.00101050810058698;0.00100233135814113;0.000994239401299281];
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.00443000000000000;0.00441097790135651;0.00439313161955009;0.00437642503057575;0.00436081938134173;0.00434627911044557;0.00433276774810818;0.00432022883937859;0.00430868163839562;0.00429804273755067;0.00428834183976756;0.00427949566053897;0.00427152281284557;0.00426435600093192;0.00425799816795863;0.00425239110428559;0.00424752772004798;0.00424339634824492;0.00423992967066933;0.00423716250925937;0.00423500233202908;0.00423347834345757;0.00423254758193795;0.00423217293597428;0.00423238172882072];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0697085366826356;0.138881522876914;0.207453895535517;0.275395483031137;0.342697127061280;0.409355822862658;0.475367737349098;0.540734087140786;0.605451099008955;0.669520116667068;0.732938080676862;0.795705674782135;0.857820866766872;0.919283746769250;0.980092601258398;1.04024722032501;1.09974719177297;1.15859036230949;1.21677784577630;1.27430644410246;1.33117716664699;1.38738839363353;1.44293873635485;1.49782921724347];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000260744110198658;0.000516938667824006;0.000770283115366340;0.00102160685326603;0.00127116469332545;0.00151902671867861;0.00176525606188035;0.00200983412779941;0.00225280578852487;0.00249415506431260;0.00273389360591882;0.00297201880955431;0.00320852691012397;0.00344341803454862;0.00367668704681515;0.00390833336701094;0.00413835632846670;0.00436675111502447;0.00459352021795610;0.00481865702955093;0.00504216363530130;0.00526403643387893;0.00548427235586919;0.00570287366112160];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;1.28915665512514e-05;2.56208759021105e-05;3.81987487725662e-05;5.06350513046481e-05;6.29379634243770e-05;7.51151444917463e-05;8.71775334200859e-05;9.91214055831169e-05;0.000110961172973100;0.000122691703224832;0.000134325859421633;0.000145860766835309;0.000157305407370460;0.000168659384162855;0.000179929321219478;0.000191116046634964;0.000202220388803647;0.000213248928508481;0.000224198261683510;0.000235076319246200;0.000245880598080347;0.000256614256140880;0.000267280062167769;0.000277875977682533];
 
 State_diff_8=State8-State_data_8;
 
 %%%%%
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors [norm gives sqrt(sum from 1 to N
 % of (Diff#)^2)]
 % normalized by norm of the data [because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; gives least squares percentage
 % error so each piece weighted evenly].
 
 
 % Objective function value we wish to minimize; want value=fval(x) to be small  when run MultiStart
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4);

 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);

 % For testing purposes with states and data sets
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
f=zeros(8,1);
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


end

