%File name: heroin_multistart_final.m 

clf;
clear all;

%Make sure check what data points included! 

% We wish to estimate the parameter vector (10 parameters)
% x =[m,theta_1,epsilon,gamma,sigma,b,P0,A0,H0,R0]
% Ranges on each of the parameters 

%Taking out estimating theta_1 and seeing fit 
%LowerBounds=[-0.1   0.000001     0.8      0.001        0.001    0.1   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.0001      8         0.01         0.1      1       0.5     0.5     0.5     ];


%THE ONE from 4/25/19 when not using 2016 heroin data point 
%LowerBounds=[-0.1   0.000001     0.8      0.001        0.001    0.1   0.0001  0.0001  0.00001  ];
%UpperBounds=[0.1     0.0001      8         0.1         0.1      1       0.5     0.5     0.5     ];


LowerBounds=[-0.1   0.000001     0.8      0.0001        0.0001    0.1   0.0001  0.0001  0.00001  ];
UpperBounds=[0.1     0.0001       8         0.1          0.1        1     0.5     0.5     0.5     ];

%OKAY!!!-taking out 2016 value for heroin fval=.1192 sigma, gamma, theta
%LowerBounds=[-0.1   0.00001     0.8      0.001       0.0001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.0001      8        0.01         0.1      1       0.5     0.5     0.5     ];

%%THIS ONE- 4/24/18
%LowerBounds=[-0.1   0.000001     0.8      0.001       0.0001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.0001      8           0.01         0.1      1       0.5     0.5     0.5     ];

%OR THIS ONE WITH HIGHER SIGMA LOWER BOUND (BEST)
%LowerBounds=[-0.1   0.000001     0.8      0.001       0.001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.0001      8           0.01         0.1      1       0.5     0.5     0.5     ];

%LowerBounds=[-0.1   0.000001     0.8      0.001       0.0001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.00001      8        0.01         0.1      1       0.5     0.5     0.5     ];


%Take out theta_1
%LowerBounds=[-0.1       0.8     0.001      0.00001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1         8      0.01         0.1      1       0.5     0.5     0.5     ];

%Runs all converge fval=.2775
%LowerBounds=[-0.1   0.000001     0.8      0.001       0.00001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.00001      8        0.01         0.1      1       0.5     0.5     0.5     ];


%LowerBounds=[-0.1   0.00001   0.8     0.001       0.00001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.1       8       0.5          0.1      1       0.5     0.5     0.5     ];


%Result going to go with from 4/10; bounds have reasons behind them and low enough objective function value 
%LowerBounds=[-0.1   0.00001   0.8    0.0001       0.00001   0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.1       8       0.5          0.1      1       0.5     0.5     0.5     ];

%Result COULD go with from 4/21; convergence happens on more runs 
%LowerBounds=[-0.1   0.0001   0.8     0.0001       0.00001     0.01   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.1       8      0.099        0.099      1       0.5     0.5     0.5     ];


%Good result on 4/8/19 since don't want to estimate P0 and use calculated
%value instead (gives about same objective value as letting multistart estimate P0) 
%LowerBounds=[-0.1   0.00001   0.2    0.0001    0.00000001   0.000001   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.1       6       0.2          0.1          1       0.5     0.5  0.5     ];

%REALLY GOOD 4/10 (moved epsilon to range in opioid paper)
%LowerBounds=[-0.1   0.00001   0.8    0.0001    0.00000001   0.000001   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.1       8       0.2          0.1          1       0.5     0.5  0.5     ];

%ALSO GOOD 4/10
%LowerBounds=[-0.1   0.00001   0.8    0.0001       0.00001   0.000001   0.0001  0.0001  0.0001  ];
%UpperBounds=[0.1     0.1       8       0.2          0.1          1       0.5     0.5  0.5     ];


%Best we can do on 4/8/19 for realistic numbers although theta_1 still hits upper bound; includes estimating P0
%LowerBounds=[-0.1   0.00001   0.2   0.0001    0.000001     0.0001  0.0001  0.0001 0.0001  0.0001 ];
%UpperBounds=[0.1     0.1       6       0.2        0.1          1       0.5     0.5     0.5    0.5   ];
%x=[-0.0120382648602391,0.0997531214323014,2.10018581425436,0.000112314487481875,0.00268044692122784,0.269298366811826,0.0963523859585739,0.00754523097086079,0.00118983997840967,0.000641822830099002];
%fval=0.1819

% Best result for objective function (on 4/4/19) BUT R0 just way too high so need to adjust
%LowerBounds=[-0.1   0.0001    0.2     0.0001    0.00000001    0.000001  0.0001  0.0001  0.0001  0.0001 ];
%UpperBounds=[0.1      1       6       0.2         0.1          1        0.5     0.5     0.5     0.5   ];

%For testing purposes when change code
%LowerBounds=[-0.1  0.000001    0.1     0.00001    0.01    0.001   0.000001 0.000001 0.000001 0.00001  ];
%UpperBounds=[0.1      0.2        6      0.1        1        1       0.5      0.5      0.5     0.5  ];

%Gives best result with alpha time dependent when not estimating P0 and A0,
%only H0 and R0
%LowerBounds=[-0.1  0.001    0.1     0.00001    3     0.001     0.000001 0.00001];
%UpperBounds=[0.1      2        6      0.01       15          1        0.1     0.2  ];

%Get  very few runs that converge when make bounds much much wider
%LowerBounds=[-1  0.00001    0.00001   0.0000000001    3     0.000001     0.00001 0.00001];
%UpperBounds=[1      2        6            1           9        1        0.5     0.5  ];

%Run 4 from my paper log
%LowerBounds=[-0.1  0.001    0.1     0.0000001    0.0001    0.001  0.000001   0.0001];
%UpperBounds=[0.1      2        6      0.01          4        1          0.1     4];

%Run 7 from my paper log 
%LowerBounds=[-0.1  0.001    0.1     0.0000001    0.001    0.001   0.00001  0.000001   1];
%UpperBounds=[0.1      2        6      0.01          3        1      0.2      0.1     10];

%Run 14 from my paper log
%LowerBounds=[-0.1  0.001    0.1     0.0000001    0.001    0.001     0.000001   1];
%UpperBounds=[0.1      2        6      0.01          3        1        0.1     10];

%Run 16 from my paper log 
%LowerBounds=[-0.1  0.0000001   0.1    0.0000001    0.0000001    0.001   0.00001  0.0000001   ];
%UpperBounds=[0.1      1        5         0.01          1          1       0.2       5];

%Run 17 from my paper log 
%LowerBounds=[-0.1  0.0000001   0.1    0.1    0.0000001    0.001   0.00001  0.0000001   ];
%UpperBounds=[0.1      1        5        2       0.0001          1       0.2    0.0001];

%Run 18 from my paper log 
%LowerBounds=[-0.1  0.0000001   0.1    0.00001    0.0000001    0.001   0.00001  0.0000001   ];
%UpperBounds=[0.1      1        10        2       0.0001          1       0.3    0.1];

%Come back from my paper log 
%LowerBounds=[-0.1  0.0000001   0.001    0.00001    0.0000001    0.001   0.00001  0.0000001   ];
%UpperBounds=[0.1      1        10        2       0.0001          1       0.3    5];

%Without estimating R0 from my paper log 
%LowerBounds=[-1  0.000001   0.001    0.00001    0.00000001     0.000001   0.00000001];
%UpperBounds=[1     0.1        4       0.2          0.1            1          5];

%Run 19 great when splitting up sigma into sigma_A and sigma_H, from my paper log 
%LowerBounds=[-1  0.000001   0.001    0.0001    0.00000001     0.000001  0.00001  0.00001 0.00001 0.00001 0.00001 ];
%UpperBounds=[1     0.1        4       0.2          0.1            1        5       0.5      0.5     0.5     0.5   ];
%RESULTS of Run 19: 
%x= -0.0368    0.0025    2.1339    0.0001    0.0000    0.6012    0.0042    0.0878    0.0076    0.0009    0.4780
%fval=0.1417

%Run 20 from my paper log 
%LowerBounds=[-1  0.000001   0.001    0.000001    0.00000001     0.000001  0.00001  0.00001 0.00001 0.00001 0.00001 ];
%UpperBounds=[1     0.1        4       0.2          0.1            1        5       0.5  0.5   0.5   0.5   ];
%RESULTS of Run 20:
%x = -0.0188    0.0549    3.9640    0.0000    0.0000    0.3985    0.0071    0.0530    0.0077    0.0009    0.2151
%fval=0.1392

%Run 21 from my paper log 
%LowerBounds=[-1  0.000001   0.001    0.000001    0.00000001     0.000001  0.00001  0.00001 0.00001 0.00001 0.00001 ];
%UpperBounds=[1     0.1        6       0.2          0.1            1        5         0.5     0.5     0.5     0.5   ];
%RESULTS of Run 21:
%x = -0.0119    0.0969    2.1429    0.0000    0.0001    0.2729    0.0870    0.0964    0.0076    0.0008    0.0124
%fval=0.1437

%Run 23 few converge though because of -1 to 1 for m, from my paper log 
%LowerBounds=[-1   0.000001   0.001    0.0001    0.00000001     0.000001  0.00001 0.00001 0.00001 0.00001 ];
%UpperBounds=[1      1        6       0.2          0.1            1          0.5     0.5     0.5     0.5   ];

%Run 23 with Lenhart in office, good enough, from my paper log 
%LowerBounds=[-0.1   0.0001     1    0.000001    0.00000001     0.000001 0.0001 0.0001 0.0001 0.0001 ];
%UpperBounds=[0.1      1        6       0.2          0.1          1       0.5     0.5     0.5    0.5   ];
%x= -0.0210    0.2562    3.5215    0.0001    0.0000    0.4609    0.0648    0.0077    0.0009    0.3288
%fval=0.1350

%Run 23b made lower bound on epsilon lower for bigger range, one with Lenhart in office, from my paper log 
%LowerBounds=[-0.1   0.0001     0.2    0.000001    0.00000001     0.000001   0.0001   0.0001   0.0001  0.0001 ];
%UpperBounds=[0.1      1        6       0.2          0.1              1       0.5       0.5      0.5     0.5   ];
%x=[-0.0232938669941276,0.274491197674059,3.80589997869509,5.50775141421631e-06,9.35412160715760e-07,0.507882966712332,0.0609677343429615,0.00775635360045256,0.000856017464737244,0.379999731834762]
%fval=0.1326

%Run 24 from my paper log 
%LowerBounds=[-0.1   0.0001     1    0.000001    0.00000001     0.000001 0.00001 0.0001 0.0001 0.0001 ];
%UpperBounds=[0.1      1        6       0.2          0.1          1       0.5     0.5     0.5    0.5   ];
%RESULT from Run 24
% [-0.0196064313678204,0.257382363769185,4.55226111301852,6.21116506528310e-06,1.15975751417718e-06,0.451198543349247,0.0510786877893240,0.00774917078493417,0.000853637843509391,0.299319438962484]
% fval= 0.1331


%Run 26b increasing LB on gamma again, from my paper log 
%LowerBounds=[-0.1   0.0001    0.2     0.001    0.00000001    0.000001  0.0001  0.0001  0.0001  0.0001 ];
%UpperBounds=[0.1      1       6       0.2         0.1          1        0.5     0.5     0.5     0.5   ];

%Run 25c, increasing UB on theta_1, from my paper log 
%LowerBounds=[-0.1   0.0001     1     0.0001    0.00000001    0.000001  0.0001  0.0001  0.0001  0.0001 ];
%UpperBounds=[0.1     0.5       6       0.2         0.1          1        0.5     0.5     0.5     0.5   ];
%x=[-0.0134957161092346,0.204489901975793,5.78352777097860,0.000130040100191768,3.90947534199090e-05,0.331340127979528,0.0385444335706383,0.00769974455765799,0.000857642856768981,0.0422983312355828]
%fval=0.1359


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
numstartpoints=30;

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
theta_2=3*x(2);
sigma=x(5);
zeta=0.0214;
theta_3=16*x(2);
nu=0.0155;
omega=0.0000000001;
b=x(6);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Print optimal parameter solution and objective function value in command
% window when completed 
format short 
x
fval

% Final time N; will run 2013-2018 where t=0 represents 2013
% and t=5 represents 2018, with spacing (T-0)/((N+1)-1)=1 between the points:
N = 5; 
tspan=linspace(0,N,N+1);


% Initial conditions
P0=0.0710;
A0=x(7);
H0=x(8);%0.0001;%0.00136;%x(7);
R0=x(9);
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
  for i=1:N+1
      sum(i)=y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5);
  end
  
 %%% For testing purposes: states and corresponding simulated data 
 State1=y(:,1);
 State_data_1=[0.915220000000000;0.909363647939647;0.912195235212623;0.915271066237107;0.918375056138264;0.921495555512809];
     
 
 % Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1))
 plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 legend('Proportion of susceptibles')%,'Proportion of susceptibles (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
 
 
 State2=y(:,2);
 State_data_2=[0.0710000000000000;0.0769912098590027;0.0742870075803353;0.0713354446896955;0.0683531714139834;0.0653522800264156];
 
 
 % Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2))
 plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 legend('Proportion of prescription users')%,'Proportion of prescription users (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})

 State3=y(:,3);
 State_data_3=[0.00737000000000000;0.00722303269879493;0.00708084842482764;0.00693951595311553;0.00679885789635725;0.00665879294390987];
 
 
 % Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3))
 plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend('Proportion of opioid addicts')%,'Proportion of opioid addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
 
 State4=y(:,4);
 State_data_4=[0.00173000000000000;0.00160625207987268;0.00149140244485008;0.00138478933359258;0.00128583422886192;0.00119399639447003];
 
 % Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4))
 plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend('Proportion of heroin/fentanyl addicts')%,'Proportion of heroin/fentanyl addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})

 State5=y(:,5);
 State_data_5=[0.00468000000000000;0.00481585742213024;0.00494550633579949;0.00506918378463544;0.00518708032080239;0.00529937512066341];
 
 % Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5))
 plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 legend('Proportion of stably recovered addicts')%,'Proportion of stably recovered addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})

 
 State6=y(:,6);
 State_data_6=[0;0.262183673545998;0.513812866430539;0.755489669560203;0.987164212037442;1.20877102276257];
 
 % Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6))
 plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend('Proportion that enter P at some point during the year')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})

 State7=y(:,7);
 State_data_7=[0;0.000135965936828635;0.000271146594528454;0.000401682477506223;0.000527424046729980;0.000648314109552405];
 
 % Simulated data for L and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7))
 plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend('Proportion that enter A at some point during the year')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
 
 
 State8=y(:,8);
 State_data_8=[0;1.10245742605816e-06;2.17546374550268e-06;3.21771718984808e-06;4.22894458150900e-06;5.20904900937005e-06];
 
 % Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8))
 plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend('Proportion that enter H at some point during the year')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018'})
 %%%
 
 %%% Data points we are interested in 
 
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 %Data1=[0.399384466780766;0.476721593432771;0.469765087997695;0.449866817308604;0.427737148099884];
 
 % Actual Data for years 2013-2017
 Data1=[1825910./5517176; 1805325./5559006; 1800613./5602117; 1744766./5651993; 1620955./5708586];
 %Data1=[0.333183673545998;0.328620402743543;0.315963810709999;0.303009987166935;0.289959982139115];
 
 % Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 plot(t(1:end-1),Estim1, 'o')
 plot(t(1:end-1), Data1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
 
 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 %Data2=[0.00709261474856600;0.0106675766257930;0.0130723402928730;0.0148284456654410;0.0162165492691030];
 
 % Actual Data for years 2013-2017 
 Data2=[43418./5517176; 42928./5559006; 42816./5602117; 37464./5651993; 34805./5708586];
 %Data2=[0.00750596593682863;0.00735821335649475;0.00721138430780541;0.00706525752233929;0.00691974795917967];
 
 
 % Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 plot(t(1:end-1),Estim2, 'o')
 plot(t(1:end-1), Data2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [0 1 2 3 4])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017'})
 

 %Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8); 
 Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8); 
 %Data3=[0.00116527288223448;0.00120952017524577;0.00118883157707289];
 
 % Actual Data for years 2014-2016
 Data3=[7560./5559006; 7560./5602117; 10260./5651993];
 %Data3=[0.00160732508619213;0.00149244469829442;0.00138580056098424];
 
 %Data3=[7560./5559006; 7560./5602117];
 % Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(11)
 hold all
 plot(t(2:4),Estim3, 'o')
 plot(t(2:4), Data3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 1 2 3 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 
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
theta_2=3*z(2); 
sigma=z(5);
zeta=0.0214;
theta_3=16*z(2);
nu=0.0155;
omega=0.0000000001;
b=z(6);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Final time N; will run 2013-2018 where t=0 represents 2013
% and t=5 represents 2018, with spacing (T-0)/((N+1)-1)=1 between the points
N = 5; 
tspan=linspace(0,N,N+1);

% Initial conditions
P0=0.0710;
A0=z(7);
H0=z(8);%0.0001;%0.00136;%x(7);
R0=z(9);
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

 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  

 % Actual proportions of population (updated 3/12/19) that were non-addicted prescription opioid users at some point
 % during the year for 2013-2017 
 % (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 
 Data1=[1825910./5517176; 1805325./5559006; 1800613./5602117; 1744766./5651993; 1620955./5708586];
  %Data1=[0.333183673545998;0.328620402743543;0.315963810709999;0.303009987166935;0.289959982139115];
  
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
 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 % When testing all points with simulated data
 %Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 % Actual proportions of population (updated 3/12/19) that were opioid addicted individuals in
 % the population at some point during the year in 2014 and 2015 
 % (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 Data2=[43418./5517176; 42928./5559006; 42816./5602117; 37464./5651993; 34805./5708586];
 %Data2=[0.00750596593682863;0.00735821335649475;0.00721138430780541;0.00706525752233929;0.00691974795917967];
 
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
 Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8);  
 % When testing all points with simulated data 
 %Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);
 
 %Data3=[0.00160732508619213;0.00149244469829442;0.00138580056098424];
 
 % Actual proportion (updated 3/12/19) of heroin addicted individuals in the population at some point during the year 
 % in 2014, 2015, and 2016
 % (total number of heroin addicted individuals in 2014, 2015, and 2016 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 Data3=[7560./5559006; 7560./5602117; 10260./5651993];
 %Data3=[7560./5559006; 7560./5602117];
  
 % Data simulated when testing codes 
 %Data3=[0.00116527288223448;0.00120952017524577;0.00118883157707289];
 
 % The difference between estimated value and data
 
 Diff3=Estim3-Data3;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.915220000000000;0.909363647939647;0.912195235212623;0.915271066237107;0.918375056138264;0.921495555512809];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.0710000000000000;0.0769912098590027;0.0742870075803353;0.0713354446896955;0.0683531714139834;0.0653522800264156];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00737000000000000;0.00722303269879493;0.00708084842482764;0.00693951595311553;0.00679885789635725;0.00665879294390987];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.00173000000000000;0.00160625207987268;0.00149140244485008;0.00138478933359258;0.00128583422886192;0.00119399639447003];
 
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.00468000000000000;0.00481585742213024;0.00494550633579949;0.00506918378463544;0.00518708032080239;0.00529937512066341];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.262183673545998;0.513812866430539;0.755489669560203;0.987164212037442;1.20877102276257];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000135965936828635;0.000271146594528454;0.000401682477506223;0.000527424046729980;0.000648314109552405];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;1.10245742605816e-06;2.17546374550268e-06;3.21771718984808e-06;4.22894458150900e-06;5.20904900937005e-06];
 
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
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8);
 % Objective function value we wish to minimize; want value=fval(x) to be small  when run MultiStart
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);



 
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

