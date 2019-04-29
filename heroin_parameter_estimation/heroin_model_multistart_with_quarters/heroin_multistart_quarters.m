%File name: heroin_multistart_quarters.m 

clf;
clear all;

%LowerBounds=[-0.1      0.001     0.001    0.001     1     0.0001    0.0001      0.0001     0.01      0.001       0.01        0.01   0.01     0.01    0.1    0.001   ];
%UpperBounds=[0.1       0.01      0.01     0.01      3      0.001     0.001       0.001       0.1       0.01       0.1         0.1    0.1       0.1    0.3    0.01    ];
 

%            [m       beta_A     beta_P    theta_1  epsilon   mu       mu_A        mu_H      gamma    theta_2     sigma       zeta  theta_3      nu       b     R0    
LowerBounds=[-0.1     0.00001    0.00001   0.00001     1     0.0001    0.0001      0.001      0.0001    0.001     0.01        0.01   0.001      0.001    0.01    0.001 ];
UpperBounds=[0.1       0.0001      0.0001   0.0001      3      0.001     0.001       0.01      0.001      0.01     0.1         0.09   0.01       0.1       0.2    0.01 ];
 
%Wider bounds
%LowerBounds=[-0.1     0.00001    0.00001   0.00001     1     0.0001    0.0001      0.001      0.0001    0.001     0.001        0.001   0.0001      0.001    0.01    0.00001 ];
%UpperBounds=[0.1       0.1         0.1       0.1      8       0.1       0.1         0.1      0.1          0.1     0.1            0.1     0.1       0.1       0.2    0.01 ];
 




%            [m       beta_A     beta_P    theta_1  epsilon   mu       mu_A        mu_H      gamma    theta_2     sigma       zeta  theta_3      nu       b          P0      A0     H0           R0];
%LowerBounds=[-0.1     0.0001    0.0001     0.001     1     0.0001    0.0001      0.001      0.0001    0.0001     0.001        0.01   0.001      0.001    0.01      0.01   0.0001   0.0001     0.001 ];
%UpperBounds=[0.1       0.001      0.001     0.1      3      0.001     0.001       0.01      0.001      0.001     0.01         0.09   0.01       0.1       0.2       0.1    0.01    0.01        0.01];
 

%Go one by one with altering bounds: first realistic ones that can be
%changed up/down and then others if needed NOTE: P0 seems to be off (make
%it higher and see what happens with rest of things, make it like .15) 
%LowerBounds=[-0.1     0.00001    0.00001    0.001     3     0.0001    0.0001      0.000001     0.00001      0.1       0.0001        0.1   0.2        0.001    0.001      0.1   0.0001   0.0001     0.001 ];
%UpperBounds=[0.1       0.001      0.001     0.1       5      0.01     0.01       0.0001       0.0001       0.2        0.01          0.3   0.3        0.1       0.1       0.3    0.01    0.01        0.01];
 

%Bounds that converge well with minimizing all 4 Diffs 
%LowerBounds=[-0.1     0.000001    0.000001    0.001     5     0.01    0.01      0.000001     0.00001      0.1       0.0001        0.1   0.2        0.001    0.001     0.00001   0.0001   0.0001     0.001 ];
%UpperBounds=[0.1       0.0001       0.0001     0.1      8      0.1     0.1       0.0001       0.0001       0.2        0.01        0.3   0.3        0.1       0.1       0.001      0.2    0.01        0.01];
 
%Pretty good bounds when minimizing Diff4 only 
%LowerBounds=[-0.1  0.1  0.001  0.0001    1     0.001   0.001   0.001   0.001    0.001   0.1      0.01    0.01       0.1       0.1       0.01    0.01     0.01      0.01 ];
%UpperBounds=[0.1   0.4    0.1   0.01     3      0.01    0.01    0.01    0.1     0.1     0.5      0.1      0.1       0.3       0.3       0.3      0.1       0.2         0.2    ];


%Bounds I wanted as realistic but with decreasing order of magnitude range
%LowerBounds=[-0.1     0.0001    0.0001    0.00001      1     0.0001    0.0001      0.0001     0.001      0.001       0.01        0.001   0.001        0.1      0.1     0.001   0.0001   0.0001     0.0001 ];
%UpperBounds=[0.1      0.1        0.1       0.001      5      0.01       0.01        0.01      0.5         0.1        0.5          0.1     0.2        0.5       0.5       0.5      0.2    0.1        0.1];
 
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
numstartpoints=10;

% Runs MultiStart with numstartpoints to find a solution or multiple local solutions to problem; 
% solutions contains the distinct local minima found during the run
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

% x vector to estimate, alpha=m*t+b 
m=x(1);
beta_A=x(2); 
beta_P=x(3); 
theta_1=x(4);
epsilon=x(5);
mu=x(6);  
mu_A=x(7);   
mu_H=x(8);
gamma=x(9);   
theta_2=x(10);
sigma=x(11);
zeta=x(12);
theta_3=x(13);
nu=x(14);
omega=0.0000000001;
b=x(15);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Print optimal parameter solution and objective function value in command
% window when completed 
format short 
x
fval

% Final time N; will run 2013-2018 where t=0 represents 2013
% and t=5 represents 2018, with spacing (T-0)/((N+1)-1)=1 between the points:
N = 24; 
tspan=linspace(0,N,N+1);


% Initial conditions
P0=0.12;%x(16);%0.12;%x(16);%0.0710;
A0=0.003;%x(17);%0.003;%x(17);
H0=0.0006;%x(18);%0.0006;%x(18);%0.0001;%0.00136;%x(7);
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
  alpha=x(1)*t+x(15);
  
   % Making sure S+P+A+H+R=1
  for i=1:N+1
      total(i)=y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5);
  end
  
 %%% For testing purposes: states and corresponding simulated data 
 State1=y(:,1);
 %State_data_1=[0.786760000000000;0.895080745860360;0.901924830346184;0.905069849396411;0.908137177042235;0.911231997643852;0.914378610631427;0.917577035637392;0.920832262684462;0.924140649483092;0.927502296899898;0.930917095573708;0.934384411942796;0.937906705890588;0.941483788276043;0.945116467103051;0.948805724069555;0.952552051009700;0.956355460016689;0.960216644894470;0.964136078269128;0.968114241872428;0.972151789271917;0.976249701711764;0.980408510395841];
 
 
 % Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1))
 %plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 legend('Proportion of susceptibles')%,'Proportion of susceptibles (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
 
 State2=y(:,2);
 %State_data_2=[0.200000000000000;0.0909481525493455;0.0836513997725067;0.0800935206004890;0.0766440206327074;0.0731973067177299;0.0697288104612356;0.0662384437829699;0.0627210220737682;0.0591800019843818;0.0556152256127668;0.0520265690436615;0.0484145481754339;0.0447765513782161;0.0411126188983823;0.0374218111890252;0.0337030212317659;0.0299556183187020;0.0261795036887717;0.0223738480348019;0.0185380936036518;0.0146716575474121;0.0107737896087228;0.00684344682625655;0.00287999678287193];
 
 
 % Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2))
 %plot(t(1:end), State_data_2, 'x')
 xlabel('Year')
 ylabel('Prescription Users')
 legend('Proportion of prescription users')%,'Proportion of prescription users (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})

 State3=y(:,3);
 %State_data_3=[0.00760000000000000;0.00829704878967021;0.00870267430126883;0.00905828329941436;0.00937450301473490;0.00965269200017913;0.00989442412473266;0.0101002338899966;0.0102715450528115;0.0104096884790280;0.0105150275438714;0.0105891001530028;0.0106327393039710;0.0106467927074651;0.0106322619292071;0.0105899570928910;0.0105206317835471;0.0104251672719983;0.0103040776705839;0.0101582082742494;0.00998808627740331;0.00979433685522646;0.00957757634588856;0.00933824323329237;0.00907696847090869];
 
 
 % Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3))
 %plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend('Proportion of opioid addicts')%,'Proportion of opioid addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
 State4=y(:,4);
 %State_data_4=[0.00121000000000000;0.00112681092166828;0.00104927460360314;0.000977059818811785;0.000909745840400286;0.000847058950455061;0.000788626668987576;0.000734318368299567;0.000683788417178302;0.000636722979143746;0.000593041200880438;0.000552393991735842;0.000514584882404989;0.000479448353007741;0.000446755083059968;0.000416339889477986;0.000388055729672108;0.000361713882829795;0.000337232311223003;0.000314434482745690;0.000293231186199367;0.000273506782310314;0.000255146539436483;0.000238080755639956;0.000222191814765932];
 
 % Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4))
 %plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend('Proportion of heroin/fentanyl addicts')%,'Proportion of heroin/fentanyl addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 State5=y(:,5);
 %State_data_5=[0.00443000000000000;0.00454724188048252;0.00467182097860501;0.00480128688737465;0.00493455347269826;0.00507094469061387;0.00520952811771746;0.00534996832772595;0.00549138178122815;0.00563293708715384;0.00577440875872846;0.00591484125637144;0.00605371571489895;0.00619050169074296;0.00632457583293143;0.00645542474411386;0.00658256720267669;0.00670544953236559;0.00682372632661423;0.00693686432594808;0.00704451067421613;0.00714625695172395;0.00724169824189503;0.00733052748007980;0.00741233254232748];
 
 % Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5))
 %plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 legend('Proportion of stably recovered addicts')%,'Proportion of stably recovered addicts (simulated) data' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})

 
 State6=y(:,6);
 %State_data_6=[0;0.254333896046208;0.507557947619061;0.750828745235986;0.983758573114804;1.20628779321856;1.41833707721303;1.61987397892858;1.81081466092476;1.99107303551528;2.16062204429173;2.31935240699286;2.46721054299674;2.60410705347866;2.72997154384002;2.84471979137727;2.94826400541727;3.04052494696663;3.12140763321612;3.19083011612309;3.24869732900780;3.29491832747158;3.32940049486518;3.35204113923119;3.36274919069678];
 
 % Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6))
 %plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend('Proportion that enter P at some point in each quarter')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 State7=y(:,7);
 %State_data_7=[0;0.00100779946271746;0.00174327368150607;0.00244340232074780;0.00311681908985485;0.00376362920275078;0.00438365397230402;0.00497678834646893;0.00554277579184399;0.00608133215598109;0.00659236989809272;0.00707553557113954;0.00753065599310841;0.00795743530786329;0.00835564534480326;0.00872500175669082;0.00906520507909812;0.00937599160623996;0.00965702218400600;0.00990802210443387;0.0101286537004824;0.0103185988417532;0.0104775333294769;0.0106050889170413;0.0107009474056781];
 
 % Simulated data for L and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7))
 %plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend('Proportion that enter A at some point in each quarter')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 
 
 State8=y(:,8);
 %State_data_8=[0;4.26501816191890e-06;8.16701023874181e-06;1.17859981066272e-05;1.51647113241712e-05;1.83272545883561e-05;2.13003873839925e-05;2.40927260868802e-05;2.67275118175865e-05;2.92226007008250e-05;3.15819745663623e-05;3.38239691802428e-05;3.59601048702499e-05;3.79943060596506e-05;3.99389057739059e-05;4.18001749552355e-05;4.35828254466968e-05;4.52953475429103e-05;4.69396975929900e-05;4.85231546324224e-05;5.00491557095012e-05;5.15221483446924e-05;5.29467906066497e-05;5.43265077898365e-05;5.56658925104229e-05];
 
 % Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8))
 %plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend('Proportion that enter H at some point in each quarter')%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24])
 set(gca, 'fontsize',10)
 xtickangle(90)
 set(gca,'XLim',[0 N])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018', 'Q1 2019'})
 %%%
 
 %%% Data points we are interested in 


 % Yearly simulation of individuals in P class throughout year for years
 % 2013-2017
 Estim1=[sum(y(1:4,2)+y(2:5,6)-y(1:4,6)); sum(y(5:8,2)+y(6:9,6)-y(5:8,6)); sum(y(9:12,2)+y(10:13,6)-y(9:12,6));...
         sum(y(13:16,2)+y(14:17,6)-y(13:16,6)); sum(y(17:20,2)+y(18:21,6)-y(17:20,6))];
 
 % Actual Data for years 2013-2017
 Data1=[1825910./5517176; 1805325./5559006; 1800613./5602117; 1744766./5651993; 1620951./5708586];
 %Data1=[1.43845164603714;1.11286466940459;0.885938700786563;0.652778992061590;0.412645314864569];
 
 % Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 scatter(0:1:4, Estim1,'o')
 scatter(0:1:4, Data1,'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
 

 % Yearly simulation of individuals in A class throughout year for years
 % 2013-2017
 Estim2=[sum(y(1:4,3)+y(2:5,7)-y(1:4,7)); sum(y(5:8,3)+y(6:9,7)-y(5:8,7)); sum(y(9:12,3)+y(10:13,7)-y(9:12,7));...
         sum(y(13:16,3)+y(14:17,7)-y(13:16,7)); sum(y(17:20,3)+y(18:21,7)-y(17:20,7))];
 
 
 % Actual Data for years 2013-2017 
 Data2=[43418./5517176; 42928./5559006; 42816./5602117; 37464./5651993; 34816./5708586];
 %Data2=[0.0367748254802082;0.0414478097316324;0.0437732414299781;0.0440363001195239;0.0424715336217629];
 
 % Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 scatter(0:1:4, Estim2,'o')
 scatter(0:1:4, Data2,'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 3 4])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017'})
 


 % Yearly simulation of individuals in H class throughout year for years
 % 2014-2016
 
 Estim3=[sum(y(5:8,4)+y(6:9,8)-y(5:8,8)); sum(y(9:12,4)+y(10:13,8)-y(9:12,8));...
         sum(y(13:16,4)+y(14:17,8)-y(13:16,8))];
 
 % Actual Data for years 2014-2016
 Data3=[7560./5559006; 7560./5602117; 10260./5651993];
 %Data3=[0.00329131262863590;0.00247517918199099;0.00186475092852713];
 
  
 % Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(11)
 hold all
 scatter(1:1:3, Estim3,'o')
 scatter(1:1:3, Data3,'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('ODE solution', 'Data')
 set(gca, 'xtick', [ 0 1 2 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 
 
 Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);


 
 % Actual Data for years 2013-2018 INCLUDES addicts 
 %Data4=[856000./5517176; 863000./5517176; 870000./5517176; 850000./5517176;...
        %830000./5559006; 860000./5559006; 871000./5559006; 840000./5559006;...
        %815000./5602117; 855000./5602117; 856000./5602117; 850000./5602117;...
        %820000./5651993; 810000./5651993; 805000./5651993; 800000./5651993;...
        %790000./5708586; 780000./5708586; 750000./5708586; 713000./5708586;...
        %698000./5779971; 695000./5779971; 650000./5779971; 631000./5779971];
        
 %Actual Data for years 2013-2018 EXCLUDES addicts 
 Data4=[836027./5517176; 843027./5517176; 850027./5517176; 830027./5517176;...
        810038./5559006; 840038./5559006; 851038./5559006; 820038./5559006;...
        795305./5602117; 835305./5602117; 836305./5602117; 830305./5602117;...
        802767./5651993; 792767./5651993; 787767./5651993; 782767./5651993;...
        769633./5708586; 759633./5708586; 729633./5708586; 692633./5708586;...
        679000./5779971; 676000./5779971; 631000./5779971; 612000./5779971];
 %Data4= [0.454333896046208;0.344172204122199;0.326922197389432;0.313023348479306;0.299173240736468;0.285246590712193;0.271265712176788;0.257179125779146;0.242979396664295;0.228729010760832;0.214345588313890;0.199884705047545;0.185311058657349;0.170641041739577;0.155860866435641;0.140966025229023;0.125963962781128;0.110838304568182;0.0956019865957489;0.0802410609195099;0.0647590920674275;0.0491538249410182;0.0334144339747313;0.0175514982918497];
 
 
 figure(12)
 hold all
 plot(t(1:end-1),Estim4, 'o')
 plot(t(1:end-1), Data4, 'x')
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
 
 
 
function value = HeroinModel_ODE15s(z)

%Parameters
m=z(1);
beta_A=z(2); 
beta_P=z(3); 
theta_1=z(4);
epsilon=z(5);
mu=z(6);  
mu_A=z(7);   
mu_H=z(8);
gamma=z(9);   
theta_2=z(10); 
sigma=z(11);
zeta=z(12);
theta_3=z(13);
nu=z(14);
omega=0.0000000001;
b=z(15);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Final time N; will run 2013-2018 where t=0 represents 2013
% and t=5 represents 2018, with spacing (T-0)/((N+1)-1)=1 between the points
N = 24; 
tspan=linspace(0,N,N+1);

% Initial conditions
P0=0.12;%z(16);%0.12;%x(16);%0.0710;
A0=0.003;%z(17);%0.0035;%x(17);
H0=0.0006;%z(18);%0.0006;%x(18);%0.0001;%0.00136;%x(7);
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

 Estim1=[sum(y(1:4,2)+y(2:5,6)-y(1:4,6)); sum(y(5:8,2)+y(6:9,6)-y(5:8,6)); sum(y(9:12,2)+y(10:13,6)-y(9:12,6));...
         sum(y(13:16,2)+y(14:17,6)-y(13:16,6)); sum(y(17:20,2)+y(18:21,6)-y(17:20,6))];
 


 % Actual proportions of population (updated 3/12/19) that were non-addicted prescription opioid users at some point
 % during the year for 2013-2017 
 % (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 
 Data1=[1825910./5517176; 1805325./5559006; 1800613./5602117; 1744766./5651993; 1620955./5708586];
 %Data1=[1.43845164603714;1.11286466940459;0.885938700786563;0.652778992061590;0.412645314864569];
 
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
 
  Estim2=[sum(y(1:4,3)+y(2:5,7)-y(1:4,7)); sum(y(5:8,3)+y(6:9,7)-y(5:8,7)); sum(y(9:12,3)+y(10:13,7)-y(9:12,7));...
         sum(y(13:16,3)+y(14:17,7)-y(13:16,7)); sum(y(17:20,3)+y(18:21,7)-y(17:20,7))];
 
 % When testing all points with simulated data
 %Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 % Actual proportions of population (updated 3/12/19) that were opioid addicted individuals in
 % the population at some point during the year in 2014 and 2015 
 % (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 Data2=[43418./5517176; 42928./5559006; 42816./5602117; 37464./5651993; 34805./5708586];
 %Data2=[0.0367748254802082;0.0414478097316324;0.0437732414299781;0.0440363001195239;0.0424715336217629];
 
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

  
 Estim3=[sum(y(5:8,4)+y(6:9,8)-y(5:8,8)); sum(y(9:12,4)+y(10:13,8)-y(9:12,8));...
         sum(y(13:16,4)+y(14:17,8)-y(13:16,8))];
 

 
 % Actual proportion (updated 3/12/19) of heroin addicted individuals in the population at some point during the year 
 % in 2014, 2015, and 2016
 % (total number of heroin addicted individuals in 2014, 2015, and 2016 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 
 Data3=[7560./5559006; 7560./5602117; 10260./5651993];
 %Data3=[0.00329131262863590;0.00247517918199099;0.00186475092852713];
 
  
 % Data simulated when testing codes 
 %Data3=[0.00116527288223448;0.00120952017524577;0.00118883157707289];
 
 % The difference between estimated value and data
 
 Diff3=Estim3-Data3;

 
  
 Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);
 %Data1=[0.399384466780766;0.476721593432771;0.469765087997695;0.449866817308604;0.427737148099884];
 %FILL IN 2018 total pop/data points and COMPARE TO YEARS that we know and
 %use 2nd graph for 2014 on!!
 
 
 % Actual Data for years 2013-2018 INDCLUDES ADDICTS
 %Data4=[856000./5517176; 863000./5517176; 870000./5517176; 850000./5517176;...
        %830000./5559006; 860000./5559006; 871000./5559006; 840000./5559006;...
        %815000./5602117; 855000./5602117; 856000./5602117; 850000./5602117;...
        %820000./5651993; 810000./5651993; 805000./5651993; 800000./5651993;...
        %790000./5708586; 780000./5708586; 750000./5708586; 713000./5708586;...
        %698000./5779971; 695000./5779971; 650000./5779971; 631000./5779971];
    
     %Actual Data for years 2013-2018 EXCLUDES addicts 
 Data4=[836027./5517176; 843027./5517176; 850027./5517176; 830027./5517176;...
        810038./5559006; 840038./5559006; 851038./5559006; 820038./5559006;...
        795305./5602117; 835305./5602117; 836305./5602117; 830305./5602117;...
        802767./5651993; 792767./5651993; 787767./5651993; 782767./5651993;...
        769633./5708586; 759633./5708586; 729633./5708586; 692633./5708586;...
        679000./5779971; 676000./5779971; 631000./5779971; 612000./5779971];
 %Data4=[0.454333896046208;0.344172204122199;0.326922197389432;0.313023348479306;0.299173240736468;0.285246590712193;0.271265712176788;0.257179125779146;0.242979396664295;0.228729010760832;0.214345588313890;0.199884705047545;0.185311058657349;0.170641041739577;0.155860866435641;0.140966025229023;0.125963962781128;0.110838304568182;0.0956019865957489;0.0802410609195099;0.0647590920674275;0.0491538249410182;0.0334144339747313;0.0175514982918497];
 
 
 Diff4=Estim4-Data4;
 
 
 
                  
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
 
 
 
 

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.786760000000000;0.895080745860360;0.901924830346184;0.905069849396411;0.908137177042235;0.911231997643852;0.914378610631427;0.917577035637392;0.920832262684462;0.924140649483092;0.927502296899898;0.930917095573708;0.934384411942796;0.937906705890588;0.941483788276043;0.945116467103051;0.948805724069555;0.952552051009700;0.956355460016689;0.960216644894470;0.964136078269128;0.968114241872428;0.972151789271917;0.976249701711764;0.980408510395841];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.200000000000000;0.0909481525493455;0.0836513997725067;0.0800935206004890;0.0766440206327074;0.0731973067177299;0.0697288104612356;0.0662384437829699;0.0627210220737682;0.0591800019843818;0.0556152256127668;0.0520265690436615;0.0484145481754339;0.0447765513782161;0.0411126188983823;0.0374218111890252;0.0337030212317659;0.0299556183187020;0.0261795036887717;0.0223738480348019;0.0185380936036518;0.0146716575474121;0.0107737896087228;0.00684344682625655;0.00287999678287193];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00760000000000000;0.00829704878967021;0.00870267430126883;0.00905828329941436;0.00937450301473490;0.00965269200017913;0.00989442412473266;0.0101002338899966;0.0102715450528115;0.0104096884790280;0.0105150275438714;0.0105891001530028;0.0106327393039710;0.0106467927074651;0.0106322619292071;0.0105899570928910;0.0105206317835471;0.0104251672719983;0.0103040776705839;0.0101582082742494;0.00998808627740331;0.00979433685522646;0.00957757634588856;0.00933824323329237;0.00907696847090869];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.00121000000000000;0.00112681092166828;0.00104927460360314;0.000977059818811785;0.000909745840400286;0.000847058950455061;0.000788626668987576;0.000734318368299567;0.000683788417178302;0.000636722979143746;0.000593041200880438;0.000552393991735842;0.000514584882404989;0.000479448353007741;0.000446755083059968;0.000416339889477986;0.000388055729672108;0.000361713882829795;0.000337232311223003;0.000314434482745690;0.000293231186199367;0.000273506782310314;0.000255146539436483;0.000238080755639956;0.000222191814765932];
 
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.00443000000000000;0.00454724188048252;0.00467182097860501;0.00480128688737465;0.00493455347269826;0.00507094469061387;0.00520952811771746;0.00534996832772595;0.00549138178122815;0.00563293708715384;0.00577440875872846;0.00591484125637144;0.00605371571489895;0.00619050169074296;0.00632457583293143;0.00645542474411386;0.00658256720267669;0.00670544953236559;0.00682372632661423;0.00693686432594808;0.00704451067421613;0.00714625695172395;0.00724169824189503;0.00733052748007980;0.00741233254232748];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.254333896046208;0.507557947619061;0.750828745235986;0.983758573114804;1.20628779321856;1.41833707721303;1.61987397892858;1.81081466092476;1.99107303551528;2.16062204429173;2.31935240699286;2.46721054299674;2.60410705347866;2.72997154384002;2.84471979137727;2.94826400541727;3.04052494696663;3.12140763321612;3.19083011612309;3.24869732900780;3.29491832747158;3.32940049486518;3.35204113923119;3.36274919069678];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.00100779946271746;0.00174327368150607;0.00244340232074780;0.00311681908985485;0.00376362920275078;0.00438365397230402;0.00497678834646893;0.00554277579184399;0.00608133215598109;0.00659236989809272;0.00707553557113954;0.00753065599310841;0.00795743530786329;0.00835564534480326;0.00872500175669082;0.00906520507909812;0.00937599160623996;0.00965702218400600;0.00990802210443387;0.0101286537004824;0.0103185988417532;0.0104775333294769;0.0106050889170413;0.0107009474056781];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;4.26501816191890e-06;8.16701023874181e-06;1.17859981066272e-05;1.51647113241712e-05;1.83272545883561e-05;2.13003873839925e-05;2.40927260868802e-05;2.67275118175865e-05;2.92226007008250e-05;3.15819745663623e-05;3.38239691802428e-05;3.59601048702499e-05;3.79943060596506e-05;3.99389057739059e-05;4.18001749552355e-05;4.35828254466968e-05;4.52953475429103e-05;4.69396975929900e-05;4.85231546324224e-05;5.00491557095012e-05;5.15221483446924e-05;5.29467906066497e-05;5.43265077898365e-05;5.56658925104229e-05];
 
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
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5);
 
 % Objective function value we wish to minimize; want value=fval(x) to be small  when run MultiStart
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4);
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 %value=norm(Diff4,2)./norm(Data4);



 
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

