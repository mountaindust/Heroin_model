%File name: heroin_multistart_MLE_parameter_estimation.m 

%clf;
%clear all;

% Realistic parameter bounds
%           [m      betaA     betaP   theta1   epsilon  gamma   theta2   sigma    zeta    theta3    nu       b     P0        A0       H0       R0       c      d         e       p     var_1    p2]   
LowerBounds=[-0.1  0.00001   0.00001   0.05      1     0.005     0.1     0.1    0.0001     10      0.0001   0.1   0.001   0.0001   0.00001  0.00001   -0.1   0.0001   0.001    0.001  0.000001   0.001];  
UpperBounds=[-0.001  0.01     0.01      0.3      5       0.1     0.4     2       0.2      20       0.2     0.5    0.7      0.05      0.01     0.5    -0.001  0.008     0.1     0.99    0.1       0.99];  

% Initial starting points for parameters, starting in the middle of each of the ranges
%Will crash 
%xstart=0.5*(LowerBounds + UpperBounds); 
%Selected randomly in intervals just to get running, works
xstart=[-0.01 0.001 0.0001 0.1 1.1 0.08 0.2 0.2 0.1 15 0.2 0.49 0.2 0.001 0.0001 0.0001 -0.002 0.0003 0.01 0.1 0.001 0.1];

% Create MultiStart problem using optimization function fmincon;
% x0 is xstart, objective is what we are trying to maximize which comes from 
% value = HeroinModel_ODE45(z) = fval(x) as output
problem=createOptimProblem('fmincon','objective',@HeroinModel_ODE15s,...
         'x0', xstart,...
         'lb',LowerBounds,...
         'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

% Define a multistart problem; results are reported after each local solver run, as well as final summary
ms=MultiStart('Display', 'iter'); 

% Number of times I want to run optimization scheme
numstartpoints=20000;

%Set up parallel processing for hermione, parpool(32). If doing parallel processing on
%local computer, do parpool. 
parpool(32)


% Runs MultiStart with numstartpoints to find a solution or multiple local solutions to problem; 
% solutions contains the distinct local minima found during the run
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

% x vector to estimate, alpha=m*t+b then slope of c starting at 2016Q2, mu_A=d*t+e
m=x(1);
beta_A=x(2);
beta_P=x(3);
theta_1=x(4);
epsilon=x(5);
mu=0.00710; 
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
d=x(18);
e=x(19);
p=x(20);
var_1=x(21);
p2=x(22);

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e,p,var_1,p2];


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

%Plotting section
% Yearly simulation of individuals in P class at all during the year for years
 % 2013-2017
 
%  Estim1=[pars(21)*(y(1,2)+y(5,6)-y(1,6)); pars(21)*(y(5,2)+y(9,6)-y(5,6)); pars(21)*(y(9,2)+y(13,6)-y(9,6));...
%          pars(21)*(y(13,2)+y(17,6)-y(13,6)); pars(21)*(y(17,2)+y(21,6)-y(17,6)); pars(21)*(y(21,2)+y(25,6)-y(21,6))];
%      
%  % Actual Data for years 2013-2018
%  Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
%  
%  %Testing Data 
%  %Data1=[0.331992810413221;0.323700348480663;0.321231125203553;0.312639103866186;0.285267888975045;0.251873011646405];
%  
%  % Data points from proportion that is in P at some point in the year and corresponding ODE solution points 
%  figure(1)
%  hold all
%  z1 = linspace(0,5,6); %defines mesh where going to plot Estim1, Data1 values 
%  scatter(z1, Estim1, 100,'o');
%  scatter(z1, Data1, 100,'x');
%  set(gca, 'fontsize',10)
%  xlabel('Year')
%  ylabel('Proportion in P') % at some point during the year
%  legend({'Model simulation', 'Data'},'FontSize', 14)
%  set(gca, 'xtick', [ 0 1 2 3 4 5 ])
%  set(gca, 'fontsize',10)
%  set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
%  
% 
%  % Yearly simulation of individuals in A class at all during the year for years
%  % 2013-2017
%  Estim2=[pars(19)*(y(1,3)+y(5,7)-y(1,7)); pars(19)*(y(5,3)+y(9,7)-y(5,7)); pars(19)*(y(9,3)+y(13,7)-y(9,7));...
%         pars(19)*(y(13,3)+y(17,7)-y(13,7)); pars(19)*(y(17,3)+y(21,7)-y(17,7)); pars(19)*(y(21,3)+y(25,7)-y(21,7))];
%  
%  
%  % Actual Data for years 2013-2018
%  Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
%  
%  %Testing Data
%  %Data2=[0.00798220556723511;0.00493886234233757;0.00323004248927816;0.00206292144694655;0.00124601472385274;0.000729472026712258];
%  
%  
%   % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
%  figure(2)
%  hold all
%  z2 = linspace(0,5,6);
%  scatter(z2, Estim2, 100,'o');
%  scatter(z2, Data2, 100,'x');
%  set(gca, 'fontsize',10)
%  xlabel('Year')
%  ylabel('Proportion in A') % at some point during the year
%  legend({'Model simulation', 'Data'},'FontSize', 14) 
%  set(gca, 'xtick', [ 0 1 2 3 4 5 ])
%  set(gca, 'fontsize',10)
%  set(gca,'xticklabel',{'2013','2014','2015','2016', '2017', '2018'})
% 
% 
% 
%  % Yearly simulation of individuals in H class at all during the year for years
%  % 2014-2016
%  
%  Estim3=[pars(19)*(y(5,4)+y(9,8)-y(5,8)); pars(19)*(y(9,4)+y(13,8)-y(9,8)); pars(19)*(y(13,4)+y(17,8)-y(13,8))];
%  
%  % Actual Data for years 2014-2016
%  Data3=[7560./5559702; 7560./5602187; 10260./5648259];
%  
%  %Testing Data
%  %Data3=[0.00133732150290409;0.00231840400714378;0.00377884713025854];
%  
%  
%  % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
%  figure(3)
%  hold all
%  z3 = linspace(0,2,3);
%  scatter(z3, Estim3, 100,'o');
%  scatter(z3, Data3, 100,'x');
%  set(gca, 'fontsize',10)
%  xlabel('Year')
%  ylabel('Proportion in H') % at some point during the year
%  legend({'Model simulation', 'Data'},'FontSize', 14)
%  set(gca, 'xtick', [ 0 1 2 ])
%  set(gca, 'fontsize',10)
%  set(gca,'xticklabel',{'2014', '2015', '2016'})
% 
%  
%  
%  % Yearly simulation of individuals in P class at all during the quarters in years
%  % 2013-2018
%  Estim4=pars(21)*(y(1:24,2)+y(2:25,6)-y(1:24,6));
% 
%         
%  %Actual Data for years 2013-2018 
%  Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
%         833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
%         827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
%         832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
%         775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
%         688451./5754509; 683498./5754509; 641894./5754509; 625054./5754509];
%  
%  %Testing Data
% %Data4=[0.159345888525904;0.155369684395518;0.153226232495256;0.152038272755859;0.151303670490909;0.150774641564431;0.150353591755673;0.149995131706456;0.149672168321776;0.149358805511151;0.149050068737102;0.148734280370870;0.148085068864361;0.147016117783976;0.145043027658001;0.142252916511528;0.138990114824242;0.135356694964831;0.131519862311692;0.127622405412348;0.123717988037982;0.119806623029534;0.115881867009583;0.111940266417211];
%  
%  % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
%  figure(4)
%  hold all
%  z4 = linspace(0,23,24);
%  scatter(z4, Estim4, 100,'o');
%  scatter(z4, Data4, 100,'x');
%  set(gca, 'fontsize',10)
%  xlabel('Quarter')
%  ylabel('Proportion in P') % at some point during the year
%  legend({'Model simulation', 'Data'},'FontSize', 14)
%  set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23])
%  set(gca, 'fontsize',10)
%  xtickangle(90)
%  set(gca,'XLim',[0 23])
%  set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
%                        'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
%                        'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
%                        'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
%                        'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
%                        'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
%                    
%                    
%                    
%                    
% % Yearly simulation of individuals overdosing from A class during the year for years
% % 2013-2017 (anyone in A class throughout the year times mu_A)
%  
%  Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];
%      
%  % Actual Data for years 2013-2017
%  Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
%  
%  %Testing Data
%  %Data5=[6.93893503425710e-05;6.66363564939064e-05;6.36911749420292e-05;6.03701269735515e-05;5.65022143567912e-05];
%  
%  
%   % Data points from proportion that is in A at some point and overdoses in the year and corresponding ODE solution points 
%  figure(5)
%  hold all
%  z5 = linspace(0,3,4); %defines mesh where going to plot Estim5, Data5 values 
%  scatter(z5, Estim5, 100,'o');
%  scatter(z5, Data5, 100,'x');
%  set(gca, 'fontsize',10)
%  xlabel('Year')
%  ylabel('Proportion overdose from A') % at some point during the year
%  legend({'Model simulation', 'Data'},'FontSize', 14)
%  set(gca, 'xtick', [ 0 1 2 3])
%  set(gca, 'fontsize',10)
%  set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})                   
%                    
%  
% % Yearly simulation of individuals overdosing from H class during the year for years
% % 2013-2017 (anyone in H class throughout the year times mu_H)
%  
%  Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
%         y(17,10)-y(13,10); y(21,10)-y(17,10)];
%     
%  % Actual Data for years 2013-2017
%  Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
%  
%  %Testing Data
%  %Data6=[2.99015652964405e-05;4.35881487901691e-05;6.30645587490766e-05;9.02990840950082e-05;0.000126571634010397];
%  
%   % Data points from proportion that is in H at some point and overdoses in the year and corresponding ODE solution points 
%  figure(6)
%  hold all
%  z6 = linspace(0,4,5); %defines mesh where going to plot Estim6, Data6 values 
%  scatter(z6, Estim6, 100,'o');
%  scatter(z6, Data6, 100,'x');
%  set(gca, 'fontsize',10)
%  xlabel('Year')
%  ylabel('Proportion overdose from H') % at some point during the year
%  legend({'Model simulation', 'Data'},'FontSize', 14)
%  set(gca, 'xtick', [ 0 1 2 3 4 ])
%  set(gca, 'fontsize',10)
%  set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})    
%  
%  
% W=0:0.0001:1;
% A1 = betapdf(W,pars(19)*((y(1,3)+y(5,7)-y(1,7))*5519417-1),(1-pars(19))*((y(5,3)+y(9,7)-y(5,7))*5519417-1));
% A2 = betapdf(W,pars(19)*((y(5,3)+y(9,7)-y(5,7))*5559702-1),(1-pars(19))*((y(5,3)+y(9,7)-y(5,7))*5559702-1));
% A3 = betapdf(W,pars(19)*((y(9,3)+y(13,7)-y(9,7))*5602187-1),(1-pars(19))*((y(9,3)+y(13,7)-y(9,7))*5602187-1));
% A4 = betapdf(W,pars(19)*((y(13,3)+y(17,7)-y(13,7))*5648259-1),(1-pars(19))*((y(13,3)+y(17,7)-y(13,7))*5648259-1));
% A5 = betapdf(W,pars(19)*((y(17,3)+y(21,7)-y(17,7))*5702475-1),(1-pars(19))*((y(17,3)+y(21,7)-y(17,7))*5702475-1));
% A6 = betapdf(W,pars(19)*((y(21,3)+y(25,7)-y(21,7))*5754509-1),(1-pars(19))*((y(21,3)+y(25,7)-y(21,7))*5754509-1));
% 
% %Plotting A beta distributions for each year
% figure(7)
% plot(W,A1,'Color','g','LineWidth',2)
% hold on
% plot(W,A2,'LineStyle','-.','Color','b','LineWidth',2)
% plot(W,A3,'LineStyle',':','Color','r','LineWidth',2)
% plot(W,A4,'LineStyle','-','Color','m','LineWidth',2)
% plot(W,A5,'LineStyle',':','Color','y','LineWidth',2)
% plot(W,A6,'LineStyle',':','LineWidth',2)
% legend('A proportion in 2013','A proportion in 2014', 'A proportion in 2015','A proportion in 2016','A proportion in 2017','A proportion in 2018','Location','NorthEast');
% hold off
% 
% 
% %Plotting H beta distributions for each year
% W=0:0.0001:1;
% H1 = betapdf(W,pars(19)*((y(5,4)+y(9,8)-y(5,8))*5559702-1),(1-pars(19))*((y(5,4)+y(9,8)-y(5,8))*5559702-1));
% H2 = betapdf(W,pars(19)*((y(9,4)+y(13,8)-y(9,8))*5602187-1),(1-pars(19))*((y(9,4)+y(13,8)-y(9,8))*5602187-1));
% H3 = betapdf(W,pars(19)*((y(13,4)+y(17,8)-y(13,8))*5648259-1),(1-pars(19))*((y(13,4)+y(17,8)-y(13,8))*5648259-1));
% 
% figure(8)
% plot(W,H1,'Color','g','LineWidth',2)
% hold on
% plot(W,H2,'LineStyle','-.','Color','b','LineWidth',2)
% plot(W,H3,'LineStyle',':','Color','r','LineWidth',2)
% legend('H proportion in 2014','H proportion in 2015', 'H proportion in 2016','Location','NorthEast');
% hold off
% 
% 
% %Plotting A overdose gamma distributions for each year 
% %Overdoses added together quarters 1-4 for 2013
% overdoses_year1 = (pars(17)*1+pars(18))*(y(2,9)-y(1,9))+(pars(17)*2+pars(18))*(y(3,9)-y(2,9))+(pars(17)*3+pars(18))*(y(4,9)-y(3,9))+(pars(17)*4+pars(18))*(y(5,9)-y(4,9));
% %Overdoses added together quarters 5-8 for 2014
% overdoses_year2 = (pars(17)*5+pars(18))*(y(6,9)-y(5,9))+(pars(17)*6+pars(18))*(y(7,9)-y(6,9))+(pars(17)*7+pars(18))*(y(8,9)-y(7,9))+(pars(17)*8+pars(18))*(y(9,9)-y(8,9));
% %Overdoses added together quarters 9-12 for 2015
% overdoses_year3 = (pars(17)*9+pars(18))*(y(10,9)-y(9,9))+(pars(17)*10+pars(18))*(y(11,9)-y(10,9))+(pars(17)*11+pars(18))*(y(12,9)-y(11,9))+(pars(17)*12+pars(18))*(y(13,9)-y(12,9));
% %Overdoses added together quarters 13-16 for 2016
% overdoses_year4 = (pars(17)*13+pars(18))*(y(14,9)-y(13,9))+(pars(17)*14+pars(18))*(y(15,9)-y(14,9))+(pars(17)*15+pars(18))*(y(16,9)-y(15,9))+(pars(17)*16+pars(18))*(y(17,9)-y(16,9));
% 
% Aoverdose1 = gampdf(W,(overdoses_year1).^2./(pars(20)^2),(pars(20)^2)./overdoses_year1);
% Aoverdose2 = gampdf(W,(overdoses_year2).^2./(pars(20)^2),(pars(20)^2)./overdoses_year2);
% Aoverdose3 = gampdf(W,(overdoses_year3).^2./(pars(20)^2),(pars(20)^2)./overdoses_year3);
% Aoverdose4 = gampdf(W,(overdoses_year4).^2./(pars(20)^2),(pars(20)^2)./overdoses_year4);
% 
% figure(9)
% plot(W,Aoverdose1,'Color','g','LineWidth',2)
% hold on
% plot(W,Aoverdose2,'LineStyle','-.','Color','b','LineWidth',2)
% plot(W,Aoverdose3,'LineStyle',':','Color','r','LineWidth',2)
% plot(W,Aoverdose4,'LineStyle',':','LineWidth',2)
% legend('A overdose proportion in 2013','A overdose proportion in 2015', 'A overdose proportion in 2016','A overdose proportion in 2017','Location','NorthEast');
% hold off



function value = HeroinModel_ODE15s(z)

% Parameters
m=z(1);
beta_A=z(2);
beta_P=z(3);
theta_1=z(4);
epsilon=z(5);
mu=0.00710; 
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
d=z(18);
e=z(19);
p=z(20);
var_1=z(21);
p2=z(22);

% Parameter vector
pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e,p,var_1,p2];

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
  

% COMPARING MODEL ESTIMATES TO DATA-distribution functions for each type of data, evaluated at each data point (48 total) 

%Prescribed users yearly (non-addicted), 2013-2018: beta distributed 
P_prop1 = betalike([pars(21)*((y(1,2)+y(5,6)-y(1,6))*5519417-1),(1-pars(21))*((y(1,2)+y(5,6)-y(1,6))*5519417-1)],1825910./((y(1,2)+y(5,6)-y(1,6))*5519417));
P_prop2 = betalike([pars(21)*((y(5,2)+y(9,6)-y(5,6))*5559702-1),(1-pars(21))*((y(5,2)+y(9,6)-y(5,6))*5559702-1)],1805325./((y(5,2)+y(9,6)-y(5,6))*5559702));
P_prop3 = betalike([pars(21)*((y(9,2)+y(13,6)-y(9,6))*5602187-1),(1-pars(21))*((y(9,2)+y(13,6)-y(9,6))*5602187-1)],1800614./((y(9,2)+y(13,6)-y(9,6))*5602187));
P_prop4 = betalike([pars(21)*((y(13,2)+y(17,6)-y(13,6))*5648259-1),(1-pars(21))*((y(13,2)+y(17,6)-y(13,6))*5648259-1)],1744766./((y(13,2)+y(17,6)-y(13,6))*5648259));
P_prop5 = betalike([pars(21)*((y(17,2)+y(21,6)-y(17,6))*5702475-1),(1-pars(21))*((y(17,2)+y(21,6)-y(17,6))*5702475-1)],1620955./((y(17,2)+y(21,6)-y(17,6))*5702475));
P_prop6 = betalike([pars(21)*((y(21,2)+y(25,6)-y(21,6))*5754509-1),(1-pars(21))*((y(21,2)+y(25,6)-y(21,6))*5754509-1)],1455093./((y(21,2)+y(25,6)-y(21,6))*5754509)); 

%Prescription opioid addicts yearly (non heroin-addicted), 2013-2018: beta distributed
A_prop1 = betalike([pars(19)*((y(1,3)+y(5,7)-y(1,7))*5519417-1),(1-pars(19))*((y(5,3)+y(9,7)-y(5,7))*5519417-1)],43418./((y(1,3)+y(5,7)-y(1,7))*5519417));
A_prop2 = betalike([pars(19)*((y(5,3)+y(9,7)-y(5,7))*5559702-1),(1-pars(19))*((y(5,3)+y(9,7)-y(5,7))*5559702-1)],42928./((y(5,3)+y(9,7)-y(5,7))*5559702));
A_prop3 = betalike([pars(19)*((y(9,3)+y(13,7)-y(9,7))*5602187-1),(1-pars(19))*((y(9,3)+y(13,7)-y(9,7))*5602187-1)],42816./((y(9,3)+y(13,7)-y(9,7))*5602187));
A_prop4 = betalike([pars(19)*((y(13,3)+y(17,7)-y(13,7))*5648259-1),(1-pars(19))*((y(13,3)+y(17,7)-y(13,7))*5648259-1)],37464./((y(13,3)+y(17,7)-y(13,7))*5648259));
A_prop5 = betalike([pars(19)*((y(17,3)+y(21,7)-y(17,7))*5702475-1),(1-pars(19))*((y(17,3)+y(21,7)-y(17,7))*5702475-1)],34805./((y(17,3)+y(21,7)-y(17,7))*5702475));
A_prop6 = betalike([pars(19)*((y(21,3)+y(25,7)-y(21,7))*5754509-1),(1-pars(19))*((y(21,3)+y(25,7)-y(21,7))*5754509-1)],31244./((y(21,3)+y(25,7)-y(21,7))*5754509));


%Heroin addicts yearly, 2014-2016: beta distributed 
H_prop1 = betalike([pars(19)*((y(5,4)+y(9,8)-y(5,8))*5559702-1),(1-pars(19))*((y(5,4)+y(9,8)-y(5,8))*5559702-1)],7560./((y(5,4)+y(9,8)-y(5,8))*5559702));
H_prop2 = betalike([pars(19)*((y(9,4)+y(13,8)-y(9,8))*5602187-1),(1-pars(19))*((y(9,4)+y(13,8)-y(9,8))*5602187-1)],7560./((y(9,4)+y(13,8)-y(9,8))*5602187));
H_prop3 = betalike([pars(19)*((y(13,4)+y(17,8)-y(13,8))*5648259-1),(1-pars(19))*((y(13,4)+y(17,8)-y(13,8))*5648259-1)],10260./((y(13,4)+y(17,8)-y(13,8))*5648259));

%Prescribed opioid users quarterly, 2013Q1-2018Q4: beta distributed 
%2013 Q1-Q4
PQ_prop1 = betalike([pars(21)*((y(1,2)+y(2,6)-y(1,6))*5519417-1),(1-pars(21))*((y(1,2)+y(2,6)-y(1,6))*5519417-1)],847077./((y(1,2)+y(2,6)-y(1,6))*5519417));
PQ_prop2 = betalike([pars(21)*((y(2,2)+y(3,6)-y(2,6))*5519417-1),(1-pars(21))*((y(2,2)+y(3,6)-y(2,6))*5519417-1)],860931./((y(2,2)+y(3,6)-y(2,6))*5519417));
PQ_prop3 = betalike([pars(21)*((y(3,2)+y(4,6)-y(3,6))*5519417-1),(1-pars(21))*((y(3,2)+y(4,6)-y(3,6))*5519417-1)],864889./((y(3,2)+y(4,6)-y(3,6))*5519417));
PQ_prop4 = betalike([pars(21)*((y(4,2)+y(5,6)-y(4,6))*5519417-1),(1-pars(21))*((y(4,2)+y(5,6)-y(4,6))*5519417-1)],847077./((y(4,2)+y(5,6)-y(4,6))*5519417));
%2014 Q1-Q4
PQ_prop5 = betalike([pars(21)*((y(5,2)+y(6,6)-y(5,6))*5559702-1),(1-pars(21))*((y(5,2)+y(6,6)-y(5,6))*5559702-1)],833223./((y(5,2)+y(6,6)-y(5,6))*5559702));
PQ_prop6 = betalike([pars(21)*((y(6,2)+y(7,6)-y(6,6))*5559702-1),(1-pars(21))*((y(6,2)+y(7,6)-y(6,6))*5559702-1)],851035./((y(6,2)+y(7,6)-y(6,6))*5559702));
PQ_prop7 = betalike([pars(21)*((y(7,2)+y(8,6)-y(7,6))*5559702-1),(1-pars(21))*((y(7,2)+y(8,6)-y(7,6))*5559702-1)],861921./((y(7,2)+y(8,6)-y(7,6))*5559702));
PQ_prop8 = betalike([pars(21)*((y(8,2)+y(9,6)-y(8,6))*5559702-1),(1-pars(21))*((y(8,2)+y(9,6)-y(8,6))*5559702-1)],841140./((y(8,2)+y(9,6)-y(8,6))*5559702));
%2015 Q1-Q4
PQ_prop9 = betalike([pars(21)*((y(9,2)+y(10,6)-y(9,6))*5602187-1),(1-pars(21))*((y(9,2)+y(10,6)-y(9,6))*5602187-1)],827285./((y(9,2)+y(10,6)-y(9,6))*5602187));
PQ_prop10 = betalike([pars(21)*((y(10,2)+y(11,6)-y(10,6))*5602187-1),(1-pars(21))*((y(10,2)+y(11,6)-y(10,6))*5602187-1)],852025./((y(10,2)+y(11,6)-y(10,6))*5602187));
PQ_prop11 = betalike([pars(21)*((y(11,2)+y(12,6)-y(11,6))*5602187-1),(1-pars(21))*((y(11,2)+y(12,6)-y(11,6))*5602187-1)],855983./((y(11,2)+y(12,6)-y(11,6))*5602187));
PQ_prop12 = betalike([pars(21)*((y(12,2)+y(13,6)-y(12,6))*5602187-1),(1-pars(21))*((y(12,2)+y(13,6)-y(12,6))*5602187-1)],845098./((y(12,2)+y(13,6)-y(12,6))*5602187));
%2016 Q1-Q4
PQ_prop13 = betalike([pars(21)*((y(13,2)+y(14,6)-y(13,6))*5648259-1),(1-pars(21))*((y(13,2)+y(14,6)-y(13,6))*5648259-1)],832085./((y(13,2)+y(14,6)-y(13,6))*5648259));
PQ_prop14 = betalike([pars(21)*((y(14,2)+y(15,6)-y(14,6))*5648259-1),(1-pars(21))*((y(14,2)+y(15,6)-y(14,6))*5648259-1)],821189./((y(14,2)+y(15,6)-y(14,6))*5648259));
PQ_prop15 = betalike([pars(21)*((y(15,2)+y(16,6)-y(15,6))*5648259-1),(1-pars(21))*((y(15,2)+y(16,6)-y(15,6))*5648259-1)],793453./((y(15,2)+y(16,6)-y(15,6))*5648259));
PQ_prop16 = betalike([pars(21)*((y(16,2)+y(17,6)-y(16,6))*5648259-1),(1-pars(21))*((y(16,2)+y(17,6)-y(16,6))*5648259-1)],775622./((y(16,2)+y(17,6)-y(16,6))*5648259));
%2017 Q1-Q4
PQ_prop17 = betalike([pars(21)*((y(17,2)+y(18,6)-y(17,6))*5702475-1),(1-pars(21))*((y(17,2)+y(18,6)-y(17,6))*5702475-1)],775622./((y(17,2)+y(18,6)-y(17,6))*5702475));
PQ_prop18 = betalike([pars(21)*((y(18,2)+y(19,6)-y(18,6))*5702475-1),(1-pars(21))*((y(18,2)+y(19,6)-y(18,6))*5702475-1)],764726./((y(18,2)+y(19,6)-y(18,6))*5702475));
PQ_prop19 = betalike([pars(21)*((y(19,2)+y(20,6)-y(19,6))*5702475-1),(1-pars(21))*((y(19,2)+y(20,6)-y(19,6))*5702475-1)],739961./((y(19,2)+y(20,6)-y(19,6))*5702475));
PQ_prop20 = betalike([pars(21)*((y(20,2)+y(21,6)-y(20,6))*5702475-1),(1-pars(21))*((y(20,2)+y(21,6)-y(20,6))*5702475-1)],641894./((y(20,2)+y(21,6)-y(20,6))*5702475));
%2018 Q1-Q4
PQ_prop21 = betalike([pars(21)*((y(21,2)+y(22,6)-y(21,6))*5754509-1),(1-pars(21))*((y(21,2)+y(22,6)-y(21,6))*5754509-1)],688451./((y(21,2)+y(22,6)-y(21,6))*5754509));
PQ_prop22 = betalike([pars(21)*((y(22,2)+y(23,6)-y(22,6))*5754509-1),(1-pars(21))*((y(22,2)+y(23,6)-y(22,6))*5754509-1)],683498./((y(22,2)+y(23,6)-y(22,6))*5754509));
PQ_prop23 = betalike([pars(21)*((y(23,2)+y(24,6)-y(23,6))*5754509-1),(1-pars(21))*((y(23,2)+y(24,6)-y(23,6))*5754509-1)],641894./((y(23,2)+y(24,6)-y(23,6))*5754509));
PQ_prop24 = betalike([pars(21)*((y(24,2)+y(25,6)-y(24,6))*5754509-1),(1-pars(21))*((y(24,2)+y(25,6)-y(24,6))*5754509-1)],625054./((y(24,2)+y(25,6)-y(24,6))*5754509));


%Prescription opioid overdoses, 2013-2016: gamma distributed; since muA is
%time-dependent, must add up overdoses each quarter since that's our
%linspace

%Overdoses added together quarters 1-4 for 2013
overdoses_year1 = (pars(17)*1+pars(18))*(y(2,9)-y(1,9))+(pars(17)*2+pars(18))*(y(3,9)-y(2,9))+(pars(17)*3+pars(18))*(y(4,9)-y(3,9))+(pars(17)*4+pars(18))*(y(5,9)-y(4,9));
%Overdoses added together quarters 5-8 for 2014
overdoses_year2 = (pars(17)*5+pars(18))*(y(6,9)-y(5,9))+(pars(17)*6+pars(18))*(y(7,9)-y(6,9))+(pars(17)*7+pars(18))*(y(8,9)-y(7,9))+(pars(17)*8+pars(18))*(y(9,9)-y(8,9));
%Overdoses added together quarters 9-12 for 2015
overdoses_year3 = (pars(17)*9+pars(18))*(y(10,9)-y(9,9))+(pars(17)*10+pars(18))*(y(11,9)-y(10,9))+(pars(17)*11+pars(18))*(y(12,9)-y(11,9))+(pars(17)*12+pars(18))*(y(13,9)-y(12,9));
%Overdoses added together quarters 13-16 for 2016
overdoses_year4 = (pars(17)*13+pars(18))*(y(14,9)-y(13,9))+(pars(17)*14+pars(18))*(y(15,9)-y(14,9))+(pars(17)*15+pars(18))*(y(16,9)-y(15,9))+(pars(17)*16+pars(18))*(y(17,9)-y(16,9));

A_overdose1 = gamlike([(overdoses_year1).^2./(pars(20)^2),(pars(20)^2)./overdoses_year1],351./5519417);
A_overdose2 = gamlike([(overdoses_year2).^2./(pars(20)^2),(pars(20)^2)./overdoses_year2],360./5559702);
A_overdose3 = gamlike([(overdoses_year3).^2./(pars(20)^2),(pars(20)^2)./overdoses_year3],377./5602187);
A_overdose4 = gamlike([(overdoses_year4).^2./(pars(20)^2),(pars(20)^2)./overdoses_year4],381./5648259);


%Heroin/fentanyl overdoses, 2013-2017: gamma distributed 
H_overdose1 = gamlike([(pars(7)*(y(5,10)-y(1,10))).^2./(pars(20)^2),(pars(20)^2)./(pars(7)*(y(5,10)-y(1,10)))],112./5519417);
H_overdose2 = gamlike([(pars(7)*(y(9,10)-y(5,10))).^2./(pars(20)^2),(pars(20)^2)./(pars(7)*(y(9,10)-y(5,10)))],201./5559702);
H_overdose3 = gamlike([(pars(7)*(y(13,10)-y(9,10))).^2./(pars(20)^2),(pars(20)^2)./(pars(7)*(y(13,10)-y(9,10)))],344./5602187);
H_overdose4 = gamlike([(pars(7)*(y(17,10)-y(13,10))).^2./(pars(20)^2),(pars(20)^2)./(pars(7)*(y(17,10)-y(13,10)))],488./5648259);
H_overdose5 = gamlike([(pars(7)*(y(21,10)-y(17,10))).^2./(pars(20)^2),(pars(20)^2)./(pars(7)*(y(21,10)-y(17,10)))],702./5702475);


%Minimizing this because each term is the beta negative log-likelihood function
value = P_prop1+P_prop2+P_prop3+P_prop4+P_prop5+P_prop6+...
        A_prop1+A_prop2+A_prop3+A_prop4+A_prop5+A_prop6+...
        H_prop1+H_prop2+H_prop3+...
        PQ_prop1+PQ_prop2+PQ_prop3+PQ_prop4+PQ_prop5+PQ_prop6+...
        PQ_prop7+PQ_prop8+PQ_prop9+PQ_prop10+PQ_prop11+PQ_prop12+...
        PQ_prop13+PQ_prop14+PQ_prop15+PQ_prop16+PQ_prop17+PQ_prop18+...
        PQ_prop19+PQ_prop20+PQ_prop21+PQ_prop22+PQ_prop23+PQ_prop24+...
        A_overdose1+A_overdose2+A_overdose3+A_overdose4+...
        H_overdose1+H_overdose2+H_overdose3+H_overdose4+H_overdose5;

% Works with original bounds on parameters because proportion correctly
% surveyed stays less than 1 
% value = A_prop1+A_prop2+A_prop3+A_prop4+A_prop5+A_prop6+...
%         H_prop1+H_prop2+H_prop3+...
%         A_overdose1+A_overdose2+A_overdose3+A_overdose4+...
%         H_overdose1+H_overdose2+H_overdose3+H_overdose4+H_overdose5;         

%Want to know how many times getting value = NaN, even though run continues
%and eventually converges because if it gets NAN, it just goes back and
%takes a shorter step 

% if isnan(value)
%      value
% end 

%Want to know if getting value = NaN anywhere in the process for each run;
%if so, will output value of 0
% if isnan(value)
%      value = 0 
% end 

end

           
function alpha = a(t,pars)
if  t<=3.25 
    alpha = pars(1)*t+pars(15);
else 
    alpha = pars(1)*3.25+pars(15)-pars(16)*3.25+pars(16)*t;
end
end

           
function mu_A = muA(t,pars)
    mu_A = pars(17)*t+pars(18);
end


function f = HeroinModel(t,y,pars)
f=zeros(10,1);
f(1)=-a(t,pars)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+muA(t,pars))*y(3)+(pars(6)+pars(7))*y(4);
f(2)=a(t,pars)*y(1)-pars(5)*y(2)-pars(8)*y(2)-pars(9)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(8)*y(2)+(pars(10)*y(5)*y(3))/(y(3)+y(4)+pars(14))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(11)*y(3)-pars(12)*y(3)*y(4)-pars(6)*y(3)-muA(t,pars)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(9)*y(2)*y(4)+pars(12)*y(3)*y(4)+(pars(10)*y(5)*y(4))/(y(3)+y(4)+pars(14))-pars(13)*y(4)-(pars(6)+pars(7))*y(4);
f(5)=pars(11)*y(3)+pars(13)*y(4)-(pars(10)*y(5)*y(3))/(y(3)+y(4)+pars(14))-(pars(10)*y(5)*y(4))/(y(3)+y(4)+pars(14))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time;
% i.e. individuals who enter the P class at any time from S (used in
% Estim1, Estim4) 
f(6) = a(t,pars)*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
% i.e. individuals who enter the A class at any time (used in Estim2)
f(7) = pars(8)*y(2)+(pars(10)*y(5)*y(3))/(y(3)+y(4)+pars(14))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction over time; 
% i.e. individuals who enter the H class at any time (used in Estim3)
f(8) = pars(4)*y(1)*y(4)+pars(9)*y(2)*y(4)+pars(12)*y(3)*y(4)+(pars(10)*y(5)*y(4))/(y(3)+y(4)+pars(14));

%J' ODE to calculate number of prescription opioid addict overdoses over
%time; i.e. individuals who overdose at any time (used in Estim5)
f(9) = muA(t,pars)*y(3);

%K' ODE to calculate number of heroin/fentanayl addict overdoses over
%time; i.e. individuals who overdose at any time (used in Estim6)
f(10) = pars(7)*y(4);
end


