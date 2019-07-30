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
numstartpoints=100;

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
mu_A=0.0109;   
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
 State_data_1=[0.837700000000000;0.838534996478993;0.839061867855981;0.839405698053262;0.839657227344269;0.839873266893214;0.840064498265746;0.840248630799731;0.840424390178778;0.840600913830506;0.840776725221554;0.840953449561807;0.841267196084226;0.841682851304061;0.842354764682153;0.843722511684880;0.845580974731544;0.847734412745936;0.849976013563451;0.852296427867473;0.854604011778378;0.856945453031114;0.859307761250540;0.861705891610171;0.864136103015468];
 
 
 % Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1),'k-','LineWidth',3)
 plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State2=y(:,2);
 State_data_2=[0.0950000000000000;0.0941897164482818;0.0936883992203288;0.0933707168788075;0.0931457543702420;0.0929566028758150;0.0927925667344050;0.0926359241367061;0.0924879545363399;0.0923395288020814;0.0921921343157709;0.0920441580105594;0.0917598202035540;0.0913741484475699;0.0907329583133100;0.0893976273587691;0.0875741695414704;0.0854587639471082;0.0832587126759282;0.0809833952317272;0.0787246900031921;0.0764358443142316;0.0741299173005789;0.0717919639776442;0.0694257743772307];
 
 
 % Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2),'b-','LineWidth',3)
 plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State3=y(:,3);
 State_data_3=[0.00650000000000000;0.00640147360863710;0.00630841892464260;0.00622080504056483;0.00613837648161604;0.00606061358468768;0.00598730250468231;0.00591782054101086;0.00585211810398835;0.00578955787287069;0.00572998876009984;0.00567307574375391;0.00561794611083936;0.00556486117130198;0.00551356890987475;0.00546271434321223;0.00541139508652284;0.00535920227942272;0.00530567313828903;0.00525081343727614;0.00519442638130915;0.00513663227446000;0.00507740917024204;0.00501679786056049;0.00495479351211911];
 
 % Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3),'r-','LineWidth',3)
 plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 State4=y(:,4);
 State_data_4=[0.000800000000000000;0.000823896466700484;0.000848818115045734;0.000874822053842100;0.000901945989947395;0.000930222798685114;0.000959663173388433;0.000990296149434153;0.00102212444885004;0.00105518015190150;0.00108947220227729;0.00112501867916510;0.00116185177990351;0.00119996122665827;0.00123934780031855;0.00128001182997199;0.00132193592415434;0.00136512744171482;0.00140959377968833;0.00145537027214857;0.00150247386593604;0.00155092997851498;0.00160075862147620;0.00165198059483702;0.00170461707270081];
 
 % Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
 plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 

 State5=y(:,5);
 State_data_5=[0.0600000000000000;0.0600499169974244;0.0600924958839867;0.0601279579732353;0.0601566958133574;0.0601792938468531;0.0601959693209655;0.0602073283722037;0.0602134127310209;0.0602148193415499;0.0602116794991715;0.0602042980035911;0.0601931858203510;0.0601781778492830;0.0601593602932229;0.0601371347820321;0.0601115247153813;0.0600824935851436;0.0600500068420274;0.0600139931900415;0.0599743979683702;0.0599311403972087;0.0598841536511476;0.0598333659497914;0.0597787120151490];
 
 % Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5),'Color', [0.7,0,0.7],'LineWidth',3)
 plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State6=y(:,6);
 State_data_6=[0;0.0586137833088671;0.117169747574732;0.175652646158159;0.234052252100597;0.292362567480044;0.350582516052390;0.408710332644865;0.466746137653744;0.524688971873682;0.582538974881024;0.640295963083095;0.697700548583259;0.754819323680538;0.811432226441668;0.866718158876736;0.920517612670827;0.972849693891269;1.02372316861514;1.07313756831715;1.12109293990422;1.16757721508125;1.21258231124827;1.25609718151729;1.29811272855400];
 
 % Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6),'LineWidth',3)
 plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend({'Proportion that enter P at some point during the year'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State7=y(:,7);
 State_data_7=[0;0.000451718570797464;0.000901231205390013;0.00134868403100729;0.00179417278989685;0.00223771264669495;0.00267928901873512;0.00311883151675718;0.00355632735593305;0.00399167813082840;0.00442484717718664;0.00485577092873549;0.00528401148866773;0.00570968832569186;0.00613262069844567;0.00655177926242629;0.00696624637342801;0.00737553410451814;0.00777910002703949;0.00817682701218087;0.00856841517230218;0.00895384981114964;0.00933298941441715;0.00970574901776474;0.0100720034655782];
 
 % Simulated data for L and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7),'LineWidth',3)
 plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend({'Proportion that enter A at some point during the year'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 State8=y(:,8);
 State_data_8=[0;4.57308948932932e-05;9.31174285469629e-05;0.000142260648144533;0.000193235018789079;0.000246116424132701;0.000300931412810888;0.000357758821776690;0.000416606513505489;0.000477562542178301;0.000540653827889325;0.000605932148797455;0.000673508784665315;0.000743349406411762;0.000815470056475564;0.000889947120919525;0.000966795349471101;0.00104605328007820;0.00112776159574155;0.00121199127247183;0.00129879365911876;0.00138823118217508;0.00148036035237268;0.00157523903848484;0.00167292619086782];
 
 % Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8),'LineWidth',3)
 plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend({'Proportion that enter H at some point during the year'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) %correspond to the actual t values from t vector that I want to label
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})

 
 State9=y(:,9);
 State_data_9=[0;1.74493826209240e-05;3.45609275579770e-05;5.13493724982908e-05;6.78432439771329e-05;8.40394435967946e-05;9.99686807193538e-05;0.000115623556443435;0.000131034858309164;0.000146194877393502;0.000161128093589656;0.000175828073271007;0.000190310509006145;0.000204580391685494;0.000218642971067753;0.000232500964757226;0.000246157146311893;0.000259612975678528;0.000272868830134806;0.000285924740705594;0.000298780340965029;0.000311434922007524;0.000323887561397308;0.000336137093840454;0.000348182168683961];
 
 % Simulated data for J and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 plot(t,y(:,9),'LineWidth',3)
 plot(t(1:end), State_data_9, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('J(t)')
 legend({'Proportion that overdose from A'}, 'FontSize', 11)%,'data simulated' )
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) %correspond to the actual t values from t vector that I want to label
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})


 State10=y(:,10);
 State_data_10=[0;1.08932978164076e-05;2.22023690014102e-05;3.39519496113532e-05;4.61852744094627e-05;5.88976341735230e-05;7.21444158389913e-05;8.59128270882070e-05;0.000100270741817981;0.000115202868627527;0.000130778609500791;0.000146982322853206;0.000163874870208086;0.000181471234905930;0.000199802083985234;0.000218902389452747;0.000238804059254965;0.000259537459076326;0.000281133046239279;0.000303621626248499;0.000327034332785866;0.000351402469833880;0.000376757380842180;0.000403130448089299;0.000430553021650780];
 
 % Simulated data for K and corresponding ODE solution plotted on top 
 figure(10)
 hold all
 plot(t,y(:,10),'LineWidth',3)
 plot(t(1:end), State_data_10, 'x')
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
 %Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
 
 %Testing Data 
 Data1=[0.329052252100597;0.325839639923389;0.323442365465854;0.314576884291122;0.288149496774865;0.255744478652972];
 
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
 %Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
 
 %Testing Data
 Data2=[0.00829417278989685;0.00790053104765223;0.00757980223672303;0.00730018099559964;0.00701356388539701;0.00669801467458518];
 
 
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
 %Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 
 %Testing Data
 Data3=[0.00112531748466380;0.00127902672000986;0.00145513834470930];
 
 
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
 %Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
  %      833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
   %     827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
    %    832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
     %   775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
      %  688502./5754509; 683722./5754509; 641942./5754509; 625162./5754509];
 
 %Testing Data
 Data4=[0.153613783308867;0.152745680714147;0.152171297803756;0.151770322821245;0.151456069749689;0.151176551448161;0.150920383326880;0.150671729145585;0.150430788756278;0.150189531809423;0.149949122517842;0.149448743510723;0.148878595300833;0.147987051208701;0.146018890748377;0.143197081152860;0.139906250761913;0.136332238670983;0.132673112377937;0.128938766818795;0.125208965180218;0.121440940481250;0.117644787569601;0.113807511014357];
 
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
 
 Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9); y(21,9)-y(17,9)];
     
 % Actual Data for years 2013-2017
 %Data5=[348./5519417; 381./5559702; 463./5602187; 551./5648259; 381./5702475];
 
 %Testing Data
 Data5=[6.78432439771329e-05;6.31916143320307e-05;5.92756506969818e-05;5.58466373057473e-05;5.26231946531359e-05];
 
 
  % Data points from proportion that is in A at some point and overdoses in the year and corresponding ODE solution points 
 figure(15)
 hold all
 z5 = linspace(0,4,5); %defines mesh where going to plot Estim5, Data5 values 
 scatter(z5, Estim5, 100,'o');
 scatter(z5, Data5, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})                   
                   
 
 % Yearly simulation of individuals overdosing from H class during the year for years
% 2013-2017 (anyone in H class throughout the year times mu_H)
 
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
        y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual Data for years 2013-2017
 %Data6=[116./5519417; 216./5559702; 374./5602187; 554./5648259; 811./5702475];
 
 %Testing Data
 Data6=[4.61852744094627e-05;5.40854674085182e-05;6.36041283901053e-05;7.49291890468790e-05;8.82302735309003e-05];
 
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
mu=0.00868;  
mu_A=0.0109;   
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
 
 %Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];

 %Testing Data
 Data1=[0.329052252100597;0.325839639923389;0.323442365465854;0.314576884291122;0.288149496774865;0.255744478652972];
 
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
 
 %Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
 
 %Testing Data
 Data2=[0.00829417278989685;0.00790053104765223;0.00757980223672303;0.00730018099559964;0.00701356388539701;0.00669801467458518];
 
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
 
 
 %Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 
 %Testing Data
 Data3=[0.00112531748466380;0.00127902672000986;0.00145513834470930];
 

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
 %Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
  %      833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
   %     827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
    %    832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
     %   775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
      %  688502./5754509; 683722./5754509; 641942./5754509; 625162./5754509];
      
 %Testing Data
 Data4=[0.153613783308867;0.152745680714147;0.152171297803756;0.151770322821245;0.151456069749689;0.151176551448161;0.150920383326880;0.150671729145585;0.150430788756278;0.150189531809423;0.149949122517842;0.149448743510723;0.148878595300833;0.147987051208701;0.146018890748377;0.143197081152860;0.139906250761913;0.136332238670983;0.132673112377937;0.128938766818795;0.125208965180218;0.121440940481250;0.117644787569601;0.113807511014357];
 

 %The difference between estimated value and data 
 Diff4=Estim4-Data4;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from A class at some point during the year for
 % years 2013-2017, Estim5 is a column vector
 Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9); y(21,9)-y(17,9)];
    
 % Actual proportion of addicted prescription opioid users that overdose
 % each year 2013-2017
 
 %Data5=[348./5519417; 381./5559702; 463./5602187; 551./5648259; 381./5702475];
 
 %Testing Data
 Data5=[6.78432439771329e-05;6.31916143320307e-05;5.92756506969818e-05;5.58466373057473e-05;5.26231946531359e-05];

 
 %The difference between estimated value and data
 Diff5=Estim5-Data5;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from H class at some point during the year for
 % years 2013-2017, Estim6 is a column vector
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10); y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual proportion of heroin addicts that overdose each year 2013-2017
 %Data6=[116./5519417; 216./5559702; 374./5602187; 554./5648259; 811./5702475];
 
 %Testing Data
 Data6=[4.61852744094627e-05;5.40854674085182e-05;6.36041283901053e-05;7.49291890468790e-05;8.82302735309003e-05];
 
 

 %The difference between estimated value and data 
 Diff6=Estim6-Data6;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.837700000000000;0.838534996478993;0.839061867855981;0.839405698053262;0.839657227344269;0.839873266893214;0.840064498265746;0.840248630799731;0.840424390178778;0.840600913830506;0.840776725221554;0.840953449561807;0.841267196084226;0.841682851304061;0.842354764682153;0.843722511684880;0.845580974731544;0.847734412745936;0.849976013563451;0.852296427867473;0.854604011778378;0.856945453031114;0.859307761250540;0.861705891610171;0.864136103015468];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.0950000000000000;0.0941897164482818;0.0936883992203288;0.0933707168788075;0.0931457543702420;0.0929566028758150;0.0927925667344050;0.0926359241367061;0.0924879545363399;0.0923395288020814;0.0921921343157709;0.0920441580105594;0.0917598202035540;0.0913741484475699;0.0907329583133100;0.0893976273587691;0.0875741695414704;0.0854587639471082;0.0832587126759282;0.0809833952317272;0.0787246900031921;0.0764358443142316;0.0741299173005789;0.0717919639776442;0.0694257743772307];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00650000000000000;0.00640147360863710;0.00630841892464260;0.00622080504056483;0.00613837648161604;0.00606061358468768;0.00598730250468231;0.00591782054101086;0.00585211810398835;0.00578955787287069;0.00572998876009984;0.00567307574375391;0.00561794611083936;0.00556486117130198;0.00551356890987475;0.00546271434321223;0.00541139508652284;0.00535920227942272;0.00530567313828903;0.00525081343727614;0.00519442638130915;0.00513663227446000;0.00507740917024204;0.00501679786056049;0.00495479351211911];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.000800000000000000;0.000823896466700484;0.000848818115045734;0.000874822053842100;0.000901945989947395;0.000930222798685114;0.000959663173388433;0.000990296149434153;0.00102212444885004;0.00105518015190150;0.00108947220227729;0.00112501867916510;0.00116185177990351;0.00119996122665827;0.00123934780031855;0.00128001182997199;0.00132193592415434;0.00136512744171482;0.00140959377968833;0.00145537027214857;0.00150247386593604;0.00155092997851498;0.00160075862147620;0.00165198059483702;0.00170461707270081];
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.0600000000000000;0.0600499169974244;0.0600924958839867;0.0601279579732353;0.0601566958133574;0.0601792938468531;0.0601959693209655;0.0602073283722037;0.0602134127310209;0.0602148193415499;0.0602116794991715;0.0602042980035911;0.0601931858203510;0.0601781778492830;0.0601593602932229;0.0601371347820321;0.0601115247153813;0.0600824935851436;0.0600500068420274;0.0600139931900415;0.0599743979683702;0.0599311403972087;0.0598841536511476;0.0598333659497914;0.0597787120151490];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0586137833088671;0.117169747574732;0.175652646158159;0.234052252100597;0.292362567480044;0.350582516052390;0.408710332644865;0.466746137653744;0.524688971873682;0.582538974881024;0.640295963083095;0.697700548583259;0.754819323680538;0.811432226441668;0.866718158876736;0.920517612670827;0.972849693891269;1.02372316861514;1.07313756831715;1.12109293990422;1.16757721508125;1.21258231124827;1.25609718151729;1.29811272855400];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000451718570797464;0.000901231205390013;0.00134868403100729;0.00179417278989685;0.00223771264669495;0.00267928901873512;0.00311883151675718;0.00355632735593305;0.00399167813082840;0.00442484717718664;0.00485577092873549;0.00528401148866773;0.00570968832569186;0.00613262069844567;0.00655177926242629;0.00696624637342801;0.00737553410451814;0.00777910002703949;0.00817682701218087;0.00856841517230218;0.00895384981114964;0.00933298941441715;0.00970574901776474;0.0100720034655782];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;4.57308948932932e-05;9.31174285469629e-05;0.000142260648144533;0.000193235018789079;0.000246116424132701;0.000300931412810888;0.000357758821776690;0.000416606513505489;0.000477562542178301;0.000540653827889325;0.000605932148797455;0.000673508784665315;0.000743349406411762;0.000815470056475564;0.000889947120919525;0.000966795349471101;0.00104605328007820;0.00112776159574155;0.00121199127247183;0.00129879365911876;0.00138823118217508;0.00148036035237268;0.00157523903848484;0.00167292619086782];
 
 State_diff_8=State8-State_data_8;
 
 % Comparing simulated data for proportion of individuals that overdose out of A throughout the year and the model output 
 State9=y(:,9);
 State_data_9=[0;1.74493826209240e-05;3.45609275579770e-05;5.13493724982908e-05;6.78432439771329e-05;8.40394435967946e-05;9.99686807193538e-05;0.000115623556443435;0.000131034858309164;0.000146194877393502;0.000161128093589656;0.000175828073271007;0.000190310509006145;0.000204580391685494;0.000218642971067753;0.000232500964757226;0.000246157146311893;0.000259612975678528;0.000272868830134806;0.000285924740705594;0.000298780340965029;0.000311434922007524;0.000323887561397308;0.000336137093840454;0.000348182168683961];
 
 State_diff_9=State9-State_data_9;
 
 % Comparing simulated data for proportion of individuals that overdose out of H the year and the model output 
 State10=y(:,10);
 State_data_10=[0;1.08932978164076e-05;2.22023690014102e-05;3.39519496113532e-05;4.61852744094627e-05;5.88976341735230e-05;7.21444158389913e-05;8.59128270882070e-05;0.000100270741817981;0.000115202868627527;0.000130778609500791;0.000146982322853206;0.000163874870208086;0.000181471234905930;0.000199802083985234;0.000218902389452747;0.000238804059254965;0.000259537459076326;0.000281133046239279;0.000303621626248499;0.000327034332785866;0.000351402469833880;0.000376757380842180;0.000403130448089299;0.000430553021650780];
 
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
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4);

 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5);

 %value=norm(Diff5,2)./norm(Data5);
 
 % For testing purposes with states and data sets
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8)+norm(State_diff_9,2)./norm(State_data_9)+norm(State_diff_10,2)./norm(State_data_10);
  
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

%J' ODE to calculate number of prescription opioid addict overdoses over
%time; i.e. individuals who overdose at any time (used in Estim5)
f(9) = pars(7)*y(3);

%K' ODE to calculate number of heroin/fentanayl addict overdoses over
%time; i.e. individuals who overdose at any time (used in Estim6)
f(10) = pars(8)*y(4);
end

