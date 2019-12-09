%File name: heroin_multistart_alphapiecewise_muAlinear_Hweighted.m 

clf;
clear all;

% Realistic parameter bounds
%           [m      betaA     betaP   theta1   epsilon  gamma   theta2   sigma    zeta    theta3    nu       b     P0        A0       H0       R0       c      d         e ]
LowerBounds=[-0.1  0.00001   0.00001   0.05      1       0.005     0.1     0.1    0.0001    10      0.0001   0.1   0.001   0.0001   0.00001  0.00001   -0.1   0.0001   0.001];
UpperBounds=[-0.001  0.01     0.01      0.3      5       0.1       0.4     2       0.2      20       0.2     0.5   0.35      0.01    0.002     0.1    -0.001  0.008     0.1 ];

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

% x vector to estimate, alpha=m*t+b then slope of c, muA=d*t+e
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

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e];


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
 State1=y(:,1);
 State_data_1=[0.815027000000000;0.820092606374491;0.823083462782946;0.824953101284708;0.826266867578232;0.827295589605025;0.828151064299051;0.828888056709765;0.829550987768877;0.830160968116141;0.830731476779983;0.831267853739846;0.831775277792388;0.832369337439226;0.833320273070873;0.834935193709164;0.837008222792739;0.839430514198263;0.841948852372022;0.844480533456716;0.846969875070211;0.849438372705049;0.851927568637527;0.854422077354162;0.856929766092053];
 
 
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
 State_data_2=[0.102000000000000;0.0979710775121495;0.0958934843546083;0.0948328729097304;0.0942413078659304;0.0938618315943667;0.0935931706678741;0.0933881201714721;0.0932086415480486;0.0930379528768021;0.0928659717199934;0.0926897694985381;0.0925055692402903;0.0921990472272650;0.0915016558954936;0.0901080865418251;0.0882267833312419;0.0859689788179024;0.0835905667670046;0.0811767119270997;0.0787855157583665;0.0763980831664409;0.0739761603517100;0.0715351944035062;0.0690707812789624];
 
 
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
 State_data_3=[0.00649000000000000;0.00557005732990452;0.00481680348275666;0.00419948784946296;0.00368993913670148;0.00326506233113871;0.00290653429869580;0.00260036450503723;0.00233437975220614;0.00209967546473880;0.00188914545290367;0.00169753992521592;0.00152187388950920;0.00136034861566794;0.00121144950900786;0.00107448268278180;0.000949202799976213;0.000835191611386870;0.000732334699928782;0.000640322519674607;0.000558716326499436;0.000486843058871295;0.000423549302360444;0.000368831898934075;0.000321520803955404];
 
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
 State_data_4=[0.000483000000000000;0.000534496414891422;0.000594788339622348;0.000666162971638433;0.000751058321370202;0.000852060426969876;0.000971869686648249;0.00111314220199549;0.00127841009146950;0.00146988181803081;0.00168929926139188;0.00193796039917839;0.00221634817354167;0.00252402447348422;0.00285986262535913;0.00322192044711546;0.00360708817907164;0.00401259934329640;0.00443513747456256;0.00487182414548497;0.00531996492795778;0.00577689249206152;0.00623965326272952;0.00670819686641716;0.00717943119012899];
 
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
 State_data_5=[0.0760000000000000;0.0758605363903707;0.0757181401581843;0.0755727825422331;0.0754244617688271;0.0752731721118128;0.0751188937667554;0.0749615876046472;0.0748011915335650;0.0746376442409857;0.0744708526854361;0.0743007137953781;0.0741271137386892;0.0739498985058705;0.0737689344758305;0.0735840366549660;0.0733949495332789;0.0732015805058689;0.0730037037245750;0.0728011210647986;0.0725936294717000;0.0723810121475027;0.0721631109778996;0.0719397688995477;0.0717108477396862];
 
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
 State_data_6=[0;0.0573785967656248;0.114899216796172;0.172458177048749;0.229992810413221;0.287471538007361;0.344879711515249;0.402208900648780;0.459451851027953;0.516603914123558;0.573661332995451;0.630621293386233;0.687474334683457;0.744056080485544;0.799917568498351;0.854488881064813;0.907607869309353;0.959146772567659;1.00913894843627;1.05764053365004;1.10464897495316;1.15017978552844;1.19420650624074;1.23673115265295;1.27773647084120];
 
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
 State_data_7=[0;0.000389180315868358;0.000768681957694587;0.00113696941281152;0.00149220556723511;0.00183233557154486;0.00215520270684619;0.00245876082349312;0.00274112877287119;0.00300082027204611;0.00323685511164484;0.00344873162519481;0.00363679150994321;0.00380215949821286;0.00394618418180325;0.00407068773786603;0.00417783906738057;0.00426967774615484;0.00434844110257091;0.00441614739842455;0.00447465099125710;0.00452571645067030;0.00457095856478861;0.00461034539543194;0.00464540669146992];
 
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
 State_data_8=[0;5.90916012789079e-05;0.000127805695289644;0.000208581158732387;0.000304039989683831;0.000416988963681621;0.000550389660184624;0.000707206150409667;0.000890303171217715;0.00110226663171848;0.00134523863722652;0.00162093940444064;0.00193029708689200;0.00227331589116862;0.00264929484663589;0.00305670553802587;0.00349279604360887;0.00395513629435763;0.00444069629353317;0.00494681083696165;0.00547099664371668;0.00601074864799582;0.00656324201049048;0.00712842187941653;0.00770330889649219];
 
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
 State_data_9=[0;0.000903758167271174;0.00168656213145445;0.00236835354310935;0.00296712411586597;0.00349725715482611;0.00396980839989720;0.00439294638427643;0.00477370847766048;0.00511716207066631;0.00542749017322415;0.00570793333808240;0.00596076204304586;0.00618797828570501;0.00639154819297779;0.00657320612796351;0.00673462451727724;0.00687750964240158;0.00700350347436477;0.00711423674824049;0.00721131150775328;0.00729643289345413;0.00737153156102500;0.00743658890090713;0.00749349282907433];
 
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
 State_data_10=[0;6.44888715589666e-06;1.35999452770244e-05;2.15819615108868e-05;3.05511552712059e-05;4.06949492781865e-05;5.22351061475257e-05;6.54331153632863e-05;8.05713950830950e-05;9.79704256083125e-05;0.000117970028717223;0.000140928717288436;0.000167224509424851;0.000197232966179268;0.000231315541576287;0.000269823554879986;0.000313060837551684;0.000361313026170702;0.000414823354129494;0.000473772864045376;0.000538340835828419;0.000608665043975220;0.000684854794873369;0.000766906286158579;0.000854915255102489];
 
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
 %Data1=[0.331992810413221;0.323700348480663;0.321231125203553;0.312639103866186;0.285267888975045;0.251873011646405];
 
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
 %Data2=[0.00798220556723511;0.00493886234233757;0.00323004248927816;0.00206292144694655;0.00124601472385274;0.000729472026712258];
 
 
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
 %Data3=[0.00133732150290409;0.00231840400714378;0.00377884713025854];
 
 
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
%Data4=[0.159345888525904;0.155369684395518;0.153226232495256;0.152038272755859;0.151303670490909;0.150774641564431;0.150353591755673;0.149995131706456;0.149672168321776;0.149358805511151;0.149050068737102;0.148734280370870;0.148085068864361;0.147016117783976;0.145043027658001;0.142252916511528;0.138990114824242;0.135356694964831;0.131519862311692;0.127622405412348;0.123717988037982;0.119806623029534;0.115881867009583;0.111940266417211];
 
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
 %Data5=[6.93893503425710e-05;6.66363564939064e-05;6.36911749420292e-05;6.03701269735515e-05;5.65022143567912e-05];
 
 
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
 %Data6=[2.99015652964405e-05;4.35881487901691e-05;6.30645587490766e-05;9.02990840950082e-05;0.000126571634010397];
 
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

fprintf('muA values')

disp(muA(0,pars))
disp(muA(1,pars))
disp(muA(2,pars))
disp(muA(3,pars))
disp(muA(4,pars))
disp(muA(5,pars))
disp(muA(6,pars))



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

% Parameter vector
pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e];

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
 %Data1=[0.331574325986979;0.322519449361526;0.319556848622964;0.310384043536721;0.283643697595965;0.250408236361685];
 
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
 %Data2=[0.00802951524269548;0.00770633035801977;0.00735449191012147;0.00695263159030065;0.00648268985580833;0.00594322842197851];
 
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
 %Data3=[0.00108078107424111;0.00155735427751202;0.00221294953998018];
 

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
 %Data4=[0.159345888525904;0.155369684395518;0.153226232495256;0.152038272755859;0.151303670490909;0.150774641564431;0.150353591755673;0.149995131706456;0.149672168321776;0.149358805511151;0.149050068737102;0.148734280370870;0.148085068864361;0.147016117783976;0.145043027658001;0.142252916511528;0.138990114824242;0.135356694964831;0.131519862311692;0.127622405412348;0.123717988037982;0.119806623029534;0.115881867009583;0.111940266417211];
 

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
 %Data5=[6.93893503425710e-05;6.66363564939064e-05;6.36911749420292e-05;6.03701269735515e-05;5.65022143567912e-05];
 
 
 %The difference between estimated value and data
 Diff5=Estim5-Data5;
 
 
 % Yearly output from the model as a proportion of the population that
 % overdoses from H class at some point during the year for
 % years 2013-2017, Estim6 is a column vector
 Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10); y(17,10)-y(13,10); y(21,10)-y(17,10)];
    
 % Actual proportion of heroin addicts that overdose each year 2013-2017
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 
 %Testing Data
 %Data6=[2.99015652964405e-05;4.35881487901691e-05;6.30645587490766e-05;9.02990840950082e-05;0.000126571634010397];
 
 

 %The difference between estimated value and data 
 Diff6=Estim6-Data6;

 %%% For testing purposes: states and corresponding simulated data 
 
 % Comparing simulated data for susceptibles to output of model for
 % susceptibles 
 State1=y(:,1);
 State_data_1=[0.815027000000000;0.819243442292438;0.821530683367519;0.822806182288711;0.823598962135632;0.824173041976803;0.824631376776187;0.825021102687401;0.825368945864479;0.825702832886299;0.826027387994084;0.826348807646423;0.826670523172988;0.827204385975363;0.828056538080662;0.829475221404466;0.831376706313333;0.833598255818893;0.836012519622903;0.838509164267455;0.841016034961814;0.843522233544351;0.846033461705639;0.848551996583345;0.851081705812981];
 
 State_diff_1= State1-State_data_1;
 
 % Comparing simulated data for prescription users to output of model for
 % prescription users
 State2=y(:,2);
 State_data_2=[0.102000000000000;0.0979350465762521;0.0957985200591501;0.0946721855501568;0.0940264994589885;0.0935970720286697;0.0932806571536079;0.0930298569736657;0.0928176876851886;0.0926159633432694;0.0924198005369926;0.0922227104376727;0.0920209693481170;0.0916025808567509;0.0908612882369540;0.0895492181874396;0.0877505212491147;0.0856282651335113;0.0833101794659122;0.0809069353177249;0.0784909818795858;0.0760735486454169;0.0736491476490123;0.0712158118381970;0.0687700013978874];
 
 State_diff_2=State2-State_data_2;
 
 
 % Comparing simulated data for opioid addicts to output of model for
 % opioid addicts
 State3=y(:,3);
 State_data_3=[0.00649000000000000;0.00642840137268246;0.00636638767597314;0.00630407914981803;0.00624141527973456;0.00617820957990189;0.00611420802401626;0.00604913926917130;0.00598272984936470;0.00591470799628756;0.00584480387567860;0.00577275835700289;0.00569831384038206;0.00562121141462316;0.00554118661222267;0.00545795824310944;0.00537122324430698;0.00528093115653817;0.00518698675766366;0.00508937223818101;0.00498815906110979;0.00488350882916058;0.00477558419516409;0.00466457361702389;0.00455066846033045];
 
 State_diff_3=State3-State_data_3;
 
 % Comparing simulated data for heroin addicts to output of model for
 % heroin addicts
 State4=y(:,4);
 State_data_4=[0.000483000000000000;0.000532573368261227;0.000586268739188888;0.000644770469111269;0.000708661356861379;0.000778504302881531;0.000854864279527708;0.000938313465224011;0.00102944506752207;0.00112885153329098;0.00123715490797918;0.00135500976375805;0.00148307990015417;0.00162192324774424;0.00177205259461763;0.00193356551021345;0.00210659966010781;0.00229096738537788;0.00248661042924719;0.00269340711228656;0.00291119462634154;0.00313969683403486;0.00337869547256479;0.00362784906197863;0.00388677658914036];
 
 State_diff_4=State4-State_data_4;
 
 
 % Comparing simulated data for stably recovered individuals to output of model for
 % stably recovered individuals 
 State5=y(:,5);
 State_data_5=[0.0760000000000000;0.0758605363903707;0.0757181401581843;0.0755727825422331;0.0754244617688271;0.0752731721118128;0.0751188937667554;0.0749615876046472;0.0748011915335650;0.0746376442409857;0.0744708526854361;0.0743007137953781;0.0741271137386892;0.0739498985058705;0.0737689344758305;0.0735840366549660;0.0733949495332789;0.0732015805058689;0.0730037037245750;0.0728011210647986;0.0725936294717000;0.0723810121475027;0.0721631109778996;0.0719397688995477;0.0717108477396862];
 
 State_diff_5=State5-State_data_5;
 
 % Comparing simulated data for proportion of individuals entering P throughout the year and the model output 
 State6=y(:,6);
 State_data_6=[0;0.0573458885259043;0.114780526345170;0.172208238781276;0.229574325986979;0.286851497018899;0.344029066554661;0.401102001156726;0.458067275889516;0.514921756526104;0.571664598693985;0.628294866894094;0.684806436827292;0.740870536343535;0.796284073270760;0.850465812691807;0.903169511015896;0.954409104591023;1.00413753442234;1.05234721726812;1.09906268736275;1.14428969352114;1.18802276790526;1.23025548726583;1.27097994184484];
 
 State_diff_6=State6-State_data_6;
 
 % Comparing simulated data for proportion of individuals entering A throughout the year and the model output 
 State7=y(:,7);
 State_data_7=[0;0.000391415611247909;0.000778444323641322;0.00116116519535126;0.00153951524269548;0.00191330449097840;0.00228226239100783;0.00264608182861675;0.00300443032098068;0.00335697764837222;0.00370336294517546;0.00404322486445433;0.00437619238173745;0.00470184783457662;0.00501978965584632;0.00532953818779343;0.00563051013165605;0.00592254778008270;0.00620529466859068;0.00647849730674416;0.00674197674315740;0.00699561400301962;0.00723935598383820;0.00747316986089003;0.00769704610402612];
 
 State_diff_7=State7-State_data_7;
 
 % Comparing simulated data for proportion of individuals entering H throughout the year and the model output 
 State8=y(:,8);
 State_data_8=[0;5.71593159343135e-05;0.000119203854875563;0.000186890728781261;0.000260877969620559;0.000341813799215656;0.000430357164553689;0.000527182816573776;0.000632997687000288;0.000748508096864143;0.000874466227767260;0.00101166364903441;0.00116090689699023;0.00132292027690361;0.00149836952118418;0.00168752939122077;0.00189077653681624;0.00210799552038908;0.00233932289751394;0.00258480317795645;0.00284444143224080;0.00311814577369257;0.00340584049713199;0.00370733043699128;0.00402236744260534];
 
 State_diff_8=State8-State_data_8;
 
 % Comparing simulated data for proportion of individuals that overdose out of A throughout the year and the model output 
 State9=y(:,9);
 State_data_9=[0;1.75995816570893e-05;3.50324406263730e-05;5.22960416938377e-05;6.93893503425710e-05;8.63111795585811e-05;0.000103059748660548;0.000119632455033441;0.000136025706836477;0.000152235941259567;0.000168258219724049;0.000184087152748587;0.000199716881778507;0.000215139755304727;0.000230349023780927;0.000245335633100410;0.000260087008752058;0.000274598662896371;0.000288859484720705;0.000302859700642183;0.000316589223108849;0.000330037112269178;0.000343195201984157;0.000356055262949137;0.000368609952484885];
 
 State_diff_9=State9-State_data_9;
 
 % Comparing simulated data for proportion of individuals that overdose out of H the year and the model output 
 State10=y(:,10);
 State_data_10=[0;6.44104278922220e-06;1.35301173192590e-05;2.13289986160856e-05;2.99015652964405e-05;3.93202616591573e-05;4.96648792677383e-05;6.10225099382901e-05;7.34897140866096e-05;8.71631794807317e-05;0.000102153401548173;0.000118576701224108;0.000136554272835686;0.000156227381068224;0.000177725887206501;0.000201200240605014;0.000226853356930694;0.000254746557559656;0.000285045252898473;0.000317890550131164;0.000353424990941092;0.000391805018096980;0.000433151673792301;0.000477589189996024;0.000525230376985628];
 
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
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 
 %Put more weight on heroin data to see if it STILL gives such as drastic
 %projected increase in heroin use in future 
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+2*norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);
 

 
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
    alpha = pars(1)*t+pars(15);
else 
    alpha = pars(1)*3.25+pars(15)-pars(16)*3.25+pars(16)*t;
    %alpha = pars(17)*t+pars(18);
    %alpha = pars(1)*t+pars(16);
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


