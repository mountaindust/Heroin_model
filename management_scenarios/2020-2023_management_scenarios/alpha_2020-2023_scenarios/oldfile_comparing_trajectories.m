%File name: heroin_multistart_alphapiecewise_muAlinear_testing.m
clf;
clear all;

m=-0.00559565027929907;
beta_A=0.000878431642350708;
beta_P=6.54313717400116e-05;
theta_1=0.222457489109919;
epsilon=2.52786996559537;
mu=0.00710; 
mu_H=0.0466; 
gamma=0.00505079477762453;
theta_2=0.236479520411597;
sigma=0.101518004918260;
zeta=0.198182427387906;
theta_3=19.7264083013258;
nu=0.000531263148928530;
omega=0.0000000001;
b=0.270110337915851;
c=-0.0269690987063522;
d=0.000977482526657751;
e=0.00883138792481281;

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e];

% Final time and last entry of tspan is # of equally spaced points from 0 to N (quarterly linspace)
N = 3;
tspan=linspace(0,N,13);
%For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);


% Initial Conditions
P0=0.0585;
A0=0.0037;
H0=0.00597;
R0=0.00751;
S0=1-P0-A0-H0-R0;
X0=1.443875288499938;
L0=0.006431117471779;
M0=0.006361097357661;
J0=4.568512848210284e-04;
K0=7.308319929066670e-04;
initials = [S0;P0;A0;H0;R0;X0;L0;M0;J0;K0];

[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

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

format short
fprintf('alpha values')
disp(a(0,pars))
disp(a(1,pars))
disp(a(2,pars))
disp(a(3,pars))
% disp(a(3.25,pars))
% disp(a(4,pars))
% disp(a(5,pars))
% disp(a(6,pars))
% disp(a(7,pars))
% disp(a(8,pars))
% disp(a(9,pars))
% disp(a(10,pars))

format long
fprintf('muA values')
disp(muA(0,pars))
disp(muA(1,pars))
disp(muA(2,pars))
disp(muA(3,pars))
% disp(muA(3.25,pars))
% disp(muA(4,pars))
% disp(muA(5,pars))
% disp(muA(6,pars))
% disp(muA(7,pars))
% disp(muA(8,pars))
% disp(muA(9,pars))
% disp(muA(10,pars))
 
 A_baseline_alpha_case=[0.00370000000000000;0.00353278337029652;0.00336289307560543;0.00319097482980213;0.00301766402796066;0.00284378632891392;0.00267026100551299;0.00249770644297795;0.00232715473139335;0.00215947534975907;0.00199553547623715;0.00183611325313993;0.00168192319325566];
 H_baseline_alpha_case=[0.00597000000000000;0.00646253334841020;0.00698258162521298;0.00753001934779802;0.00810629377051649;0.00871221154764143;0.00934807558313681;0.0100140013779746;0.0107101117664536;0.0114363784036355;0.0121928145645927;0.0129794511679071;0.0137963409589377];
 
 A_25percent_alpha_case=[0.00370000000000000;0.00353267040767489;0.00336238928474401;0.00318968332757997;0.00301514719584218;0.00283955422999254;0.00266379884441428;0.00248877771301105;0.00231554824945441;0.00214487141341584;0.00197761384549922;0.00181450857447343;0.00165638654016796];
 H_25percent_alpha_case=[0.00597000000000000;0.00646252562871979;0.00698213581555402;0.00752931844143701;0.00810535995274669;0.00871100819807622;0.00934652875810367;0.0100120066686606;0.0107074079365409;0.0114327074419068;0.0121878730733820;0.0129728826579485;0.0137877434894481];
 
 A_50percent_alpha_case=[0.00370000000000000;0.00353280392347552;0.00336233468380589;0.00318892577674831;0.00301315616123394;0.00283579021459923;0.00265775076305633;0.00248000877817790;0.00230356335342014;0.00212938188714762;0.00195838905036281;0.00179145596662234;0.00162939725352773];
 H_50percent_alpha_case=[0.00597000000000000;0.00646073977077021;0.00697854297962477;0.00752471520320035;0.00810023226051878;0.00870540350069781;0.00934038930763936;0.0100051708090444;0.0106996910829246;0.0114238593390079;0.0121775733566048;0.0129607303266466;0.0137732386906761];
 
  
 %compare all 3 alpha cases
 figure(1)
 hold all
 connect1=interp1(t,A_baseline_alpha_case,t);
 plot(t,connect1,'k-','LineWidth',3)
 connect2=interp1(t,A_25percent_alpha_case,t);
 plot(t,connect2,'g-','LineWidth',3)
 connect3=interp1(t,A_50percent_alpha_case,t);
 plot(t,connect3,'r-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Opioid addicts')
 set(gca, 'xtick', [ 0 1 2 3 ])
 set(gca,'xticklabel',{'2020', '2021', '2022', '2023'})   
                   
 figure(2)
 hold all
 plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 set(gca, 'xtick', [ 0 1 2 3 ])
 set(gca,'xticklabel',{'2020', '2021', '2022', '2023'})         
 
      
 
 