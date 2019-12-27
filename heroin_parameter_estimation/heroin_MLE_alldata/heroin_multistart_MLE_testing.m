%File name: heroin_multistart_alphapiecewise_muAlinear_testing.m
clf;
clear all;

m=-0.012143181981601;
beta_A=0.000020764377750;
beta_P=0.000010811238536;
theta_1=0.050025981088221;
epsilon=2.537290448818093;
mu=0.00710; 
mu_H=0.0466; 
gamma=0.005000757797158;
theta_2=0.306352161471387;
sigma=1.999359146133607;
zeta=0.000113666496202;
theta_3=10.972422802831584;
nu=0.000142212638234;
omega=0.0000000001;
b=0.499999876215983;
c=-0.055790243975095;
d=0.007998591530345;
e=0.099987210261792;
p=0.549159307571905;
var_1=0.000141332963870;
p2=0.589828041076585;

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e,p,var_1,p2];

% Final time and last entry of tspan is # of equally spaced points from 0 to N (quarterly linspace)
N = 6;
%tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
tspan=linspace(0,N,3000);


% Initial Conditions
P0= 0.157602004252730;
A0=0.008521073242881;
H0=0.000933614135117;
R0=0.006482116886049;
S0=1-P0-A0-H0-R0;
X0=0;
L0=0;
M0=0;
J0=0;
K0=0;
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
disp(a(3.25,pars))
disp(a(4,pars))
disp(a(5,pars))
disp(a(6,pars))

format long
fprintf('muA values')
disp(muA(0,pars))
disp(muA(1,pars))
disp(muA(2,pars))
disp(muA(3,pars))
disp(muA(3.25,pars))
disp(muA(4,pars))
disp(muA(5,pars))
disp(muA(6,pars))

% Comment out if don't need objective function value 

Estim1=[pars(21)*(y(1,2)+y(5,6)-y(1,6)); pars(21)*(y(5,2)+y(9,6)-y(5,6)); pars(21)*(y(9,2)+y(13,6)-y(9,6));...
        pars(21)*(y(13,2)+y(17,6)-y(13,6)); pars(21)*(y(17,2)+y(21,6)-y(17,6)); pars(21)*(y(21,2)+y(25,6)-y(21,6))];
     
Estim2=[pars(19)*(y(1,3)+y(5,7)-y(1,7)); pars(19)*(y(5,3)+y(9,7)-y(5,7)); pars(19)*(y(9,3)+y(13,7)-y(9,7));...
        pars(19)*(y(13,3)+y(17,7)-y(13,7)); pars(19)*(y(17,3)+y(21,7)-y(17,7)); pars(19)*(y(21,3)+y(25,7)-y(21,7))];

Estim3=[pars(19)*(y(5,4)+y(9,8)-y(5,8)); pars(19)*(y(9,4)+y(13,8)-y(9,8)); pars(19)*(y(13,4)+y(17,8)-y(13,8))];

Estim4=pars(21)*(y(1:24,2)+y(2:25,6)-y(1:24,6));

Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];

Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
         y(17,10)-y(13,10); y(21,10)-y(17,10)];

 Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
 Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
        833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
        827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
        832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
        775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
        688451./5754509; 683498./5754509; 641894./5754509; 625054./5754509];
 Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 

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


 
 for i=1:21;
    continuous1(i)=pars(21)*(y(i,2)+y(i+4,6)-y(i,6));
 end 
 
 for i=1:21;
    continuous2(i)=pars(19)*(y(i,3)+y(i+4,7)-y(i,7));
 end 
 
 for i=1:9;
    continuous3(i)=pars(19)*(y(i+4,4)+y(i+8,8)-y(i+4,8));
 end 
 
 for i=1:24;
    continuous4(i)=pars(21)*(y(i,2)+y(i+1,6)-y(i,6));
 end 
 
 for i=1:13;
    continuous5(i)=y(i+4,9)-y(i,9);
 end 
 
 for i=1:17;
     continuous6(i)=y(i+4,10)-y(i,10);
 end 
 
 
 
 figure(1)
 hold all
 plot(t,y(:,1),'k-','LineWidth',3)
 set(gca, 'FontSize',16)
 xlabel('Year')
 ylabel('Susceptibles')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
           
 figure(2)
 hold all
 plot(t,y(:,2),'b-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Prescription Users')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
           
          
 figure(3)
 hold all
 plot(t,y(:,3),'r-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Opioid addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
          
                   
 figure(4)
 hold all
 plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})         
 
      
 figure(5)
 hold all
 plot(t,y(:,5),'Color', [0.7,0,0.7],'LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) 
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 
 
 
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





