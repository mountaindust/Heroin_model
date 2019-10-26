%File name: heroin_multistart_alphapiecewise_muAconstant_testing.m
clf;
clear all;

%Parameters
%slope of alpha 
m=-0.00759829901019442;
beta_A=6.01126909926855e-05;
beta_P=1.37312211291756e-05;
theta_1=0.284401807266879;
epsilon=2.56273677041097;
mu=0.00710; 
mu_A=0.00884;
mu_H=0.0466; 
gamma=0.00500339664720886;
theta_2=0.355814628540795;
sigma=0.100184804571766;
zeta=0.100402513311739;
theta_3=18.0635609734696;
nu=0.000181332211254903;
omega=0.0000000001;
%y-intercept of alpha 
b=0.275305858322444;
c=-0.0264943422500030;


pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c];
 

% Final time and last entry of tspan is # of equally spaced points from 0 to N 
N = 6;
%tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);
% For smooth data plots; 73 marks to represent 72 months in 6 years 
tspan=linspace(0,N,73);



% Initial Conditions
P0=0.0922293025105621;
A0=0.00763634664877567;
H0=0.000470731911097036;
R0=1.85917132509588e-05;
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
  
  %alpha=m*t+b;
  
% Making sure S+P+A+H+R=1
  total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);

format short
disp(a(0,pars))
disp(a(1,pars))
disp(a(2,pars))
disp(a(3,pars))
disp(a(3.25,pars))
disp(a(4,pars))
disp(a(5,pars))
disp(a(6,pars))

% Comment out if don't need objective function value 

Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6); y(21,2)+y(25,6)-y(21,6)];
 
Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7); y(21,3)+y(25,7)-y(21,7)];      

Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8); y(13,4)+y(17,8)-y(13,8)];
 
Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);

Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];

Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
        y(17,10)-y(13,10); y(21,10)-y(17,10)];


 Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
 Diff1=Estim1-Data1; 
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
 Diff2=Estim2-Data2;  
 Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 Diff3=Estim3-Data3;
 Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
        833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
        827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
        832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
        775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
        688451./5754509; 683498./5754509; 641894./5754509; 625054./5754509];
 Diff4=Estim4-Data4;
 Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
 Diff5=Estim5-Data5;
 Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];
 Diff6=Estim6-Data6;
 
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6)

 
 for i=1:61;
    continuous1(i)=y(i,2)+y(i+12,6)-y(i,6);
 end 
 
 for i=1:61;
    continuous2(i)=y(i,3)+y(i+12,7)-y(i,7);
 end 
 
 for i=1:25;
    continuous3(i)=y(i+12,4)+y(i+24,8)-y(i+12,8);
 end 
 
 for i=1:70;
    continuous4(i)=y(i,2)+y(i+3,6)-y(i,6);
 end 
 
 for i=1:37;
    continuous5(i)=y(i+12,9)-y(i,9);
 end 
 
 for i=1:49;
     continuous6(i)=y(i+12,10)-y(i,10);
 end 
 
 

 figure(1)
 hold all
 z1 = linspace(0,60,6); %defines mesh where going to plot Estim1, Data1 values 
 z7 = linspace(0,60,61);
 plot(z7,continuous1,'k-','LineWidth',1.3);
 %scatter(z1, Estim1, 100, 'o');
 scatter(z1, Data1, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 12 24 36 48 60])
 set(gca,'XLim',[0 60])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 

 
 % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
 figure(2)
 hold all
 z2 = linspace(0,60,6);
 z8 = linspace(0,60,61);
 plot(z8,continuous2,'k-','LineWidth',1.3);
 %scatter(z2, Estim2, 100, 'o');
 scatter(z2, Data2, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in A')
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 12 24 36 48 60])
 set(gca,'XLim',[0 60])
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017','2018'})




 % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
 figure(3)
 hold all
 z3 = linspace(0,24,3);
 z9 = linspace(0,24,25);
 plot(z9,continuous3,'k-','LineWidth',1.3);
 %scatter(z3, Estim3, 100,'o');
 scatter(z3, Data3, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in H')
 legend({'Model simulation', 'Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 12 24])
 set(gca,'XLim',[0 24])
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 

 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(4)
 hold all
 z4 = linspace(0,69,24);
 z10 = linspace(0,69,70);
 plot(z10,continuous4,'k-','LineWidth',1.3);
 %scatter(z4, Estim4, 100, 'o');
 scatter(z4, Data4, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Quarter')
 ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69])
 xtickangle(90)
 set(gca,'XLim',[0 69])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
 
  
 figure(5)
 hold all
 z5 = linspace(0,36,4); %defines mesh where going to plot Estim5, Data5 values
 z11 = linspace(0,36,37);
 plot(z11,continuous5,'k-','LineWidth',1.3);
 %scatter(z5, Estim5, 100,'o');
 scatter(z5, Data5, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion overdose from A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 12 24 36])
 set(gca,'XLim',[0 36])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})
 
 
 figure(6)
 hold all
 z6 = linspace(0,48,5);
 z12 = linspace(0,48,49);
 plot(z12,continuous6,'k-','LineWidth',1.3);
 %scatter(z6, Estim6, 100,'o');
 scatter(z6, Data6, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion overdose from H') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 12 24 36 48])
 set(gca,'XLim',[0 48])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})  
 

 

 
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


 
