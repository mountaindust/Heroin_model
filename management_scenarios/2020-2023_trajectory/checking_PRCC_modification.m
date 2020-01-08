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
%To know values each year 2020-2023 (modified alpha and muA to be t+7
%because of time-dependence) 
tspan=linspace(0,3,4);
%For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);
%To run to 2020 and know values each year
%tspan=linspace(0,7,8);
%tspan=linspace(0,7,3000);


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

format long
fprintf('muA values')
disp(muA(0,pars))
disp(muA(1,pars))
disp(muA(2,pars))
disp(muA(3,pars))
 
 
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
 
  
%{
figure(1)
           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',3)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)  
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018','2019'})
           set(gca,'box','off')
 
          
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',3)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})
           set(gca,'box','off')
           
           subplot(2,2,3);plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})
           set(gca,'box','off')
          
           subplot(2,2,4);plot(t,y(:,5),'Color',[0.7,0,0.7],'LineWidth',3)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Stably Recovered Individuals')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})  
           set(gca,'box','off')
                 
                   
 figure(2)
 hold all
 plot(t,y(:,3),'r-','LineWidth',3);
 plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3);       
 xlabel('Year')
 ylabel('Size of Addicted Populations')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})                   
 legend({'A','H'}, 'FontSize', 14)                
                   
             
 % Data points from proportion that is in P at some point in the year and corresponding ODE solution points 
 figure(8)
 hold all
 z1 = linspace(0,5,6); %defines mesh where going to plot Estim1, Data1 values 
 z7 = linspace(0,5,21);
 plot(z7,continuous1,'k-');
 %scatter(z1, Estim1, 100, 'o');
 scatter(z1, Data1, 100, 'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 5])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 

 
 % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
 figure(9)
 hold all
 z2 = linspace(0,5,6);
 z8 = linspace(0,5,21);
 plot(z8,continuous2,'k-');
 %scatter(z2, Estim2, 100, 'o');
 scatter(z2, Data2, 100, 'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A')
 legend({'Model simulation', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 5])
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017','2018'})




 % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
 figure(10)
 hold all
 z3 = linspace(0,2,3);
 z9 = linspace(0,2,9);
 plot(z9,continuous3,'k-');
 %scatter(z3, Estim3, 100,'o');
 scatter(z3, Data3, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H')
 legend({'Model simulation', 'Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 ])
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 

 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(11)
 hold all
 z4 = linspace(0,23,24);
 z10 = linspace(0,23,24);
 plot(z10,continuous4,'k-');
 %scatter(z4, Estim4, 100, 'o');
 scatter(z4, Data4, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in P')
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
 

 figure(12)
 hold all
 z5 = linspace(0,3,4); %defines mesh where going to plot Estim5, Data5 values
 z11 = linspace(0,3,13);
 plot(z11,continuous5,'k-');
 %scatter(z5, Estim5, 100,'o');
 scatter(z5, Data5, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 3])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})
 
 
 
 figure(13)
 hold all
 z6 = linspace(0,4,5); %defines mesh where going to plot Estim6, Data6 values 
 z12 = linspace(0,4,17);
 plot(z12,continuous6,'k-');
 %scatter(z6, Estim6, 100,'o');
 scatter(z6, Data6, 100,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion overdose from H') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})  
 
 
 figure(14)
 set(gcf, 'Position',  [1, 1, 700, 300])

 z1 = linspace(0,5,6); %defines mesh where going to plot Estim1, Data1 values 
 subplot(1,2,1);scatter(z1, Estim1, 80, 'o');
 hold on
 subplot(1,2,1);scatter(z1, Data1, 80, 'x');
 set(gca, 'fontsize',10)
 subplot(1,2,1);xlabel('Year')
 subplot(1,2,1);ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 10)
 set(gca, 'xtick', [ 0 1 2 3 4 5])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 
 
 z4 = linspace(0,23,24);
 subplot(1,2,2);scatter(z4, Estim4, 80, 'o');
 hold on
 subplot(1,2,2);scatter(z4, Data4, 80,'x');
 set(gca, 'fontsize',10)
 subplot(1,2,2);xlabel('Quarter')
 subplot(1,2,2);ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 10)
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23])
 xtickangle(90)
 set(gca,'XLim',[0 23])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
 
                   
 figure(15)
 set(gcf, 'Position',  [1, 1, 700, 300])  
 z2 = linspace(0,5,6);
 subplot(1,2,1);scatter(z2, Estim2, 80, 'o');
 hold on
 subplot(1,2,1);scatter(z2, Data2, 80, 'x');
 set(gca, 'fontsize',10)
 subplot(1,2,1);xlabel('Year')
 subplot(1,2,1);ylabel('Proportion in A')
 legend({'Model simulation', 'Data'},'FontSize', 10)
 set(gca, 'xtick', [ 0 1 2 3 4 5])
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017','2018'})
 
 
 z5 = linspace(0,3,4); %defines mesh where going to plot Estim5, Data5 values
 subplot(1,2,2);scatter(z5, Estim5, 80,'o');
 hold on
 subplot(1,2,2);scatter(z5, Data5, 80,'x');
 set(gca, 'fontsize',10)
 subplot(1,2,2);xlabel('Year')
 subplot(1,2,2);ylabel('Proportion overdose from A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 10, 'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 3])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})
 

 
 
 
 figure(16)
 set(gcf, 'Position',  [1, 1, 700, 300])
 
 z3 = linspace(0,2,3);
 subplot(1,2,1);scatter(z3, Estim3, 80,'o');
 hold on
 subplot(1,2,1);scatter(z3, Data3, 80,'x');
 set(gca, 'fontsize',10)
 subplot(1,2,1);xlabel('Year')
 subplot(1,2,1);ylabel('Proportion in H')
 legend({'Model simulation', 'Data'},'FontSize', 10, 'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 ])
 set(gca,'xticklabel',{'2014', '2015', '2016'})
                   
                   
                   
 z6 = linspace(0,4,5);
 subplot(1,2,2);scatter(z6, Estim6, 80,'o');
 hold on 
 subplot(1,2,2);scatter(z6, Data6, 80,'x');
 subplot(1,2,2);xlabel('Year')
 subplot(1,2,2);ylabel('Proportion overdose from H')
 legend({'Model simulation', 'Data'},'FontSize', 10, 'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'}) 
 
 
 
 
 
 
 figure(17)
 set(gcf, 'Position',  [1, 1, 700, 300])

 z1 = linspace(0,5,6); %defines mesh where going to plot Estim1, Data1 values 
 z7 = linspace(0,5,21);
 subplot(1,2,1); plot(z7,continuous1,'k-');
 hold on
 subplot(1,2,1);scatter(z1, Data1, 80, 'x');
 set(gca, 'fontsize',10)
 subplot(1,2,1);xlabel('Year')
 subplot(1,2,1);ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 10)
 set(gca, 'xtick', [ 0 1 2 3 4 5])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 set(gca,'box','off')
 
 
 z4 = linspace(0,23,24);
 z10 = linspace(0,23,24);
 subplot(1,2,2);plot(z10,continuous4,'k-');
 hold on
 subplot(1,2,2);scatter(z4, Data4, 80,'x');
 set(gca, 'fontsize',10)
 subplot(1,2,2);xlabel('Quarter')
 subplot(1,2,2);ylabel('Proportion in P')
 legend({'Model simulation', 'Data'},'FontSize', 10)
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23])
 xtickangle(90)
 set(gca,'XLim',[0 23])
 set(gca,'box','off')
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
 
 
 
 figure(18)
 set(gcf, 'Position',  [1, 1, 700, 300])
 
 z2 = linspace(0,5,6);
 z8 = linspace(0,5,21);
 subplot(1,2,1);plot(z8,continuous2,'k-');
 hold on
 subplot(1,2,1);scatter(z2, Data2, 80, 'x');
 set(gca, 'fontsize',10)
 subplot(1,2,1);xlabel('Year')
 subplot(1,2,1);ylabel('Proportion in A')
 legend({'Model simulation', 'Data'},'FontSize', 10)
 set(gca, 'xtick', [ 0 1 2 3 4 5])
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017','2018'})
 set(gca,'box','off')
 

 z5 = linspace(0,3,4); %defines mesh where going to plot Estim5, Data5 values
 z11 = linspace(0,3,13);
 subplot(1,2,2); plot(z11,continuous5,'k-');
 hold on
 subplot(1,2,2);scatter(z5, Data5, 80,'x');
 set(gca, 'fontsize',10)
 subplot(1,2,2);xlabel('Year')
 subplot(1,2,2);ylabel('Proportion overdose from A') % at some point during the year
 legend({'Model simulation', 'Data'},'FontSize', 10,'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 3 ])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})
 set(gca,'box','off')  
 
 
 
 figure(19)
 set(gcf, 'Position',  [1, 1, 700, 300])
                
 z3 = linspace(0,2,3);
 z9 = linspace(0,2,9);
 subplot(1,2,1); plot(z9,continuous3,'k-');
 hold on
 subplot(1,2,1);scatter(z3, Data3, 80,'x');
 set(gca, 'fontsize',10)
 subplot(1,2,1);xlabel('Year')
 subplot(1,2,1);ylabel('Proportion in H')
 legend({'Model simulation', 'Data'},'FontSize', 10,'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 ])
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 set(gca,'box','off')
                 
                   
 z6 = linspace(0,4,5);
 z12 = linspace(0,4,17);
 subplot(1,2,2);plot(z12,continuous6,'k-');
 hold on 
 subplot(1,2,2);scatter(z6, Data6, 80,'x');
 subplot(1,2,2);xlabel('Year')
 subplot(1,2,2);ylabel('Proportion overdose from H')
 legend({'Model simulation', 'Data'},'FontSize', 10,'Location','northwest')
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'}) 
 set(gca,'box','off')
 
 
 
 
 %}
 
 
 
 
 
 
function alpha = a(t,pars)
    alpha = pars(1)*3.25+pars(15)-pars(16)*3.25+pars(16)*(t+7);
 
end

           
function mu_A = muA(t,pars)
    mu_A = pars(17)*(t+7)+pars(18);
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





