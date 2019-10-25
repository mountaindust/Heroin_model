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
N = 6;
%tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
tspan=linspace(0,N,3000);


% Initial Conditions
P0=0.0949727450989279;
A0=0.00709742287302280;
H0=0.000464895055434927;
R0=0.00507229016950725;
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

 
 for i=1:21;
    continuous1(i)=y(i,2)+y(i+4,6)-y(i,6);
 end 
 
 for i=1:21;
    continuous2(i)=y(i,3)+y(i+4,7)-y(i,7);
 end 
 
 for i=1:9;
    continuous3(i)=y(i+4,4)+y(i+8,8)-y(i+4,8);
 end 
 
 for i=1:24;
    continuous4(i)=y(i,2)+y(i+1,6)-y(i,6);
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





