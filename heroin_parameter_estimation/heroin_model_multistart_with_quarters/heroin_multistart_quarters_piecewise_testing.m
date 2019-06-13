%File name: heroin_multistart_quarters_testing.m
clf;
clear all;

%Parameters
%slope of alpha 
m=-0.0032;%-0.00365;
beta_A=0.0041;%0.00344; 
beta_P=0.000353;%0.000353; 
theta_1=0.000493;%0.0005;
epsilon=2.47;%;2.50;
mu=0.00868; 
mu_A=0.00870;      
mu_H=0.0507;
gamma=0.0013;%0.00135;
theta_2=0.0745;%0.142; 
sigma=0.0353;%0.0771;
zeta=0.3775;%0.329;
theta_3=0.1773;%0.122; 
nu=0.0876;%0.0458;
omega=0.0000000001;
%y-intercept of alpha 
b=0.278;%0.269;
c=-0.0322;%-0.0299;


pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c];
 

% Final time and last entry of tspan is # of equally spaced points from 0 to N 
N = 6;
tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,1000);


% Initial Conditions
P0=0.0987;%0.0971;
A0=0.0062;%0.00645;
H0=0.000812;%0.000847;
R0=0.0563;%0.0220;
S0=1-P0-A0-H0-R0;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];

[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);
  
  %alpha=m*t+b;
  
% Making sure S+P+A+H+R=1
  total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);

disp(a(0,pars))
disp(a(1,pars))
disp(a(2,pars))
disp(a(3,pars))
disp(a(4,pars))
disp(a(5,pars))
disp(a(6,pars))

% Comment out if don't need objective function value 

Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6)];
 
Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7)];       

Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8); y(13,4)+y(17,8)-y(13,8)];
 
Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);



 Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475];
 Diff1=Estim1-Data1; 
 Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475];
 Diff2=Estim2-Data2;  
 Data3=[7560./5559702; 7560./5602187; 10260./5648259];
 Diff3=Estim3-Data3;
 Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
        833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
        827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
        832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
        775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
        688502./5754509; 683722./5754509; 641942./5754509; 625162./5754509];
 Diff4=Estim4-Data4;
 
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)

 
 % ODE solutions plotted separately shown all together
 figure(1)
         
           subplot(2,2,1);plot(t,y(:,2),'b-','LineWidth',3)
           subplot(2,2,1);xlabel('Year')
           subplot(2,2,1);ylabel('Prescription Users')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018','2019'})
 
 
          
           
           subplot(2,2,2);plot(t,y(:,3),'r-','LineWidth',3)
           subplot(2,2,2);xlabel('Year')
           subplot(2,2,2);ylabel('Opioid Addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})
 
           
           subplot(2,2,3);plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
           subplot(2,2,3);xlabel('Year')
           subplot(2,2,3);ylabel('Heroin/Fentanyl Addicts')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})
 
          
           subplot(2,2,4);plot(t,y(:,5) , 'Color', [0.7,0,0.7],'LineWidth',3)
           subplot(2,2,4);xlabel('Year')
           subplot(2,2,4);ylabel('Stably Recovered Individuals')
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018', '2019'})  
               
                 
 % ODE Solutions for P, A, H plotted all together
 figure(2)
           plot(t,y(:,3),'r-','LineWidth',3);
           hold all
           plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3);
           xlabel('Year')
           ylabel('Size of Addicted Populations');
           set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
           set(gca, 'fontsize',10)
           xtickangle(90)
           legend({'A','H'}, 'FontSize', 14)
           set(gca,'XLim',[0 N])
           set(gca,'xticklabel',{'2013', '2014', '2015', '2016',...
                       '2017', '2018', '2019'})
 

 
 figure(3)
 hold all
 plot(t,y(:,1),'k-','LineWidth',3)
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
           
 figure(4)
 hold all
 plot(t,y(:,2),'b-','LineWidth',3)
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
           
           
 figure(5)
 hold all
 plot(t,y(:,3),'r-','LineWidth',3)
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
          
                   
 figure(6)
 hold all
 plot(t,y(:,4),'Color', [0,0.9,0],'LineWidth',3)
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})         
 
      
 figure(7)
 hold all
 plot(t,y(:,5),'Color', [0.7,0,0.7],'LineWidth',3)
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 ]) 
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019'})
 
 %Comment out if don't need fitting plots 

              
 % Data points from proportion that is in P at some point in the year and corresponding ODE solution points 
 figure(8)
 hold all
 z1 = linspace(0,4,5); %defines mesh where going to plot Estim1, Data1 values 
 scatter(z1, Estim1, 50, 'o');
 scatter(z1, Data1, 50, 'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P')
  legend({'ODE solution', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})
 

 
  % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
 figure(9)
 hold all
 z2 = linspace(0,4,5);
 scatter(z2, Estim2, 50, 'o');
 scatter(z2, Data2, 50, 'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A')
 legend({'ODE solution', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 3 4 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2013','2014','2015','2016', '2017'})




 % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
 figure(10)
 hold all
 z3 = linspace(0,2,3);
 scatter(z3, Estim3, 50,'o');
 scatter(z3, Data3, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H')
  legend({'ODE solution', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 1 2 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 

 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(11)
 hold all
 z4 = linspace(0,23,24);
 scatter(z4, Estim4, 50, 'o');
 scatter(z4, Data4, 50,'x');
 set(gca, 'fontsize',10)
 xlabel('Quarter')
 ylabel('Proportion in P')
 legend({'ODE solution', 'Data'},'FontSize', 14)
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


end


 