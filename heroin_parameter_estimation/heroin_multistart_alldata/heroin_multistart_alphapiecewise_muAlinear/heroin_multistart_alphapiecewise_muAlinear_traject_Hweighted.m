%File name: heroin_multistart_alphapiecewise_muAlinear_testing.m
clf;
clear all;

m=-0.00959081904910168;
beta_A=0.000547665937526567;
beta_P=3.54815882817162e-05;
theta_1=0.0768026069493661;
epsilon=2.57359174497727;
mu=0.00710; 
mu_H=0.0466; 
gamma=0.00502442423818114;
theta_2=0.102162927238600;
sigma=0.218107878500014;
zeta=0.199831940842269;
theta_3=10.0287953455473;
nu=0.000485210862964808;
omega=0.0000000001;
b=0.280632334424353;
c=-0.0255287384270066;
d=0.00121423487972979;
e=0.00834122202941407;

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c,d,e];
 

% Final time and last entry of tspan is # of equally spaced points from 0 to N 
N=9; %predict trajectory
tspan=linspace(0,N,12*N+1); %number of months
%N=6; %to check fval 
%tspan=linspace(0,N,25);


% Initial Conditions
P0=0.0872022972867530;
A0=0.00771457681712701;
H0=0.000937746008611740;
R0=1.40671155678552e-05;
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
 
value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+20.*norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6)

 
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
 
 %number of prescription opioid addict overdoses at any time point
 W=muA(t,pars).*y(:,3);
 %number of heroin/fentanyl overdoses at any time point
 Z=mu_H.*y(:,4);
 
 
% Totaling the proportion of addicts
 figure(1)
 hold all
 plot(t(1:73),y(1:73,3)+y(1:73,4),'Color',[0.1 0.4 1],'LineWidth',3);
 plot(t(74:12*N+1),y(74:12*N+1,3)+y(74:12*N+1,4),'Color',[0.1 0.4 1],'LineStyle', '--','LineWidth',3);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion')
 legend({'Total Addicts'},'FontSize',14, 'Location', 'northwest')
 xlim([0 N])
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9]) 
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017', '2018','2019', '2020', '2021', '2022'})
 
 

 % Totaling the proportion of addicts
 figure(2)
 hold all
 plot(t(1:73),y(1:73,3),'Color', [0, 0, 0.5],'LineWidth',3)
 plot(t(1:73),y(1:73,4),'Color', [0,0.9,0.7],'LineWidth',3)
 plot(t(74:12*N+1),y(74:12*N+1,3),'Color', [0,0,0.5],'LineStyle', '--','LineWidth',3)
 plot(t(74:12*N+1),y(74:12*N+1,4),'Color', [0,0.9,0.7],'LineStyle', '--','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion')
 legend({'Prescription Opioid Addicts', 'Heroin and Fentanyl Addicts'},'FontSize',14, 'Location', 'northwest')
 xlim([0 N])
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9]) 
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017', '2018','2019', '2020', '2021', '2022'})
 
 
 
 figure(3)
 hold all
 plot(t(1:73),W(1:73),'Color', 'm','LineWidth',3)
 plot(t(1:73),Z(1:73),'Color', [.5 0 .5],'LineWidth',3)
 plot(t(74:12*N+1),W(74:12*N+1),'Color', 'm','LineStyle', '--','LineWidth',3)
 plot(t(74:12*N+1),Z(74:12*N+1),'Color', [.5 0 .5],'LineStyle', '--','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion')
 legend({'Prescription Opioid Addict Overdose Deaths', 'Heroin and Fentanyl Overdose Deaths'},'FontSize',14, 'Location', 'northwest')
 xlim([0 N])
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9]) 
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017', '2018','2019', '2020', '2021', '2022'})
 
 
 figure(4)
 hold all
 plot(t(1:73),W(1:73)+Z(1:73),'Color',[0.4, 0, 0.8],'LineWidth',3)
 plot(t(74:12*N+1),W(74:12*N+1)+Z(74:12*N+1),'Color',[0.4, 0, 0.8],'LineStyle', '--','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion')
 legend({'Total Addict Overdose Deaths'},'FontSize',14, 'Location', 'northwest')
 xlim([0 N])
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9]) 
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017', '2018','2019', '2020', '2021', '2022'})
 
 
 %{
 figure(5)
 hold all
 plot(t(1:73),y(1:73,9),'Color', 'm','LineWidth',3)
 plot(t(1:73),y(1:73,10),'Color', [.5 0 .5],'LineWidth',3)
 plot(t(74:12*N+1),y(74:12*N+1,9),'Color', 'm','LineStyle', '--','LineWidth',3)
 plot(t(74:12*N+1),y(74:12*N+1,10),'Color', [.5 0 .5],'LineStyle', '--','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion')
 legend({'Prescription Opioid Addict Overdose Deaths', 'Heroin and Fentanyl Overdose Deaths'},'FontSize',14, 'Location', 'northwest')
 xlim([0 N])
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9]) 
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017', '2018','2019', '2020', '2021', '2022'})
 
 
 figure(6)
 hold all
 plot(t(1:73),y(1:73,9)+y(1:73,10),'Color',[0.4, 0, 0.8],'LineWidth',3)
 plot(t(74:12*N+1),y(74:12*N+1,9)+y(74:12*N+1,10),'Color',[0.4, 0, 0.8],'LineStyle', '--','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion')
 legend({'Total Addict Overdose Deaths'},'FontSize',14, 'Location', 'northwest')
 xlim([0 N])
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9]) 
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017', '2018','2019', '2020', '2021', '2022'})
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


