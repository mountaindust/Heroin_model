clf;
clear all;

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
g=-0.038616505394549;
h=0.000977482526657751;

pars=[beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,g,h];

% Final time and last entry of tspan is # of equally spaced points from 0 to N (quarterly linspace)
N = 3;
tspan=linspace(0,N,4);
%For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);


% Initial Conditions
P0=0.0585;
A0=0.0037;
H0=0.00597;
R0=0.00751;
S0=1-P0-A0-H0-R0;
J0=0;
K0=0;

initials = [S0;P0;A0;H0;R0;J0;K0];

[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  J=y(:,6);
  K=y(:,7);
  
% Making sure S+P+A+H+R=1
  total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);

% format long
% fprintf('alpha values')
% disp(a(0,pars))
% disp(a(1,pars))
% disp(a(2,pars))
% disp(a(3,pars))
% 
% 
% format long
% fprintf('muA values')
% disp(muA(0,pars))
% disp(muA(1,pars))
% disp(muA(2,pars))
% disp(muA(3,pars))


 disp('percent decrease of alpha')
 (0.069883058240252-a(3,pars))*100/0.069883058240252
 disp('value of alpha')
 a(3,pars)
 
 disp('value of total A overdoses with this alpha value')
 J(4,1)
 disp('percent decrease of total A overdoses from baseline 2023 value with this new alpha value')
 (1.359799130226335e-04-J(4,1))*100/1.359799130226335e-04
 
 disp('value of total H overdoses with this new alpha value')
 K(4,1)
 disp('percent decrease of total H overdoses from baseline 2023 value with this new alpha value')
 (0.001334252863060-K(4,1))*100/0.001334252863060
 

 

 
 figure(1)
 hold all
 plot(t,y(:,3),'r-','LineWidth',3)
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
 
 
function alpha = a(t,pars)
    alpha =  -0.00559565027929907*3.25+0.270110337915851+0.0269690987063522*3.25-0.0269690987063522*7+pars(14)*t;

end

           
function mu_A = muA(t,pars)
    mu_A = 0.000977482526657751*7+0.00883138792481281+pars(15)*t;
end


function f = HeroinModel(t,y,pars)
f=zeros(5,1);
f(1)=-a(t,pars)*y(1)-pars(1)*y(1)*y(3)-pars(2)*y(1)*y(2)-pars(3)*y(1)*y(4)+pars(4)*y(2)+pars(5)*(y(2)+y(5))+(pars(5)+muA(t,pars))*y(3)+(pars(5)+pars(6))*y(4);
f(2)=a(t,pars)*y(1)-pars(4)*y(2)-pars(7)*y(2)-pars(8)*y(2)*y(4)-pars(5)*y(2);
f(3)=pars(7)*y(2)+(pars(9)*y(5)*y(3))/(y(3)+y(4)+pars(13))+pars(1)*y(1)*y(3)+pars(2)*y(1)*y(2)-pars(10)*y(3)-pars(11)*y(3)*y(4)-pars(5)*y(3)-muA(t,pars)*y(3);
f(4)=pars(3)*y(1)*y(4)+pars(8)*y(2)*y(4)+pars(11)*y(3)*y(4)+(pars(9)*y(5)*y(4))/(y(3)+y(4)+pars(13))-pars(12)*y(4)-(pars(5)+pars(6))*y(4);
f(5)=pars(10)*y(3)+pars(12)*y(4)-(pars(9)*y(5)*y(3))/(y(3)+y(4)+pars(13))-(pars(9)*y(5)*y(4))/(y(3)+y(4)+pars(13))-pars(5)*y(5);
f(6)=muA(t,pars)*y(3);
f(7)=pars(6)*y(4);
end








