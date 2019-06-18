%File name: heroin_multistart_quarters_testing.m

%Parameters
%slope of alpha 
m=-.0156;
beta_A=0.00235; 
beta_P=0.000141; 
theta_1=0.000507;
epsilon=2.54;
mu=0.00868; 
mu_A=0.00870;      
mu_H=0.0507;
gamma=0.00115;
theta_2=0.0370; 
sigma=0.0284;
zeta=0.265;
theta_3=3.51; 
nu=0.00657;
omega=0.0000000001;
%y-intercept of alpha 
b=0.303; 

%{
% For R_0 checking:
alpha=0.2; 
beta_A=0.000273; 
beta_P=0; 
theta_1=0.0003;
epsilon=1.5;
mu=0.00868; 
mu_A=0.00775;   
mu_H=0.0271;
gamma=0;   
theta_2=3*theta_1; 
sigma=0.7;
zeta=0.25;
theta_3=16*theta_1; 
nu=0.1;
omega=0.0000000001;
%}

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];

% Final time and last entry of tspan is # of equally spaced points from 0 to N 
N = 6;
tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,100);

% Initial Conditions
P0=0.0835;
A0=0.00671;
H0=0.000874;
R0=0.0509;
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
  
  alpha=m*t+b;
  
% Making sure S+P+A+H+R=1
  total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);

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

 % FOR CALCULATING AIC, use this objective function value (not relative error and norm is squared)
value=norm(Diff1,2)^2+norm(Diff2,2)^2+norm(Diff3,2)^2+norm(Diff4,2)^2

% sum of squares value, 37 data points, 17 parameters estimating+1 for sum
% of squares value 
disp(AIC(value,37,18))


function f = HeroinModel(t,y,pars)
f=zeros(8,1);
f(1)=-(pars(1)*t+pars(16))*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=(pars(1)*t+pars(16))*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);


% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
f(6) = (pars(1)*t+pars(16))*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));


end


% "Corrected AIC value": W is number of data points, K is number of parameters fitting +1 (since 
% least squares regression is also estimating the objective function value
% fval)
function aic=AIC(value,W,K); 
    aic = W.*log(value./W)+2.*K+2*K*(K+1)./(W-K-1);
end

