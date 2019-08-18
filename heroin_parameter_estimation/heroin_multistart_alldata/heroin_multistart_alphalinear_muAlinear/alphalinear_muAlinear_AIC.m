m=-0.016315043366142;
beta_A=0.004576784574296;
beta_P=0.001961100711192;
theta_1=5.087528802931296e-04;
epsilon=2.500714817046142;
mu=0.00868;
mu_H=0.0507;
gamma=0.003630552621414;
theta_2=0.893900278857187;
sigma=0.889702060900079;
zeta=0.518669681971285;
theta_3=1.221956397540627;
nu=4.699736045487307e-04;
omega=0.0000000001;
b=0.288729997909935;
d=0.003219684445889;
e=0.009627721745951;

pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,d,e];


% Final time N; will run from beginning of 2013 to beginning of 2019 where t=0 represents 2013
% and t=6 represents 2019, with spacing (N-0)/(25-1)=0.25 between the points to represent quarters of a year:
N = 6; 
tspan=linspace(0,N,25);


% Initial conditions
P0=0.083700671037689;
A0=0.005902511929057;
H0=4.583608179759002e-04;
R0=0.001305505480884;
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

 alpha=m*t+b;
 
 
 
  % Making sure S+P+A+H+R=1
 total=y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5);
 
Estim1=[y(1,2)+y(5,6)-y(1,6); y(5,2)+y(9,6)-y(5,6); y(9,2)+y(13,6)-y(9,6);...
         y(13,2)+y(17,6)-y(13,6); y(17,2)+y(21,6)-y(17,6); y(21,2)+y(25,6)-y(21,6)];
     
Estim2=[y(1,3)+y(5,7)-y(1,7); y(5,3)+y(9,7)-y(5,7); y(9,3)+y(13,7)-y(9,7);...
        y(13,3)+y(17,7)-y(13,7); y(17,3)+y(21,7)-y(17,7); y(21,3)+y(25,7)-y(21,7)];
    
Estim3=[y(5,4)+y(9,8)-y(5,8); y(9,4)+y(13,8)-y(9,8); y(13,4)+y(17,8)-y(13,8)];
     
Estim4=y(1:24,2)+y(2:25,6)-y(1:24,6);
      
Estim5=[y(5,9)-y(1,9); y(9,9)-y(5,9); y(13,9)-y(9,9); y(17,9)-y(13,9)];
       
Estim6=[y(5,10)-y(1,10); y(9,10)-y(5,10); y(13,10)-y(9,10);...
        y(17,10)-y(13,10); y(21,10)-y(17,10)];
 % Actual Data for years 2013-2018
Data1=[1825910./5519417; 1805325./5559702; 1800614./5602187; 1744766./5648259; 1620955./5702475; 1455093./5754509];
 
Data2=[43418./5519417; 42928./5559702; 42816./5602187; 37464./5648259; 34805./5702475; 31244./5754509];
  
Data3=[7560./5559702; 7560./5602187; 10260./5648259];
  
Data4=[847077./5519417; 860931./5519417; 864889./5519417; 847077./5519417;...
        833223./5559702; 851035./5559702; 861921./5559702; 841140./5559702;...
        827285./5602187; 852025./5602187; 855983./5602187; 845098./5602187;...
        832085./5648259; 821189./5648259; 793453./5648259; 775622./5648259;...
        775622./5702475; 764726./5702475; 739961./5702475; 706282./5702475;...
        688451./5754509; 683498./5754509; 641894./5754509; 625054./5754509];
    
Data5=[348./5519417; 381./5559702; 463./5602187; 551./5648259];
 
Data6=[116./5519417; 216./5559702; 374./5602187; 554./5648259; 811./5702475];

Diff1=Estim1-Data1;
Diff2=Estim2-Data2;
Diff3=Estim3-Data3;
Diff4=Estim4-Data4;
Diff5=Estim5-Data5;
Diff6=Estim6-Data6; 

value=norm(Diff1,2)^2+norm(Diff2,2)^2+norm(Diff3,2)^2+norm(Diff4,2)^2+norm(Diff5,2)^2+norm(Diff6,2)^2;
value2=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);

fprintf('AIC score')
% sum of squares value, 48 data points, 19 parameters estimating+1 for sum
% of squares value 
disp(AIC(value,48,20))
      
function alpha = a(t,pars)
    alpha = pars(1)*t+pars(15);
end
           
function mu_A = muA(t,pars)
    mu_A = pars(16)*t+pars(17);
end


function f = HeroinModel(t,y,pars)
f=zeros(8,1);
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



% "Corrected AIC value": W is number of data points, K is number of parameters fitting +1 (since 
% least squares regression is also estimating the objective function value
% fval)
function aic=AIC(value,W,K); 
    aic = W.*log(value./W)+2.*K+2*K*(K+1)./(W-K-1);
end


