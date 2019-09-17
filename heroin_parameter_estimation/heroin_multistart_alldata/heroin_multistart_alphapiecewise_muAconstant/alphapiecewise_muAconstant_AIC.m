%Parameters
%slope of alpha 
m=-0.003663788942920;%-0.005885298827295;
beta_A=2.558765885155569e-04;%2.928459237447347e-04;
beta_P=2.876580746292427e-05;%3.031594496251287e-05;
theta_1=0.278195789336426;%0.149588835714450;
epsilon=2.614429685431156;%2.497686906039834;
mu=0.00710; 
mu_A=0.00884;
mu_H=0.0466;      
gamma=0.005017113475167;%0.005018446561243;
theta_2=0.374414127368574;%0.995846977898647;
sigma=0.100942113632810;%0.100824905478476;
zeta=0.097881888956038;%0.173579520334630;
theta_3=18.800770308123713;%1.801581240297609;
nu=5.113823127121602e-04;%3.194810546447885e-04;
omega=0.0000000001;
%y-intercept of alpha 
b=0.266207134284557;%0.269343669604923;
c=-0.027384577579411;%-0.026374597945385;



pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b,c];
 

% Final time and last entry of tspan is # of equally spaced points from 0 to N 
N = 6;
tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);


% Initial Conditions
P0=0.095266364858678;%0.092661491071360;
A0=0.007618716405553;%0.007594784328536;
H0=4.688115365096090e-04;%7.052316770178626e-04;
R0=3.859190686154302e-05;%0.002883057248950;
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
 
 
 value=norm(Diff1,2)^2+norm(Diff2,2)^2+norm(Diff3,2)^2+norm(Diff4,2)^2+norm(Diff5,2)^2+norm(Diff6,2)^2;
 value2=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6);

 format short
 fprintf('AIC score')
 % sum of squares value, 48 data points, 18 parameters estimating+1 for sum
 % of squares value 
 disp(AIC(value,48,19))
 
 
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


 
% "Corrected AIC value": W is number of data points, K is number of parameters fitting +1 (since 
% least squares regression is also estimating the objective function value
% fval)
function aic=AIC(value,W,K); 
    aic = W.*log(value./W)+2.*K+2*K*(K+1)./(W-K-1);
end
