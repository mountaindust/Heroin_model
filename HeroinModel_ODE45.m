function value = HeroinModel_ODE45(z)
%CORRECT?


% Final time 
N = 8;
T = N;
%Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
tspan=linspace(0,T,N+1);
global value
 
%Parameters
alpha=z(1); 
 
beta=z(2); 
 
xi=z(3);
 
theta_1=z(4);
 
epsilon=z(5);
 
mu=z(6);  
 
mu_A=z(7);   
 
mu_H=z(8);
 
gamma=z(9);   
 
theta_2=z(10); 
 
sigma_A=z(11);
 
zeta=z(12);
 
theta_3=z(13);
 
sigma_H=z(14);
 
nu=z(15);

%Initials
%MADE UP VALUES IN ORDER TO RUN CODE
S0=0.7; 
P0=0.1;
A0=0.1;
H0=.05;
R0=.05; 
X0=0;
initials = [S0,P0,A0,H0,R0,X0];



[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);
%CORRECT?


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  
 % FIRST WAY: In order to count the total number of individuals in P at some point throughout the year, 
 % we want to add the total number of individuals using prescriptions at the beginning of the
 % year to those who come into P at some point during the year, so for second piece, want to integrate alpha*S 
 % from one year to the next year using trapezoidal rule. 
 
%prescript=zeros(4,8);
%for i=4:8
% prescript(i) = y(i,2)+alpha*((y(i+1,1)+y(i,1))/2)*(t(i+1)-t(i)); 
% end
 
 %yearly output from the model as a fraction
%Estim=[prescript(4),prescript(5),prescript(6),prescript(7),prescript(8)];
       

%SECOND WAY: Using ODE X'(t) in order to account for new cases coming into prescription
%class
 prescript=zeros(4,8);
 for i=4:8
 prescript(i) = y(i,2)+y(i+1,6)-y(i,6); 
 end
 
%yearly output from the model as a fraction
 Estim=[prescript(4),prescript(5),prescript(6),prescript(7),prescript(8)];
       
% Proportions of population of prescription opioid users (MADE UP FRACTIONS
% FOR NOW)
 Data=[.1 .2 .3 .25 .4];
 

% the difference between estimated value and data 
 diff= Estim-Data;

 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2)) normalized by 
 value = norm(diff,2)./norm(Data); 
 
 
