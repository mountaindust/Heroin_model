function value = HeroinModel_ODE45(z)


% Final time, note: don't want to run too long because dynamics
% can drastically change over a number of years, so here we will do 2013-2017:
N = 4; 
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

A0=z(16);

H0=z(17);

R0=z(18);


%Initials
%MADE UP VALUES IN ORDER TO RUN CODE
S0=1-0.1-z(16)-z(17)-z(18); 
P0=0.1;
A0=z(16);
H0=z(17);
R0=z(18); 
X0=0;
Z0=0;
initials = [S0,P0,A0,H0,R0,X0,Z0];



[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  Z=y(:,7);
  
  
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
       

 %SECOND WAY: In order to count the total number of individuals in P at some point throughout the year, 
 % we want to add the total number of individuals using prescriptions at the beginning of the
 % year to those who come into P at some point during the year: 
 % prescript(i) is number of prescription users at beginning of year(time step) +
 % the number of new cases that came in from X'=dy(6) ODE from that year until the beginning of the next year 
 % so adding in those who come into P at some point during the year by integrating y(6) ODE 
 % but just focusing in on the one year we care about (so have to subtract)
 % Going to run from 2013 to 2017, so "IC" will be at 2013 and want info
 % for 2014-2017 so want:
 
 prescript=zeros(1,4);
 for i=1:4
 prescript(i) = y(i,2)+y(i+1,6)-y(i,6); 
 end
 
%yearly output from the model as a fraction
 Estim1=[prescript(1),prescript(2),prescript(3),prescript(4)];
       
% Proportions of population of prescription opioid users (MADE UP FRACTIONS
% FOR NOW)
 Data1=[.1 .2 .3 .25];
 

% the difference between estimated value and data 
 diff1= Estim1-Data1;
 
 
 
 %Using ODE y(7) in order to account for new admissions coming into recovery class
 %(NOT total in recovery class), going to run from 2013-2015, so want:
 admitted=zeros(1,2);
 for i=1:2
 admitted(i) = y(i+1,7)-y(i,7);
 end
 
 %yearly output from the model as a fraction
 Estim2=[admitted(1), admitted(2)];
  
 %Proportions of population being admitted into recovery (MADE UP FRACTIONS FOR NOW)
 Data2=[.001,.0015];

 % the difference between estimated value and data 
 diff2=Estim2-Data2;
 
 %Proportion of opioid addicts in 2015
 addicts(2)=y(2,3); 
 %yearly output from the model as a fraction
 Estim3=[addicts(2)];
%Made up for now
 Data3=[.4];
 %the difference between estimated value and data
 diff3=Estim3-Data3;

 
 %Proportion of heroin/fentanyl addicts in 2015
 heroin(2)=y(2,4);
 Estim4=[heroin(2)];
 %Made up for now
 Data4=[.2];
 %the difference between estimated value and data
 diff4=Estim4-Data4;
 
 %the relative error that we are trying to minimize for ordinary least
 %squares: the sum of the squared errors (norm gives sum(diff.^2)^(1/2))
 %normalized by norm of the data
 value = norm(diff1,2)./norm(Data1)+ norm(diff2,2)./norm(Data2)+norm(diff3,2)./norm(Data3)+norm(diff4,2)./norm(Data4);