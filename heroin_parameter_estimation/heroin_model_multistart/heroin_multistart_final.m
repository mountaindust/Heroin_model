%File name: heroin_multistart_final.m 

clf

% We wish to estimate the parameter vector (7 parameters)
% x =[alpha,theta_1,epsilon,gamma,sigma,zeta,nu]
%Give ranges on each of the parameters 
LowerBounds=[0.00001 0.00001 0.00001 0.00001 0.00001 0.00001 0.00001 ];
UpperBounds=[2  2  4  2  2  2  2 ];

%Initial starting points for parameters, starting in the middle of each of the ranges
xstart=0.5*(LowerBounds + UpperBounds); 

% Create MultiStart problem using optimization function fmincon;
% x0 is xstart, objective is what we are trying to minimize which comes from 
% value = HeroinModel_ODE45(z) = fval(x) as output
problem=createOptimProblem('fmincon','objective',@HeroinModel_ODE45,...
         'x0', xstart,...
         'lb',LowerBounds,...
         'ub',UpperBounds);

problem.options=optimoptions(problem.options, 'MaxFunEvals',99999,'MaxIter',99999);

% Define a multistart problem; results are reported after each local solver run, in addition to the final summary
ms=MultiStart('Display', 'iter'); 

% Number of times I want to run optimization scheme
numstartpoints=1000;

% Manymins is a vector of solutions containing the distinct local minima found during the run;
%  runs MultiStart on numstartpoints to find a solution or multiple local solutions to problem
%[x,fval, exitflag, output, solutions]=run(ms,problem,numstartpoints);
[x,fval,exitflag,output,solutions]=run(ms,problem,numstartpoints);

alpha=x(1);
beta_A=0.00094; 
beta_P=0.00266;
theta_1=x(2);
epsilon=x(3);
mu=0.00868;  
mu_A=0.00775;   
mu_H=0.0271;
gamma=x(4);   
theta_2=3*x(2); 
sigma=x(5);
zeta=x(6);
theta_3=16*x(2);
nu=x(7);
omega=0.0000000001;

pars=[alpha,0.00094,0.00266,theta_1,epsilon,0.00868,0.00775,0.0271,gamma,theta_2,sigma,zeta,theta_3,nu,0.0000000001];

x
fval

% Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
N = 5; 
tspan=linspace(0,N,N+1);


% Initial conditions: assume for now 
% Later: although we know total number of prescription users in 2013, we do not
% know the initial number right at the start of 2013, so must be estimated;
% same with opioid addicts, heroin users, and stably recovered individuals.

S0=1-0.0538-0.0022-0.00074-0.000091;
P0=0.0538;
A0=0.0022;
H0=0.00074;
R0=0.000091;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];


[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);

 %Gives solution vector for each compartment 
  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);

  
 State1=y(:,1);
 State_data_1=[0.943169000000000;0.891413588268480;0.881203812830610;0.878597649275321;0.877380885777293;0.876427811124806];
 
 
 %Simulated data points for S and corresponding ODE solution plotted on top 
 figure(1)
 hold all
 plot(t,y(:,1))
 plot(t(1:end), State_data_1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Susceptibles')
 legend('Proportion of susceptibles simulated','Proportion of susceptibles (simulated) data' )

 State2=y(:,2);
 State_data_2=[0.0538000000000000;0.104768069960564;0.113970299917343;0.115542193082612;0.115733087680602;0.115675016908011];
 
 
 %Simulated data points for P and corresponding ODE solution plotted on top 
 figure(2)
 hold all
 plot(t,y(:,2))
 plot(t(1:end), State_data_2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Prescription Users')
 legend('Proportion of prescription users simulated','Proportion of prescription users (simulated) data' )


 State3=y(:,3);
 State_data_3=[0.00220000000000000;0.00260231845580068;0.00329371761143567;0.00402539259371014;0.00475061404614603;0.00546047944952132];
 
 
 %Simulated data points for A and corresponding ODE solution plotted on top 
 figure(3)
 hold all
 plot(t,y(:,3))
 plot(t(1:end), State_data_3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend('Proportion of opioid addicts simulated','Proportion of opioid addicts (simulated) data' )

 
 State4=y(:,4);
 State_data_4=[0.000740000000000000;0.000695209297496630;0.000692469877139593;0.000708223412451007;0.000735383907455046;0.000771298718040809];
 
 %Simulated data points for H and corresponding ODE solution plotted on top 
 figure(4)
 hold all
 plot(t,y(:,4))
 plot(t(1:end), State_data_4, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend('Proportion of heroin/fentanyl addicts simulated','Proportion of heroin/fentanyl addicts (simulated) data' )
 

 State5=y(:,5);
 State_data_5=[9.10000000000000e-05;0.000520814017756605;0.000839699755530817;0.00112654162441951;0.00140002857527287;0.00166539378469865];
 
  %Simulated data points for R and corresponding ODE solution plotted on top 
 figure(5)
 hold all
 plot(t,y(:,5))
 plot(t(1:end), State_data_5, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Stably recovered addicts')
 legend('Proportion of stably recovered addicts simulated','Proportion of stably recovered addicts (simulated) data' )


 
 State6=y(:,6);
 State_data_6=[0;0.182033383805590;0.359034822103230;0.534965567360877;0.710553170576798;0.885931758097526];
 
 %Simulated data for X and corresponding ODE solution plotted on top 
 figure(6)
 hold all
 plot(t,y(:,6))
 plot(t(1:end), State_data_6, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('X(t)')
 legend('X ODE solution','data simulated' )


 State7=y(:,7);
 State_data_7=[0;0.00102848309363892;0.00250312495050547;0.00420970261205163;0.00610431176381940;0.00817483062262488];
 
 %Simulated data for X and corresponding ODE solution plotted on top 
 figure(7)
 hold all
 plot(t,y(:,7))
 plot(t(1:end), State_data_7, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('L(t)')
 legend('L ODE solution','data simulated' )
  
 State8=y(:,8);
 State_data_8=[0;5.18916207835739e-05;0.000143072502492641;0.000253761633853627;0.000378815607538132;0.000516926137590890];
 
 %Simulated data for M and corresponding ODE solution plotted on top 
 figure(8)
 hold all
 plot(t,y(:,8))
 plot(t(1:end), State_data_8, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('M(t)')
 legend('M ODE solution','data simulated' )
 
 
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);
 Data1=[0.235833383805590;0.281769508258203;0.289901045174989;0.291129796298533;0.291111675201331];
 
 %Simulated data points from proportion that is in P at some point in the year and corresponding ODE solution plotted on top 
 figure(9)
 hold all
 plot(t(1:end-1),Estim1)
 plot(t(1:end-1), Data1, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in P at some point during the year')
 legend('Data points interested in', 'ODE solution')

 %Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 Estim2=y(2:3,3)+y(3:4,7)-y(2:3,7); 
 Data2=[0.00407696031266723;0.00500029527298183];
 
 %Simulated data points from proportion that is in A at some point in the year and corresponding ODE solution plotted on top 
 figure(10)
 hold all
% plot(t(1:end-1),Estim2)
% plot(t(1:end-1), Data2, 'x')
 plot(t(2:3),Estim2)
 plot(t(2:3), Data2, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in A at some point during the year')
 legend('Data points interested in', 'ODE solution')

 %Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);  
 Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8); 
 Data3=[0.000786390179205697;0.000803159008500578;0.000833277386135513];
 
 %Simulated data points from proportion that is in H at some point in the year and corresponding ODE solution plotted on top 
 figure(11)
 hold all
 %plot(t(1:end-1),Estim3)
 %plot(t(1:end-1), Data3, 'x')
 plot(t(2:4),Estim3)
 plot(t(2:4), Data3, 'x')
 set(gca, 'fontsize',10)
 xlabel('Year')
 ylabel('Proportion in H at some point during the year')
 legend('Data points interested in', 'ODE solution')
 
 

function value = HeroinModel_ODE45(z)

%Parameters
alpha=z(1);
beta_A=0.00094; 
beta_P=0.00266;
theta_1=z(2);
epsilon=z(3);
mu=0.00868;  
mu_A=0.00775;   
mu_H=0.0271;
gamma=z(4);   
theta_2=3*z(2); 
sigma=z(5);
zeta=z(6);
theta_3=16*z(2);
nu=z(7);
omega=0.0000000001;


pars=[alpha,0.00094,0.00266,theta_1,epsilon,0.00868,0.00775,0.0271,gamma,theta_2,sigma,zeta,theta_3,nu,0.0000000001];

% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-? where t=0 represents 2013
% and t=N represents ?; with spacing (T-0)/((N+1)-1)=1 between the points:
N = 5; 
tspan=linspace(0,N,N+1);

% Initial conditions: assumed for now; although we know total number of prescription users in 2013, we do not
% know the initial number right at the start of 2013, so must be estimated;
% same with opioid addicts, heroin users, and stably recovered individuals.
S0=1-0.0538-0.0022-0.00074-0.000091;
P0=0.0538;
A0=0.0022;
H0=0.00074;
R0=0.000091;
X0=0;
L0=0;
M0=0;
initials = [S0;P0;A0;H0;R0;X0;L0;M0];




[t,y]=ode15s(@HeroinModel,tspan,initials,[],pars);
  
 
 
 %COMPARING MODEL ESTIMATES TO DATA 
 value=0;
 %%%%%
 % In order to count the total number of individuals in P at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
  
 % To get the output from the model of the proportion of non-addicted
 % prescription opioid users in 2013 (first entry of Estim1), 
 % we take the total number of non-addiction prescription opioid users 
 % at the beginning of 2013 (P_0=IC, which is estimated because we do not
 % know the number at a particular instant at the beginning of the year)
 % and add on the number of individuals that enter the P class at any point during the year 2013, 
 % which comes fromintegrating ODE X'=dy(6) from t=0 to t=1; this gives
 % first value in Estim1. 
 
 % To get the output from the model of the proportion of non-addicted prescription opioid users in 2014, 2015
 % 2016, and 2017 (remaining entries of Estim1),
 % we take the initial number of non-addicted prescription opioid users 
 % in 2014: y(1,2), 2015: y(2,2), 2016: y(3,2), and 2017: y(4,2) and add the number of individuals that enter
 % the P class at any point during each of these years, which comes from
 % integrating ODE X'=dy(6) but just focusing in on these four specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 % for 2016, we have to subtract because integrating gives total number of new cases from t=0 to t=4, so have to 
 % subtract off the number from t=0 to t=3. 
 % for 2017, we have to subtract because integrating gives total number of new cases from t=0 to t=5, so have to 
 % subtract off the number from t=0 to t=4. 

 % Yearly output from the model as a proportion of the population in P at some point during the year for
 % 2013-final year, Estim1 is a column vector

 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  

 % Actual proportions of population (updated 3/12/19) that were non-addicted prescription opioid users at some point
 % during the year for 2013-2017 
 % (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 
 % Data1=[1799015/5517176 1778733/5559006 1771581/5602117 1719363/5651993 1595465/5708586];
 
 % Data simulated when testing codes 
 
 Data1=[0.235833383805590;0.281769508258203;0.289901045174989;0.291129796298533;0.291111675201331];
 
 % The difference between estimated value and data
 
 Diff1=Estim1-Data1; 


 %%%%%
 % In order to count the total number of individuals in A at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of opioid addicts in 2014 and 2015 (Estim2),
 % we take the initial number of opioid addicts in 2014, y(2,3), and 2015, y(3,3),
 % and add the number of individuals that enter
 % the A class at any point during the year 2014 or 2015, which comes from
 % integrating ODE L'=dy(9) but just focusing in on these two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 

 
 % Yearly output from the model as a proportion of population in A at some point during the year for
 % 2013-final year, Estim2 is a column vector
 
 %when testing all points with simulated data
 %Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 
 %only 2014 and 2015 data
 Estim2=y(2:3,3)+y(3:4,7)-y(2:3,7); 
 
 % Actual proportions of population (updated 3/12/19) that were opioid addicted individuals in
 % the population at some point during the year in 2014 and 2015 
 % (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 % Data2=[48000/5602117 42000/5651993];
 
 % Data simulated when testing codes 
 
 Data2=[0.00407696031266723;0.00500029527298183];
 
 % The difference between estimated value and data
 
 Diff2=Estim2-Data2; 


 

 %%%%%
 % In order to count the total number of individuals in H at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of heroin/fentanyl addicts in
 % 2014 and 2015 (Estim3), we take initial number of heroin/fentanyl addicts in 2014, y(2,4), 
 % and 2015, y(3,4),and add the number of individuals that enter the H class at any point
 % during the year 2014 or 2015, which comes from
 % integrating ODE M'=dy(10) but just focusing in on the two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1. 
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 

 % Yearly output from the model as a proportion of population in H at some point during the year for
 % 2013-final year, Data3 is a column vector 
 
 %when testing all points with simulated data 
 %Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);
 
 %for 2014 and 2015 only
 Estim3=y(2:4,4)+y(3:5,8)-y(2:4,8);  
 
 % Actual proportion (updated 3/12/19) of heroin addicted individuals in the population at some point during the year 
 % in 2014, 2015, and 2016
 % (total number of heroin addicted individuals in 2014, 2015, and 2016 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes 

 Data3=[0.000786390179205697;0.000803159008500578;0.000833277386135513];
 
 % The difference between estimated value and data
 
 Diff3=Estim3-Data3;

 
 %Comparing simulated data for susceptibles to output of model for
 %susceptibles 
 State1=y(:,1);
 State_data_1=[0.943169000000000;0.891413588268480;0.881203812830610;0.878597649275321;0.877380885777293;0.876427811124806];
 
 State_diff_1= State1-State_data_1;
 
 %Comparing simulated data for prescription users to output of model for
 %prescription users
 State2=y(:,2);
 State_data_2=[0.0538000000000000;0.104768069960564;0.113970299917343;0.115542193082612;0.115733087680602;0.115675016908011];
 
 State_diff_2=State2-State_data_2;
 
 
 %Comparing simulated data for opioid addicts to output of model for
 %opioid addicts
 State3=y(:,3);
 State_data_3=[0.00220000000000000;0.00260231845580068;0.00329371761143567;0.00402539259371014;0.00475061404614603;0.00546047944952132];
 
 State_diff_3=State3-State_data_3;
 
 %Comparing simulated data for heroin addicts to output of model for
 %heroin addicts
 State4=y(:,4);
 State_data_4=[0.000740000000000000;0.000695209297496630;0.000692469877139593;0.000708223412451007;0.000735383907455046;0.000771298718040809];
 
 
 State_diff_4=State4-State_data_4;
 
 
 %Comparing simulated data for stably recovered individuals to output of model for
 %stably recovered individuals 
 State5=y(:,5);
 State_data_5=[9.10000000000000e-05;0.000520814017756605;0.000839699755530817;0.00112654162441951;0.00140002857527287;0.00166539378469865];
 
 State_diff_5=State5-State_data_5;
 
 %Comparing simulated data for 6th ODE X and the model output 
 State6=y(:,6);
 State_data_6=[0;0.182033383805590;0.359034822103230;0.534965567360877;0.710553170576798;0.885931758097526];
 
 State_diff_6=State6-State_data_6;
 
 %Comparing simulated data for 7th ODE L and the model output 
 State7=y(:,7);
 State_data_7=[0;0.00102848309363892;0.00250312495050547;0.00420970261205163;0.00610431176381940;0.00817483062262488];
 
 State_diff_7=State7-State_data_7;
 
 %Comparing simulated data for 8th ODE M and the model output 
 State8=y(:,8);
 State_data_8=[0;5.18916207835739e-05;0.000143072502492641;0.000253761633853627;0.000378815607538132;0.000516926137590890];
 
 State_diff_8=State8-State_data_8;
 
 %%%%%
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sqrt(sum from 1 to N of (diff#)^2)
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; 
 % gives least squares percentage error so each piece weighted evenly)
 
 %value=norm(State_diff_1,2)+norm(State_diff_2,2)+norm(State_diff_3,2)+norm(State_diff_4,2)+norm(State_diff_5,2);
 
 %Run1
 %value = norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5);
 
 %Run2 and Run3
 %value = norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5)+norm(State_diff_6,2)./norm(State_data_6)+norm(State_diff_7,2)./norm(State_data_7)+norm(State_diff_8,2)./norm(State_data_8);
 
 %value=norm(Diff1,2)+norm(Diff2,2)+norm(Diff3,2);
 
 %Run4
 value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);
 
 %Run12
 %value = norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3);
 
 %value=sum(Diff1.^2)./sum(Data1.^2)+sum(Diff2.^2)./sum(Data2.^2)+sum(Diff3.^2)./sum(Data3.^2);
 
 %value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(State_diff_1,2)./norm(State_data_1)+norm(State_diff_2,2)./norm(State_data_2)+norm(State_diff_3,2)./norm(State_data_3)+norm(State_diff_4,2)./norm(State_data_4)+norm(State_diff_5,2)./norm(State_data_5);
 
 % Want value=f(x) to be small value when run MultiStart  
 
end
 

function f = HeroinModel(t,y,pars)
f=zeros(8,1);
f(1)=-pars(1)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=pars(1)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);


% X' ODE to calculate the number of new cases of prescription opioid use over time; i.e.
%individuals who enter the P class at any time from S (used in Estim1 in HeroinModel_ODE45.m) 
f(6) = pars(1)*y(1);

% L' ODE to calculate the number of new cases of opioid addiction over time;
%i.e. individuals who enter the A class at any time (used in Estim2 in
%HeroinModel_ODE45.m)
f(7) = pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2);

% M' ODE to calculate the number of new cases of heroin/fentanyl addiction
%over time; i.e. individuals who enter the H class at any time (used in
%Estim3 in HeroinModel_ODE45.m)
f(8) = pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15));
% Transpose column vector solution into row vector 
 % f=f';  
 

end








           
 

