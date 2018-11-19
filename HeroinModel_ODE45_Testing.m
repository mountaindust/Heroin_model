function value = HeroinModel_ODE45_Testing(z)

%%%GOING TO CHANGE FOR MY MODEL (right now it's a different model) 

% Final time 
N = 180;
T=N;
tspan=linspace(0,T,N+1);
global value 
 

% Estimated  values of parameters from the "Model_MultiStar.m"
%z0=[0.598865461058500,0.127071357078022,0.116748546762720,9.99604205662634e-05,4.78764501980746e-05,1.35331547234463e-05,1.19563837848761e-07,4.12372760163290e-06,1.50585703330823e-05,0.519590745645255,0.574934579070339];
%z0= [0.523772604010212,0.104033132979486,0.190273501576715,9.99627115279089e-05,3.65423094912493e-05,1.34581576995548e-05,1.00221316942663e-07,5.74470493795167e-06,1.30656891116669e-05,0.555312086427258,0.494415504997107];
%z0=[0.189627932102808,0.156249066705754,0.153469283449332,9.19240019288754e-05,8.44645615142794e-05,1.30317350951584e-05,7.73975660455356e-06,5.21061060135415e-06,4.09792383058129e-06,0.569060145549251,0.584974060214357];
%z0=[0.133122315522406,0.186163341466148,0.554449932048055,7.95538371925574e-05,5.43545342741910e-05,1.51685046126194e-05,4.99067412766050e-06,7.09722489103569e-06,2.69569019387619e-05,0.144429593587413,0.674200088915047];
z0=[0.406316891580689,0.101200921126464,0.162915225906188,1.07283038812341e-06,1.83877199370086e-05,1.00698823176837e-07,2.56689031764850e-05,1.00236694289975e-07,9.99506682425477e-05,0.515435842166515,0.380148738139410];

z=z0;

% Parameter

% % Intrinsic Growth rate of Anchovy Population
r1=z(1); 
% % Intrinsic Growth rate of Jelly Fish Population
r2=z(2); 
% % Intrinsic Growth rate of Zoo-Plankton Population
r3=z(3);

% % Carriying Capacity of Anchovy population
K1=3e+5;
% % Carriying Capacity of Jelly Fish population
K2=1e+4;
% % Carriying Capacity of Zoo-Plankton population
K3=4e+5;
%Initials
A0=2.0839e5;
J0=8.2297e3;
Z0=2.983e5;
initials = [A0, J0, Z0];

%  Interaction rates:

% % Predation constant of Anchovy on its prey, Zoo-Plankton
m0=z(4);   
% % Consumption constant of Anchovy by its predator,Jelly Fish
m1=z(5);
% % Predation constant of Jelly Fish on its prey, Anchovy
m2=z(6);   
% % Predation constant of Jelly Fish on its prey, Zoo-Plankton
m3=z(7);
% % Consumption constant of Zoo-Plankton by its predator, Anchovy
m4=z(8);   
% % Consumption constant of Zoo-Plankton by its predator, Jelly Fish
m5=z(9);
% % The rate of decay of Jelly Fish by ts pradator, Beroe ovata, 
m6=z(10);
% % Constant Harvesting rate 
h=z(11);




[t,y]=ode45(@(t,y) Model(t,y,z),tspan, initials);


  A=y(:,1);
  J=y(:,2);
  Z=y(:,3);
  
  
  % int of harvesting , h*A

 for i=1:N
 temp(i) =h*(((y(i+1,1)+y(i,1))/2)*(t(i+1)-t(i)));
 end
 
 % yearly landing from the model
 
 Estim=[sum(temp(1:12)),sum(temp(13:24)),sum(temp(25:36)), sum(temp(37:48)),sum(temp(49:60)),...
        sum(temp(61:72)),sum(temp(73:84)),sum(temp(85:96)),sum(temp(97:108)),sum(temp(109:120)),...
        sum(temp(121:132)),sum(temp(133:144)),sum(temp(145:156)),sum(temp(157:168)),sum(temp(169:180))];
       % sum(temp(181:192)),sum(temp(193:204)),sum(temp(205:216)),sum(temp(217:228)),sum(temp(229:240))];

 
% yearly landing of Anchovy from 1997-2016 (20 years)
% Data=[213780 195996 310801 260670 288616 336419 266069 306656 119255 212081 357089 225344 185606 203026 246390 109187 255309  71530 195350 112500];

% Yearly landing of Anchovy from 2002-2016 (15 years)
 Data=[ 336419 266069 306656 119255 212081 357089 225344 185606 203026 246390 109187 255309  71530 195350 112500];
 
% Yearly landing of Anchovy from 2002-2016 (12 years)
% Data=[119255 212081 357089 225344 185606 203026 246390 109187 255309  71530 195350 112500];
 
 
 
% the difference between estimated value and data 
 diff = (Estim-Data);
 
 value = norm(diff,2)./norm(Data) ; 
 
 
 figure(1)
 hold all
 scatter(1:1:15,Estim,'filled')
 scatter(1:1:15, Data, 'filled')
 set(gca, 'xtick', [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ])
 set(gca, 'fontsize',10)
 set(gca,'xticklabel',{'2002','2003','2004','2005','2006','2007','2008', '2009','2010', '2011', '2012','2013','2014','2015','2016'})
 ylabel('Landing')
 legend('Landing Simulated','Yearly Landing Data' )
 
 figure(2)

           subplot(2,2,1);plot(t,y(:,1),'b-','LineWidth',1)
           subplot(2,2,1);xlabel('Time')
           subplot(2,2,1);ylabel('Anchovy')
           xlim([0 , T])
           
           subplot(2,2,2);plot(t,y(:,2),'r-','LineWidth',1)
           subplot(2,2,2);xlabel('Time')
           subplot(2,2,2);ylabel('Jelly Fish')
           xlim([0 , T])
           
           subplot(2,2,3);plot(t,y(:,3) ,' g-','LineWidth',1)
           subplot(2,2,3);xlabel('Time')
           subplot(2,2,3);ylabel('Zoo-Plankton')
           xlim([0 , T])
           %ylim([0 , 1e+2])
           
           subplot(2,2,4);
           plot(t,y(:,1),'b-','LineWidth',1);
           hold on
           plot(t,y(:,2),'r-','LineWidth',1);
           hold on
           plot(t,y(:,3),' g-','LineWidth',1); 
           xlabel('time')
           ylabel('Size of Populations');
           legend('A','P','Z')
           xlim([0 , T])
           ylim([0 , 2e+5])
           legend('A','P','Z')
 