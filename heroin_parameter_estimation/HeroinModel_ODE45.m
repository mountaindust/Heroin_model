% File name: HeroinModel_ODE45.m (used to be in Heroin_model folder)
function value = HeroinModel_ODE45(z)


% Final time; don't want to run too long because dynamics can drastically change over
% a number of years, so here we will do 2013-? where t=0 represents 2013 and t=N represents ?:
N = 25; 
T = N;
% Generate N points, with spacing (T-0)/((N+1)-1)=1 between the points
tspan=linspace(0,T,N+1);
% global value
 

% Initial conditions: although we know total number of prescription users in 2013, we do not
% know the initial number right at the start of 2013, so must be estimated;
% same with opioid addicts, heroin users, and stably recovered individuals.

S0=1-z(12)-z(13)-z(14)-z(15); 
P0=z(12);
A0=z(13);
H0=z(14);
R0=z(15); 
X0=0;
L0=0;
M0=0;
initials = [S0,P0,A0,H0,R0,X0,L0,M0];


[t,y]=ode45(@(t,y) HeroinModel(t,y,z),tspan,initials);


  S=y(:,1);
  P=y(:,2);
  A=y(:,3);
  H=y(:,4);
  R=y(:,5);
  X=y(:,6);
  L=y(:,7);
  M=y(:,8);
  
 
 %COMPARING MODEL ESTIMATES TO DATA 
 
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

 
 % Yearly output from the model as a proportion of P individuals for 2013-2017 (individually typing in equations, same as below with for loop):
 %Estim1=[z(12)+y(1,6), y(1,2)+y(2,6)-y(1,6), y(2,2)+y(3,6)-y(2,6), y(3,2)+y(4,6)-y(3,6), y(4,2)+y(5,6)-y(4,6)];
 
 
 % OLD Yearly output from the model as a proportion of P individuals for
 % 2013-final year, Estim1 is a row vector
 %Estim1=zeros(1,25);
 % For 2013:
 %Estim1(1)=P0+y(2,6);  
 % For 2014-final year:
    %for i=2:25
      % Estim1(i)= y(i,2)+y(i+1,6)-y(i,6);
    %end
    
 Estim1=y(1:end-1,2)+y(2:end,6)-y(1:end-1,6);  

 % Actual proportions of population that were non-addicted prescription opioid users for
 % 2013-2017 (total number of non-addicted prescription opioid users in each year in TN that are 12 and older divided by
 % total population in TN 12 and older for each year) 
 % Data1=[1660630/5517176 1641908/5559006 1637623/5602117 1585227/5651993 1472737/5708586];
 
 % Data simulated when testing codes (Estim1 data)
 Data1=[0.231234288805174;0.279606982624230;0.288090450028217;0.289402440306315;0.289422210407325;0.289210210962975;0.288961043004265;0.288731787550339;0.288471095668373;0.288301855284399;0.288002971180094;0.287733156835231;0.287538990167825;0.287331376716734;0.287108092835211;0.286923223629713;0.286694831047265;0.286573230954602;0.286273673555622;0.286113828286653;0.285908962003143;0.285762816005408;0.285579458162195;0.285426249227171;0.285235803502895];
 
 % The difference between estimated value and data: 
 Diff1= Estim1-Data1;
 


 %%%%%
 % In order to count the total number of individuals in A at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of opioid addicts in 2014 and 2015 (Estim2),
 % we take the initial number of opioid addicts in 2014, y(1,3), and 2015, y(2,3),
 % and add the number of individuals that enter
 % the A class at any point during the year 2014 or 2015, which comes from
 % integrating ODE L'=dy(9) but just focusing in on these two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1;
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 

 
 % Yearly output from the model as a proportion of A individuals for 2014-2015 (individually typing in equations, same as below with for loop except that adds on beginning year): 
 %Estim2=[y(1,3)+y(2,7)-y(1,7), y(2,3)+y(3,7)-y(2,7)];
 
  
 % OLD Yearly output from the model as a proportion of A individuals for
 % 2013-final year, Estim2 is a row vector
% Estim2=zeros(1,25);  
 % For 2013:
 %Estim2(1)=A0+y(2,7);  
 % For 2014-final year:
   % for i=2:25
      % Estim2(i)= y(i,3)+y(i+1,7)-y(i,7);
    %end
 Estim2=y(1:end-1,3)+y(2:end,7)-y(1:end-1,7); 
 % Actual proportions of population that were opioid addicted individuals in
 % the population in 2014 and 2015 (total number of opioid addicted individuals in 2014 and 2015 in TN
 % that are 12 and older divided by the total population in TN 12 and older for each year) 
 % Data2=[42000/5651993 48000/5602117];
 
 % Data simulated when testing codes (Estim2 data)
 Data2=[0.00759310338821613;0.00797559589554223;0.00857761999932350;0.00922886976969081;0.00987804454828969;0.0105075880729010;0.0111104963696012;0.0116848964899262;0.0122284334360841;0.0127419564458735;0.0132244322194093;0.0136749129976324;0.0140937084194489;0.0144808479472825;0.0148357938333752;0.0151591436200827;0.0154498315164879;0.0157099417061239;0.0159379144327607;0.0161350350273053;0.0163013065376042;0.0164379295613826;0.0165453137891516;0.0166247937660384;0.0166769483056004];
 
 % The difference between estimated value and data
 Diff2=Estim2-Data2;
 

 %%%%%
 % In order to count the total number of individuals in H at some point throughout a certain year, 
 % we need to count the number who are in the class AT ALL during the year,
 % even if they leave or come back at some point. 
 
 % To get the output from the model of the proportion of heroin/fentanyl addicts in
 % 2014 and 2015 (Estim3), we take initial number of heroin/fentanyl addicts in 2014, y(2,3), 
 % and 2015, y(2,4),and add the number of individuals that enter the H class at any point
 % during the year 2014 or 2015, which comes from
 % integrating ODE M'=dy(10) but just focusing in on the two specific years:
 % for 2014, we have to subtract because integrating gives total number of new cases from t=0 to t=2, so have to 
 % subtract off the number from t=0 to t=1. 
 % for 2015, we have to subtract because integrating gives total number of new cases from t=0 to t=3, so have to 
 % subtract off the number from t=0 to t=2. 
 
 % Yearly output from the model as a proportion of P individuals for 2014-2016 (individually typing in equations, same as below with for loop except that adds on beginning year): 
 %Estim3=[y(1,4)+y(2,8)-y(1,8), y(2,4)+y(3,8)-y(2,8), y(3,4)+y(4,8)-y(3,8)];
 
 
 % Yearly output from the model as a proportion of H individuals for
 % 2013-final year, Estim3 is a row vector 
 %Estim3=zeros(1,25);  
 % For 2013:
% Estim3(1)=H0+y(2,8);  
 % For 2014-final year:
   % for i=2:25
    %   Estim3(i)= y(i,4)+y(i+1,8)-y(i,8);
  %  end
 Estim3=y(1:end-1,4)+y(2:end,8)-y(1:end-1,8);  
 
 % Actual proportion of heroin addicted individuals in the population in 2014 and 2015
 % Data3=[14000/5559006 14000/5602117 19000/5651993];
 
 % Data simulated when testing codes (Estim3 data)
 Data3=[0.00283624231246401;0.00286219806914482;0.00295926079526809;0.00309272264481946;0.00324987528094179;0.00342595791637886;0.00361910019025688;0.00382887564172068;0.00405497654121414;0.00429765240663133;0.00455701196578276;0.00483316661355280;0.00512632482769785;0.00543660739163244;0.00576402890032190;0.00610860924569894;0.00647013674817524;0.00684855459366240;0.00724336891745552;0.00765422492204435;0.00808049972466840;0.00852155477788744;0.00897653316390956;0.00944453916601961;0.00992446361946388];
 
 % The difference between estimated value and data
 Diff3=Estim3-Data3;
 
 
 %%%%%
 % STUDY MORE ABOUT 
 % The *relative* error that we are trying to minimize for ordinary least
 % squares: the sum of the squared errors (norm gives sqrt(sum from 1 to N of (diff#)^2)
 % normalized by norm of the data (because of difference in magnitude of
 % the data points in each estimation and the difference in the number of data points 
 % in each estimation, helpful to normalize; 
 % gives least squares percentage error so each piece weighted evenly)
 value = norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3);

 % Want value=f(x) to be small value when run MultiStart  