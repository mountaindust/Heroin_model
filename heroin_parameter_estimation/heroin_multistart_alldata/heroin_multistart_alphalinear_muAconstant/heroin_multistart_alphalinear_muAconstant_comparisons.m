%File name: heroin_multistart_alphalinear_muAconstant_comparisons (i.e.
%with alpha piecewise, muA linear case)

%parameters for alpha linear, muA constant case
m=-0.0166568947289087;
beta_A=2.01339039713503e-05;
beta_P=1.07613488853937e-05;
theta_1=0.271953012780878;
epsilon=2.55879092337285;
mu=0.00710;  
mu_A=0.00884;
mu_H=0.0466;
gamma=0.00500069447307207;
theta_2=0.361860571741226;
sigma=0.100038203446277;
zeta=0.0974935199219447;
theta_3=19.9918788391698;
nu=0.000116349238936298;
omega=0.0000000001;
b=0.291928600613385;



pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];



% Final time N; will run from beginning of 2013 to beginning of 2019 where t=0 represents 2013
% and t=6 represents 2019, with spacing (N-0)/(25-1)=0.25 between the points to represent quarters of a year:
N = 6; 
%tspan=linspace(0,N,25);
% For smooth data plots; 73 marks to represent 72 months in 6 years 
tspan=linspace(0,N,73);


% Initial conditions
P0=0.0815718347996123;
A0=0.00762380594619159;
H0=0.000467529598012306;
R0=1.14233803732730e-05;
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
    
Data5=[351./5519417; 360./5559702; 377./5602187; 381./5648259];
 
Data6=[112./5519417; 201./5559702; 344./5602187; 488./5648259; 702./5702475];

Diff1=Estim1-Data1;
Diff2=Estim2-Data2;
Diff3=Estim3-Data3;
Diff4=Estim4-Data4;
Diff5=Estim5-Data5;
Diff6=Estim6-Data6;

value=norm(Diff1,2)./norm(Data1)+norm(Diff2,2)./norm(Data2)+norm(Diff3,2)./norm(Data3)+norm(Diff4,2)./norm(Data4)+norm(Diff5,2)./norm(Data5)+norm(Diff6,2)./norm(Data6)

 
 %these are more refined Estim1-6 solutions (more points in between)
 %for the alpha linear, muA constant case 
 
 for i=1:61;
    continuous1(i)=y(i,2)+y(i+12,6)-y(i,6);
 end 
 
 for i=1:61;
    continuous2(i)=y(i,3)+y(i+12,7)-y(i,7);
 end 
 
 for i=1:25;
    continuous3(i)=y(i+12,4)+y(i+24,8)-y(i+12,8);
 end 
 
 for i=1:70;
    continuous4(i)=y(i,2)+y(i+3,6)-y(i,6);
 end 
 
 for i=1:37;
    continuous5(i)=y(i+12,9)-y(i,9);
 end 
 
 for i=1:49;
     continuous6(i)=y(i+12,10)-y(i,10);
 end 
 
 %these are the Estim1-6 solutions for alpha piecewise, muA linear case to
 %be able to compare to continuous1-6 plots (these vectors COME FROM the
 %continuous1-6 values in the
 %heroin_multistart_alphapiecewise_muAlinear_data_fits.m code
 piecewiselinear1=[0.333601728399227,0.333175915483808,0.332728321407081,0.332261709006972,0.331780106556135,0.331285822903598,0.330781805718608,0.330270243323109,0.329752542479767,0.329234185659181,0.328712747924094,0.328186832641041,0.327656427978364,0.327121471723743,0.326586343618708,0.326051779683797,0.325515258078936,0.324976778804124,0.324436161331332,0.323894939641381,0.323354987516237,0.322814340991552,0.322272988827463,0.321730919784106,0.321188343663517,0.320589637051673,0.319939210399238,0.319245188880588,0.318501689710403,0.317702830103363,0.316836163635227,0.315846061713639,0.314753472916051,0.313548164931580,0.312219905449344,0.310758462158459,0.309166356568282,0.307475419903405,0.305677177199966,0.303768676095526,0.301752121333991,0.299626275354035,0.297397599721574,0.295118252226201,0.292765158692456,0.290344669208159,0.287865048100999,0.285333968692635,0.282785076012579,0.280228764371165,0.277640818987627,0.275021803620529,0.272372181989138,0.269720261939558,0.267082564225601,0.264434879447237,0.261776468070373,0.259106558128211,0.256435796049093,0.253770992362617,0.251103097887795];
 piecewiselinear2=[0.00809714168615252,0.00805057973988540,0.00800482150852405,0.00795983904182614,0.00791569162814685,0.00787236148442237,0.00782983069567537,0.00778798352559018,0.00774668914236077,0.00770595670642028,0.00766583177338046,0.00762634058082401,0.00758744371978373,0.00754893397305952,0.00751083953935850,0.00747313815723299,0.00743590610147142,0.00739914337207379,0.00736267577740545,0.00732649171262156,0.00729047382815929,0.00725471731508819,0.00721920608533321,0.00718392405081931,0.00714882372579705,0.00711367426551655,0.00707859108475582,0.00704356521282507,0.00700857542547382,0.00697360049845156,0.00693850158188103,0.00690319956770943,0.00686771669309500,0.00683201221024708,0.00679604537137500,0.00675977542868812,0.00672306352766428,0.00668588881302380,0.00664832246457599,0.00661035028308593,0.00657196335354125,0.00653299274776481,0.00649334262549700,0.00645317318807437,0.00641245674418663,0.00637118745783799,0.00632925096247715,0.00628656639948137,0.00624325172497510,0.00619932094031522,0.00615475741009546,0.00610947997533911,0.00606342187877861,0.00601668184885476,0.00596926962755058,0.00592117093903056,0.00587233393944851,0.00582270451131562,0.00577236951205417,0.00572135050278625,0.00566961883515297];
 piecewiselinear3=[0.00103648221384980,0.00107002455559030,0.00110458603756421,0.00114016260900581,0.00117676814912476,0.00121440265792105,0.00125345156247357,0.00129373094740695,0.00133523604175588,0.00137802594941656,0.00142213383531872,0.00146759286439210,0.00151442643833486,0.00156269159662780,0.00161245390926368,0.00166373576261723,0.00171657809983273,0.00177102186405450,0.00182706283366844,0.00188482567391737,0.00194431605587443,0.00200557577826225,0.00206864663980344,0.00213357043922064,0.00220053740736959];
 piecewiselinear4=[0.155076191543844,0.154944767134601,0.154790572135144,0.154616487508320,0.154426942963866,0.154224352146962,0.154012168655031,0.153792285735652,0.153565882568073,0.153338548968647,0.153107860000117,0.152872758089242,0.152633300844241,0.152389290881063,0.152145060138647,0.151901344637531,0.151655701778615,0.151408234836889,0.151158836466601,0.150908948138783,0.150660448526742,0.150411373074332,0.150161659569235,0.149911355814042,0.149660682850652,0.149410802944591,0.149160823036355,0.148910562597915,0.148660114139058,0.148409472768446,0.148158661330400,0.147907797470448,0.147656914416299,0.147406044030684,0.147098443736653,0.146739460146971,0.146337261538609,0.145906910046768,0.145440609334673,0.144922926531459,0.144300308903285,0.143595865303935,0.142807751366255,0.141974963403193,0.141064151825807,0.140084145128680,0.139052179703703,0.137968044226293,0.136867141912649,0.135759867073107,0.134618562621669,0.133441259077928,0.132230842505049,0.131016188273953,0.129813817138450,0.128601433938379,0.127380361907575,0.126148842859350,0.124916729162127,0.123690831345500,0.122462663987046,0.121232690806031,0.120000344277639,0.118768017219373,0.117537677841835,0.116306177807936,0.115072745150877,0.113838173930581,0.112603055836636,0.111368008507527];
 piecewiselinear5=[6.34945898043319e-05,6.36319293649670e-05,6.37705045675854e-05,6.39101042922182e-05,6.40503872193195e-05,6.41912249421352e-05,6.43324573385291e-05,6.44737038163033e-05,6.46145574546938e-05,6.47551056752634e-05,6.48954643897930e-05,6.50357003572356e-05,6.51755027366857e-05,6.53132370884030e-05,6.54499566245764e-05,6.55856269128759e-05,6.57203659293981e-05,6.58541736741431e-05,6.59840422264210e-05,6.61116918698183e-05,6.62373507637732e-05,6.63604649530502e-05,6.64808258657394e-05,6.65982249299312e-05,6.67127644146214e-05,6.68242401675997e-05,6.69313403379704e-05,6.70340899539851e-05,6.71321118689878e-05,6.72250289363227e-05,6.73151843786824e-05,6.73985620127609e-05,6.74760562826899e-05,6.75473071679625e-05,6.76119546480716e-05,6.76696387025106e-05,6.77163664109258e-05];
 piecewiselinear6=[2.63770784501621e-05,2.72171725816141e-05,2.80962751557056e-05,2.90157423083189e-05,2.99612075971519e-05,3.09307091152752e-05,3.19221180796356e-05,3.29492007698966e-05,3.40320569560434e-05,3.51330345182443e-05,3.62581791598768e-05,3.74109728809894e-05,3.85998050600412e-05,3.98687880815699e-05,4.11635405470877e-05,4.24805817026176e-05,4.38318377165285e-05,4.52173085888203e-05,4.66826993809339e-05,4.81917434168724e-05,4.97253088044604e-05,5.13095403229643e-05,5.29482938399488e-05,5.46454252229785e-05,5.63964009611590e-05,5.81721237924480e-05,6.00092999430002e-05,6.19101860339372e-05,6.38797920097807e-05,6.59231278150523e-05,6.80009250244965e-05,7.01487918321322e-05,7.23636340043446e-05,7.46475554800166e-05,7.70026601980311e-05,7.94310520972710e-05,8.19578456708497e-05,8.45758078674556e-05,8.72626591700114e-05,9.00161429573951e-05,9.28312492850853e-05,9.57535668134848e-05,9.88056083906411e-05,0.000101929156791526,0.000105124994469912,0.000108387161619352,0.000111762766473875,0.000115282532940518,0.000118868748845218];
 

 figure(1)
 hold all
 z1 = linspace(0,60,6); %defines mesh where going to plot Estim1, Data1 values 
 z7 = linspace(0,60,61);
 plot(z7,continuous1,'k-','LineWidth',1.3);
 plot(z7,piecewiselinear1,'r','LineWidth',1.3);
 %scatter(z1, Estim1, 100, 'o');
 scatter(z1, Data1, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in P')
 legend({'\alpha linear, \mu_A constant', '\alpha pw linear, \mu_A linear','Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 12 24 36 48 60])
 set(gca,'XLim',[0 60])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 

 
 % Data points from proportion that is in A at some point in the year and corresponding ODE solution points 
 figure(2)
 hold all
 z2 = linspace(0,60,6);
 z8 = linspace(0,60,61);
 plot(z8,continuous2,'k-','LineWidth',1.3);
 plot(z8,piecewiselinear2,'r-','LineWidth',1.3);
 %scatter(z2, Estim2, 100, 'o');
 scatter(z2, Data2, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in A')
 legend({'\alpha linear, \mu_A constant', '\alpha pw linear, \mu_A linear','Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 12 24 36 48 60])
 set(gca,'XLim',[0 60])
 set(gca,'xticklabel',{'2013','2014','2015','2016','2017','2018'})




 % Data points from proportion that is in H at some point in the year and corresponding ODE solution points 
 figure(3)
 hold all
 z3 = linspace(0,24,3);
 z9 = linspace(0,24,25);
 plot(z9,continuous3,'k-','LineWidth',1.3);
 plot(z9,piecewiselinear3,'r-','LineWidth',1.3);
 %scatter(z3, Estim3, 100,'o');
 scatter(z3, Data3, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in H')
 legend({'\alpha linear, \mu_A constant', '\alpha pw linear, \mu_A linear','Data'},'FontSize', 14, 'Location','northwest')
 set(gca, 'xtick', [ 0 12 24])
 set(gca,'XLim',[0 24])
 set(gca,'xticklabel',{'2014', '2015', '2016'})
 

 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(4)
 hold all
 z4 = linspace(0,69,24);
 z10 = linspace(0,69,70);
 plot(z10,continuous4,'k-','LineWidth',1.3);
 plot(z10,piecewiselinear4,'r-','LineWidth',1.3);
 %scatter(z4, Estim4, 100, 'o');
 scatter(z4, Data4, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Quarter')
 ylabel('Proportion in P')
 legend({'\alpha linear, \mu_A constant', '\alpha pw linear, \mu_A linear','Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69])
 xtickangle(90)
 set(gca,'XLim',[0 69])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
 
  
 figure(5)
 hold all
 z5 = linspace(0,36,4); %defines mesh where going to plot Estim5, Data5 values
 z11 = linspace(0,36,37);
 plot(z11,continuous5,'k-','LineWidth',1.3);
 plot(z11,piecewiselinear5,'r-','LineWidth',1.3);
 %scatter(z5, Estim5, 100,'o');
 scatter(z5, Data5, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion overdose from A') % at some point during the year
  legend({'\alpha linear, \mu_A constant', '\alpha pw linear, \mu_A linear','Data'},'FontSize', 14, 'Location','southwest')
 set(gca, 'xtick', [ 0 12 24 36])
 set(gca,'XLim',[0 36])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016'})
 
 
 figure(6)
 hold all
 z6 = linspace(0,48,5);
 z12 = linspace(0,48,49);
 plot(z12,continuous6,'k-','LineWidth',1.3);
 plot(z12,piecewiselinear6,'r-','LineWidth',1.3);
 %scatter(z6, Estim6, 100,'o');
 scatter(z6, Data6, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion overdose from H') % at some point during the year
 legend({'\alpha linear, \mu_A constant', '\alpha pw linear, \mu_A linear','Data'},'FontSize', 14,'Location','northwest')
 set(gca, 'xtick', [ 0 12 24 36 48])
 set(gca,'XLim',[0 48])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017'})  
 

 
  
 figure(7)
 hold all
 z1 = linspace(0,60,6); %defines mesh where going to plot Estim1, Data1 values 
 z7 = linspace(0,60,61);
 plot(z7,continuous1,'k-','LineWidth',1.3);
 %scatter(z1, Estim1, 100, 'o');
 scatter(z1, Data1, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Proportion in P')
 legend({'\alpha linear, \mu_A constant', 'Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 12 24 36 48 60])
 set(gca,'XLim',[0 60])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017', '2018'})
 
 
 % Data points from proportion that is in P at some point in the quarter of a year and corresponding ODE solution points 
 figure(8)
 hold all
 z4 = linspace(0,69,24);
 z10 = linspace(0,69,70);
 plot(z10,continuous4,'k-','LineWidth',1.3);
 %scatter(z4, Estim4, 100, 'o');
 scatter(z4, Data4, 90,'o','MarkerFaceColor',[0.01 0.28 1], 'MarkerEdgeColor',[0.01 0.28 1]);
 set(gca, 'fontsize',16)
 xlabel('Quarter')
 ylabel('Proportion in P')
 legend({'\alpha linear, \mu_A constant','Data'},'FontSize', 14)
 set(gca, 'xtick', [ 0 3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69])
 xtickangle(90)
 set(gca,'XLim',[0 69])
 set(gca,'xticklabel',{'Q1 2013', 'Q2 2013', 'Q3 2013', 'Q4 2013',...
                       'Q1 2014', 'Q2 2014', 'Q3 2014', 'Q4 2014',...
                       'Q1 2015', 'Q2 2015', 'Q3 2015', 'Q4 2015',...
                       'Q1 2016', 'Q2 2016', 'Q3 2016', 'Q4 2016',...
                       'Q1 2017', 'Q2 2017', 'Q3 2017', 'Q4 2017',...
                       'Q1 2018', 'Q2 2018', 'Q3 2018', 'Q4 2018'})
 
  

 
         
function alpha = a(t,pars)
alpha = pars(1)*t+pars(16);
end


function f = HeroinModel(t,y,pars)
f=zeros(10,1);
f(1)=-a(t,pars)*y(1)-pars(2)*y(1)*y(3)-pars(3)*y(1)*y(2)-pars(4)*y(1)*y(4)+pars(5)*y(2)+pars(6)*(y(2)+y(5))+(pars(6)+pars(7))*y(3)+(pars(6)+pars(8))*y(4);
f(2)=a(t,pars)*y(1)-pars(5)*y(2)-pars(9)*y(2)-pars(10)*y(2)*y(4)-pars(6)*y(2);
f(3)=pars(9)*y(2)+(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))+pars(2)*y(1)*y(3)+pars(3)*y(1)*y(2)-pars(12)*y(3)-pars(13)*y(3)*y(4)-pars(6)*y(3)-pars(7)*y(3);
f(4)=pars(4)*y(1)*y(4)+pars(10)*y(2)*y(4)+pars(13)*y(3)*y(4)+(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(14)*y(4)-(pars(6)+pars(8))*y(4);
f(5)=pars(12)*y(3)+pars(14)*y(4)-(pars(11)*y(5)*y(3))/(y(3)+y(4)+pars(15))-(pars(11)*y(5)*y(4))/(y(3)+y(4)+pars(15))-pars(6)*y(5);

% X' ODE to calculate the number of new cases of prescription opioid use over time;
% i.e. individuals who enter the P class at any time from S (used in
                                                             % Estim1, Estim4)
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



