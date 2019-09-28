m=-0.016168770430769;
beta_A=6.028626431922449e-05;
beta_P=1.377512810861206e-05;
theta_1=0.140625525950935;
epsilon=2.526649250596874;
mu=0.00710;  
mu_A=0.00884;
mu_H=0.0466;
gamma=0.005003442118174;
theta_2=1.823588391956158;
sigma=0.100187240737003;
zeta=0.099396517965999;
theta_3=18.957642027116236;
nu=1.792181772640938e-04;
omega=0.0000000001;
b=0.289029299448868;





pars=[m,beta_A,beta_P,theta_1,epsilon,mu,mu_A,mu_H,gamma,theta_2,sigma,zeta,theta_3,nu,omega,b];



% Final time N; will run from beginning of 2013 to beginning of 2019 where t=0 represents 2013
% and t=6 represents 2019, with spacing (N-0)/(25-1)=0.25 between the points to represent quarters of a year:
N = 6; 
%tspan=linspace(0,N,25);
% For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);
% For smooth data plots; 73 marks to represent 72 months in 6 years 
tspan=linspace(0,N,73);


% Initial conditions
P0=0.083311348758638;
A0=0.007621840599152;
H0=4.551016351747697e-04;
R0=1.696410024717003e-05;
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
 
 piecewiselinear1=[0.333576205455922,0.332948533551826,0.332346393294431,0.331769817316041,0.331209388450098,0.330668643498609,0.330136986368282,0.329618400373071,0.329104773977755,0.328598549244399,0.328092468421844,0.327589528044105,0.327091061283114,0.326596761866611,0.326106841540710,0.325614284830568,0.325122824872527,0.324633063612023,0.324144886729343,0.323658328822512,0.323169487396035,0.322680851070402,0.322192446693094,0.321704256633001,0.321216194621781,0.320653899778063,0.320050803512866,0.319404143473015,0.318708217549912,0.317957368226896,0.317129892791311,0.316175418558624,0.315114364775579,0.313936047667188,0.312629783458462,0.311184888374410,0.309605734907321,0.307928425632358,0.306135261827283,0.304227279416610,0.302206652078159,0.300070712048836,0.297834134894508,0.295540542776882,0.293169384402421,0.290728998217721,0.288228575901356,0.285672189647857,0.283103868672025,0.280523292052817,0.277911510740975,0.275267169828702,0.272590269315998,0.269912974803048,0.267244614485013,0.264566831930624,0.261879641438743,0.259189250634676,0.256503403170428,0.253814183711539,0.251121845547804];
 piecewiselinear2=[0.00811575443808113,0.00806398225188146,0.00801337030573588,0.00796374035197106,0.00791515262163700,0.00786763643794135,0.00782113011704870,0.00777565684732437,0.00773090676354465,0.00768705471542833,0.00764396524809173,0.00760169423486230,0.00756025210126816,0.00751940701839346,0.00747931926478225,0.00743976579236534,0.00740086529144836,0.00736261654939566,0.00732487188085702,0.00728769807465693,0.00725082816586943,0.00721440645096945,0.00717841881400149,0.00714284827686090,0.00710766447525176,0.00707257378376235,0.00703773409306511,0.00700312591791188,0.00696872708775941,0.00693450112659798,0.00690029067990356,0.00686603338161363,0.00683174158473179,0.00679737030434926,0.00676287455555730,0.00672820935344713,0.00669320080044134,0.00665787140615845,0.00662226891993104,0.00658638072093730,0.00655019532206156,0.00651352940226745,0.00647632588986118,0.00643870332356205,0.00640063352654956,0.00636204253642290,0.00632282624997939,0.00628306141208440,0.00624276361790246,0.00620191638701243,0.00616039067030579,0.00611826917999747,0.00607555191608746,0.00603223953100877,0.00598820195976707,0.00594351638173859,0.00589818443081179,0.00585221555434996,0.00580556349929877,0.00575822195626523,0.00571020249411736];
 piecewiselinear3=[0.00103879401872115,0.00107222182003985,0.00110670060664646,0.00114219319608747,0.00117871937421794,0.00121632953994080,0.00125535127229703,0.00129562571008728,0.00133713526000372,0.00137994859796838,0.00142410117014150,0.00146962782229780,0.00151655974639220,0.00156497390674239,0.00161490037194226,0.00166637771004492,0.00171944847587797,0.00177414875288072,0.00183049085594584,0.00188857924556766,0.00194842829734733,0.00201008033134992,0.00207357766764055,0.00213896262628429,0.00220644495011542];
 piecewiselinear4=[0.155325253060495,0.154951319739505,0.154604885754248,0.154285343789852,0.153983383152519,0.153701520564260,0.153429418460399,0.153170194163698,0.152916692001478,0.152671129891969,0.152426243238392,0.152185028574759,0.151947980152408,0.151715554787693,0.151487827922813,0.151257784082925,0.151029156404370,0.150801963493775,0.150576591083820,0.150353046020125,0.150127443438303,0.149902270041066,0.149677264034741,0.149452637165443,0.149228320795938,0.149003239587273,0.148778122418152,0.148552854968860,0.148327471837424,0.148101999984665,0.147876619006155,0.147651035318241,0.147425204943438,0.147199083904257,0.146899360311901,0.146558695778620,0.146174301604860,0.145768805711118,0.145323821823110,0.144818485926016,0.144204874965664,0.143505608213002,0.142726279454282,0.141896988312354,0.140986526065068,0.140005771898257,0.138968978705157,0.137877612276385,0.136775701826741,0.135662539883260,0.134512388287611,0.133323671390797,0.132099254303067,0.130871302624606,0.129652142763052,0.128425450809613,0.127188818611654,0.125948453794483,0.124713801052530,0.123478762785237,0.122241984084810,0.121004023067905,0.119767691843277,0.118529894036457,0.117291308456376,0.116051949401896,0.114812475541607,0.113571180222034,0.112328185116659,0.111083743515275];
 piecewiselinear5=[6.36372881196058e-05,6.37482877954494e-05,6.38608633994359e-05,6.39763126383315e-05,6.40934140142250e-05,6.42122019347613e-05,6.43320683781718e-05,6.44532419137534e-05,6.45760377636838e-05,6.47000428540358e-05,6.48243702479074e-05,6.49493764495990e-05,6.50752144134311e-05,6.52017590874050e-05,6.53290969281743e-05,6.54558131554005e-05,6.55826606144742e-05,6.57095579812112e-05,6.58352102098388e-05,6.59604602422571e-05,6.60843469258181e-05,6.62073564963982e-05,6.63294176750341e-05,6.64504491097486e-05,6.65703774304419e-05,6.66873007906452e-05,6.68019198416603e-05,6.69140803485210e-05,6.70235500143820e-05,6.71303021973156e-05,6.72347754298907e-05,6.73341836362378e-05,6.74293835263663e-05,6.75200331920707e-05,6.76057907251452e-05,6.76863142173843e-05,6.77571157328129e-05];
 piecewiselinear6=[2.65244585561514e-05,2.73638561837871e-05,2.82233159234496e-05,2.91330060649173e-05,3.00715282125241e-05,3.10378157396937e-05,3.20251846860443e-05,3.30361477052409e-05,3.41074557069078e-05,3.52128353476776e-05,3.63353501103575e-05,3.74816268551979e-05,3.86576988290718e-05,3.99112202157973e-05,4.12092445506439e-05,4.25208185767617e-05,4.38624134141611e-05,4.52421364645530e-05,4.66972764088844e-05,4.82095027836056e-05,4.97390633242339e-05,5.13180222647876e-05,5.29507810777232e-05,5.46412414328675e-05,5.63902633066851e-05,5.81688664929883e-05,6.00073050112591e-05,6.19100339233687e-05,6.38821010059406e-05,6.59231668025444e-05,6.80041706422048e-05,7.01572070718711e-05,7.23782884372618e-05,7.46695337577950e-05,7.70330620528890e-05,7.94709923419618e-05,8.20135890914223e-05,8.46409410782949e-05,8.73397819818060e-05,9.01056567400845e-05,9.29335175765068e-05,9.58676504805741e-05,9.89158831374022e-05,0.000102035256295998,0.000105225855940816,0.000108511756943363,0.000111930853938665,0.000115438121487446,0.000118999393453470];
 

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



