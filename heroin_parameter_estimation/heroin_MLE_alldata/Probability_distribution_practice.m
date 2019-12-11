 %Probability distributions coded (FIGURE OUT IF NEED TO DO DISTRIBUTION OR
 %DENSITY FUNCTION) 

p=0.2;
H=0.002;
T=5000000;

%Distrib0 = makedist('Beta','alpha',p*(y(:,4)*T,'beta',(1-p)*(y(:,4)*T-1));
%Distrib0 = makedist('Beta','a',p*(H*T-1),'b',(1-p)*(H*T-1))

%Density function & plot
% X = 0:.001:1;
% y1 = betapdf(X,p*(H*T-1),(1-p)*(H*T-1));
% 
% figure
% plot(X,y1,'Color','r','LineWidth',2)
% hold on
% legend('a =p(HT-1), b=(1-p)(HT-1)','Location','NorthEast');
% hold off 
% 
% y2 = betapdf(0.099000000000000,p*(H*T-1),(1-p)*(H*T-1));

X = 0:.01:1;
y1 = betapdf(X,0.75,0.75);
y2 = betapdf(X,1,1);
y3 = betapdf(X,4,4);

figure
plot(X,y1,'Color','r','LineWidth',2)
hold on
plot(X,y2,'LineStyle','-.','Color','b','LineWidth',2)
plot(X,y3,'LineStyle',':','Color','g','LineWidth',2)
legend({'a = b = 0.75','a = b = 1','a = b = 4'},'Location','NorthEast');
hold off
 