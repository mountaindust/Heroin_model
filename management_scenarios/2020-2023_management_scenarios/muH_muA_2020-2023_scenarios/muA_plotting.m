t=[0:.25:10];

% disp(a(0,pars))
% disp(a(1,pars))
% disp(a(2,pars))
% disp(a(3,pars))
% disp(a(3.25,pars))
% disp(a(4,pars))
% disp(a(5,pars))
% disp(a(6,pars))
disp(muA(7))
% disp(a(8,pars))
% disp(a(9,pars))
% disp(a(10,pars))

disp(muA(10))
disp(muA25percent(10))
disp(muA50percent(10))


% disp(muA(0,pars)-muA(1,pars))
% disp(muA(1,pars)-muA(2,pars))
% disp(muA(2,pars)-muA(3,pars))
% disp(muA(3,pars)-muA(4,pars))
% disp(muA(4,pars)-muA(5,pars))
% disp(muA(5,pars)-muA(6,pars))
% disp(muA(6,pars)-muA(7,pars))
% disp(muA(7,pars)-muA(8,pars))
% disp(muA(8,pars)-muA(9,pars))
% disp(muA(9,pars)-muA(10,pars))

% disp(muA(3,pars)-muA(3.25,pars))
% disp(muA(3.25,pars)-muA(3.5,pars))
% disp(muA(3.5,pars)-muA(3.75,pars))
% disp(muA(3.75,pars)-muA(4,pars))

 muA_vec_original=[muA(0),muA(0.25),muA(0.5),muA(0.75),...
                     muA(1),muA(1.25),muA(1.5),muA(1.75),...
                     muA(2),muA(2.25),muA(2.5),muA(2.75),...
                     muA(3),muA(3.25),muA(3.5),muA(3.75),...
                     muA(4),muA(4.25),muA(4.5),muA(4.75),...
                     muA(5),muA(5.25),muA(5.5),muA(5.75),...
                     muA(6),muA(6.25),muA(6.5),muA(6.75),...
                     muA(7),muA(7.25),muA(7.5),muA(7.75),...
                     muA(8),muA(8.25),muA(8.5),muA(8.75),...
                     muA(9),muA(9.25),muA(9.5),muA(9.75),...
                     muA(10)];
                 
muA_vec_25percentcase=[muA25percent(0),muA25percent(0.25),muA25percent(0.5),muA25percent(0.75),...
                     muA25percent(1),muA25percent(1.25),muA25percent(1.5),muA25percent(1.75),...
                     muA25percent(2),muA25percent(2.25),muA25percent(2.5),muA25percent(2.75),...
                     muA25percent(3),muA25percent(3.25),muA25percent(3.5),muA25percent(3.75),...
                     muA25percent(4),muA25percent(4.25),muA25percent(4.5),muA25percent(4.75),...
                     muA25percent(5),muA25percent(5.25),muA25percent(5.5),muA25percent(5.75),...
                     muA25percent(6),muA25percent(6.25),muA25percent(6.5),muA25percent(6.75),...
                     muA25percent(7),muA25percent(7.25),muA25percent(7.5),muA25percent(7.75),...
                     muA25percent(8),muA25percent(8.25),muA25percent(8.5),muA25percent(8.75),...
                     muA25percent(9),muA25percent(9.25),muA25percent(9.5),muA25percent(9.75),...
                     muA25percent(10)];
                 
muA_vec_50percentcase=[muA50percent(0),muA50percent(0.25),muA25percent(0.5),muA25percent(0.75),...
                     muA50percent(1),muA50percent(1.25),muA50percent(1.5),muA50percent(1.75),...
                     muA50percent(2),muA50percent(2.25),muA50percent(2.5),muA50percent(2.75),...
                     muA50percent(3),muA50percent(3.25),muA50percent(3.5),muA50percent(3.75),...
                     muA50percent(4),muA50percent(4.25),muA50percent(4.5),muA50percent(4.75),...
                     muA50percent(5),muA50percent(5.25),muA50percent(5.5),muA50percent(5.75),...
                     muA50percent(6),muA50percent(6.25),muA50percent(6.5),muA50percent(6.75),...
                     muA50percent(7),muA50percent(7.25),muA50percent(7.5),muA50percent(7.75),...
                     muA50percent(8),muA50percent(8.25),muA50percent(8.5),muA50percent(8.75),...
                     muA50percent(9),muA50percent(9.25),muA50percent(9.5),muA50percent(9.75),...
                     muA50percent(10)];

 
 %to graph muA over the ten years 
 t1=0:0.25:10;
 t2=7:0.25:10;
 
 figure(1)
 hold all
 connect1=interp1(t,muA_vec_original,t1);
 plot(t1,connect1,'k-','LineWidth',3)
 connect2=interp1(t,muA_vec_25percentcase,t2);
 plot(t2,connect2,'g-','LineWidth',3)
 connect3=interp1(t,muA_vec_50percentcase,t2);
 plot(t2,connect3,'r-','LineWidth',3)
 set(gca, 'FontSize',16)
 xlabel('Year')
 ylabel('\mu_A')
 legend({'Baseline \muA', '\muA reduced 25%','\muA reduced 50%'},'FontSize',14, 'Location', 'northeast')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019','2020','2021','2022','2023'})


%baseline alpha in 2023, must shift t since function is formulated for t=0
%to t=3 and running from t=7 to t=10 for the alpha functions (2020-2023)


function muAoriginal = muA(t)
if  t<=7
    muAoriginal = 0.000977482526657751*t+0.00883138792481281;
else
    muAoriginal = 0.000977482526657751*7+0.00883138792481281+0.000977482526657751*(t-7);
end
end

%muA reduced 25 percent in 2023
function muA25 = muA25percent(t)
    muA25 = 0.000977482526657751*7+0.00883138792481281-5.730352392914993e-04*(t-7);
end

%muA reduced 50 percent in 2023
function muA50 = muA50percent(t)
    muA50 = 0.000977482526657751*7+0.00883138792481281-0.002123553005241*(t-7);
end