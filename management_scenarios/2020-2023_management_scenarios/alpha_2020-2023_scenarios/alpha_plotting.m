m=-0.00559565027929907;
b=0.270110337915851;
c=-0.0269690987063522;
t=[0:.25:10];

pars=[m,b,c];

% disp(a(0,pars))
% disp(a(1,pars))
% disp(a(2,pars))
% disp(a(3,pars))
% disp(a(3.25,pars))
% disp(a(4,pars))
% disp(a(5,pars))
% disp(a(6,pars))
disp(a(7,pars))
% disp(a(8,pars))
% disp(a(9,pars))
% disp(a(10,pars))

disp(a(10,pars))
disp(a25percent(10))
disp(a50percent(10))


% disp(a(0,pars)-a(1,pars))
% disp(a(1,pars)-a(2,pars))
% disp(a(2,pars)-a(3,pars))
% disp(a(3,pars)-a(4,pars))
% disp(a(4,pars)-a(5,pars))
% disp(a(5,pars)-a(6,pars))
% disp(a(6,pars)-a(7,pars))
% disp(a(7,pars)-a(8,pars))
% disp(a(8,pars)-a(9,pars))
% disp(a(9,pars)-a(10,pars))

% disp(a(3,pars)-a(3.25,pars))
% disp(a(3.25,pars)-a(3.5,pars))
% disp(a(3.5,pars)-a(3.75,pars))
% disp(a(3.75,pars)-a(4,pars))

 alpha_vec_original=[a(0,pars),a(0.25,pars),a(0.5,pars),a(0.75,pars),...
                     a(1,pars),a(1.25,pars),a(1.5,pars),a(1.75,pars),...
                     a(2,pars),a(2.25,pars),a(2.5,pars),a(2.75,pars),...
                     a(3,pars),a(3.25,pars),a(3.5,pars),a(3.75,pars),...
                     a(4,pars),a(4.25,pars),a(4.5,pars),a(4.75,pars),...
                     a(5,pars),a(5.25,pars),a(5.5,pars),a(5.75,pars),...
                     a(6,pars),a(6.25,pars),a(6.5,pars),a(6.75,pars),...
                     a(7,pars),a(7.25,pars),a(7.5,pars),a(7.75,pars),...
                     a(8,pars),a(8.25,pars),a(8.5,pars),a(8.75,pars),...
                     a(9,pars),a(9.25,pars),a(9.5,pars),a(9.75,pars),...
                     a(10,pars)];
                 
alpha_vec_25percentcase=[a25percent(0),a25percent(0.25),a25percent(0.5),a25percent(0.75),...
                     a25percent(1),a25percent(1.25),a25percent(1.5),a25percent(1.75),...
                     a25percent(2),a25percent(2.25),a25percent(2.5),a25percent(2.75),...
                     a25percent(3),a25percent(3.25),a25percent(3.5),a25percent(3.75),...
                     a25percent(4),a25percent(4.25),a25percent(4.5),a25percent(4.75),...
                     a25percent(5),a25percent(5.25),a25percent(5.5),a25percent(5.75),...
                     a25percent(6),a25percent(6.25),a25percent(6.5),a25percent(6.75),...
                     a25percent(7),a25percent(7.25),a25percent(7.5),a25percent(7.75),...
                     a25percent(8),a25percent(8.25),a25percent(8.5),a25percent(8.75),...
                     a25percent(9),a25percent(9.25),a25percent(9.5),a25percent(9.75),...
                     a25percent(10)];
                 
alpha_vec_50percentcase=[a50percent(0),a50percent(0.25),a25percent(0.5),a25percent(0.75),...
                     a50percent(1),a50percent(1.25),a50percent(1.5),a50percent(1.75),...
                     a50percent(2),a50percent(2.25),a50percent(2.5),a50percent(2.75),...
                     a50percent(3),a50percent(3.25),a50percent(3.5),a50percent(3.75),...
                     a50percent(4),a50percent(4.25),a50percent(4.5),a50percent(4.75),...
                     a50percent(5),a50percent(5.25),a50percent(5.5),a50percent(5.75),...
                     a50percent(6),a50percent(6.25),a50percent(6.5),a50percent(6.75),...
                     a50percent(7),a50percent(7.25),a50percent(7.5),a50percent(7.75),...
                     a50percent(8),a50percent(8.25),a50percent(8.5),a50percent(8.75),...
                     a50percent(9),a50percent(9.25),a50percent(9.5),a50percent(9.75),...
                     a50percent(10)];

 
 t1=0:0.25:10;
 t2=7:0.25:10;
 
 figure(1)
 hold all
 connect1=interp1(t,alpha_vec_original,t1);
 plot(t1,connect1,'k-','LineWidth',3)
 connect2=interp1(t,alpha_vec_25percentcase,t2);
 plot(t2,connect2,'g-','LineWidth',3)
 connect3=interp1(t,alpha_vec_50percentcase,t2);
 plot(t2,connect3,'r-','LineWidth',3)
 set(gca, 'FontSize',16)
 xlabel('Year')
 ylabel('\alpha')
 legend({'Baseline \alpha', '\alpha reduced 25%','\alpha reduced 50%'},'FontSize',14, 'Location', 'northeast')
 set(gca, 'xtick', [ 0 1 2 3 4 5 6 7 8 9 10])
 set(gca,'xticklabel',{'2013', '2014', '2015', '2016', '2017','2018', '2019','2020','2021','2022','2023'})



function alpha = a(t,pars)
if  t<=3.25 
    alpha = pars(1)*t+pars(2);
else 
    alpha = pars(1)*3.25+pars(2)-pars(3)*3.25+pars(3)*t;
end
end

%alpha reduced 25 percent
function alpha25percent = a25percent(t)
    alpha25percent = 0.38031-0.0327885*t;
end

%alpha reduced 50 percent
function alpha50percent = a50percent(t)
    alpha50percent = 0.421085-0.0386135*t;
end