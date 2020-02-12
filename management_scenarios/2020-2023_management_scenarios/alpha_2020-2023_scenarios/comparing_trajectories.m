% Final time and last entry of tspan is # of equally spaced points from 0 to N (quarterly linspace)
t=[1:0.25:4];
%For smooth plots (ONLY GOOD FOR ODE SOLUTIONS, NOT DATA/ESTIM PLOTS)
%tspan=linspace(0,N,3000);


 A_baseline_alpha_case=[0.00370000000000000;0.00353278337029652;0.00336289307560543;0.00319097482980213;0.00301766402796066;0.00284378632891392;0.00267026100551299;0.00249770644297795;0.00232715473139335;0.00215947534975907;0.00199553547623715;0.00183611325313993;0.00168192319325566];
 H_baseline_alpha_case=[0.00597000000000000;0.00646253334841020;0.00698258162521298;0.00753001934779802;0.00810629377051649;0.00871221154764143;0.00934807558313681;0.0100140013779746;0.0107101117664536;0.0114363784036355;0.0121928145645927;0.0129794511679071;0.0137963409589377];
 
 A_25percent_alpha_case=[0.00370000000000000;0.00353267040767489;0.00336238928474401;0.00318968332757997;0.00301514719584218;0.00283955422999254;0.00266379884441428;0.00248877771301105;0.00231554824945441;0.00214487141341584;0.00197761384549922;0.00181450857447343;0.00165638654016796];
 H_25percent_alpha_case=[0.00597000000000000;0.00646252562871979;0.00698213581555402;0.00752931844143701;0.00810535995274669;0.00871100819807622;0.00934652875810367;0.0100120066686606;0.0107074079365409;0.0114327074419068;0.0121878730733820;0.0129728826579485;0.0137877434894481];
 
 A_50percent_alpha_case=[0.00370000000000000;0.00353280392347552;0.00336233468380589;0.00318892577674831;0.00301315616123394;0.00283579021459923;0.00265775076305633;0.00248000877817790;0.00230356335342014;0.00212938188714762;0.00195838905036281;0.00179145596662234;0.00162939725352773];
 H_50percent_alpha_case=[0.00597000000000000;0.00646073977077021;0.00697854297962477;0.00752471520320035;0.00810023226051878;0.00870540350069781;0.00934038930763936;0.0100051708090444;0.0106996910829246;0.0114238593390079;0.0121775733566048;0.0129607303266466;0.0137732386906761];
 
  
 %compare all 3 alpha cases, output of A 
 figure(1)
 hold all
 connect1=interp1(t,A_baseline_alpha_case,t);
 plot(t,connect1,'k-','LineWidth',3)
 connect2=interp1(t,A_25percent_alpha_case,t);
 plot(t,connect2,'g-','LineWidth',3)
 connect3=interp1(t,A_50percent_alpha_case,t);
 plot(t,connect3,'r-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Opioid addicts')
 legend({'Baseline \alpha', '\alpha reduced 25%','\alpha reduced 50%'},'FontSize',14, 'Location', 'northeast')
 set(gca, 'xtick', [ 1 2 3 4 ])
 set(gca,'xticklabel',{'2020', '2021', '2022', '2023'})   
            
 %compare all 3 alpha cases, output of H
 figure(2)
 hold all
 connect4=interp1(t,H_baseline_alpha_case,t);
 plot(t,connect4,'k-','LineWidth',3)
 connect5=interp1(t,H_25percent_alpha_case,t);
 plot(t,connect5,'g-','LineWidth',3)
 connect6=interp1(t,H_50percent_alpha_case,t);
 plot(t,connect6,'r-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend({'Baseline \alpha', '\alpha reduced 25%','\alpha reduced 50%'},'FontSize',14, 'Location', 'northwest')
 set(gca, 'xtick', [ 1 2 3 4 ])
 set(gca,'xticklabel',{'2020', '2021', '2022', '2023'})      
 
 
 %{
 %ONCE ALPHA CASE IS FIGURED OUT, CAN DO CODE FOR OVERDOSE DEATHS THAT
 OCCUR FROM THE VARIOUS ALPHA CASES AND PUT VALUES IN THESE VECTORS AND
 THEN PLOT
 Aoverdose_baseline_alpha_case=[];
 Aoverdose_25percent_alpha_case=[];
 Aoverdose_50percent_alpha_case=[];
 Hoverdose_baseline_alpha_case=[];
 Hoverdose_25percent_alpha_case=[];
 Hoverdose_50percent_alpha_case=[];
 
 
 %compare all 3 alpha cases, output of A overdoses
 figure(3)
 hold all
 connect7=interp1(t,Aoverdose_baseline_alpha_case,t);
 plot(t,connect7,'k-','LineWidth',3)
 connect8=interp1(t,Aoverdose_25percent_alpha_case,t);
 plot(t,connect8,'g-','LineWidth',3)
 connect9=interp1(t,Aoverdose_50percent_alpha_case,t);
 plot(t,connect9,'r-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend({'Baseline \alpha', '\alpha reduced 25%','\alpha reduced 50%'},'FontSize',14, 'Location', 'northwest')
 set(gca, 'xtick', [ 1 2 3 4 ])
 set(gca,'xticklabel',{'2020', '2021', '2022', '2023'})      
 
 %compare all 3 alpha cases, output of H overdoses
 figure(4)
 hold all
 connect10=interp1(t,Hoverdose_baseline_alpha_case,t);
 plot(t,connect10,'k-','LineWidth',3)
 connect11=interp1(t,Hoverdose_25percent_alpha_case,t);
 plot(t,connect11,'g-','LineWidth',3)
 connect12=interp1(t,Hoverdose_50percent_alpha_case,t);
 plot(t,connect12,'r-','LineWidth',3)
 set(gca, 'fontsize',16)
 xlabel('Year')
 ylabel('Heroin/fentanyl addicts')
 legend({'Baseline \alpha', '\alpha reduced 25%','\alpha reduced 50%'},'FontSize',14, 'Location', 'northwest')
 set(gca, 'xtick', [ 1 2 3 4 ])
 set(gca,'xticklabel',{'2020', '2021', '2022', '2023'})      
 %}
 
      
 
 