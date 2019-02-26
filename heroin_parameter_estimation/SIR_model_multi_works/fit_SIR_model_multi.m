function fit_SIR_model_multi
    clf
    %for fmincon
    %A = [];
    %b = [];
    %Aeq = [];
    %beq = [];
    %lb = [0.1,0.1];
    %ub = [1,2];
    %[theta_hat,ess]=fmincon(@error_sum_of_squares,[0.9,0.2],A,b,Aeq,beq,lb,ub)
    
    problem = createOptimProblem('fmincon','objective',...
    @error_sum_of_squares,'x0',[0.9,0.2],'lb',[0.1,0.1],'ub',[2,2]);
    
    ms=MultiStart('Display', 'iter');
    [theta_hat,f] = run(ms,problem,2);

    beta=theta_hat(1);
    gamma=theta_hat(2);
    
    % integrate ODE for best fitting parameter values, so we can plot it
    pars=[beta,gamma,763];
    
    tspan=[0,13];	  
    y0=[760;3];       % take initial conditions to be S(0)=760, I(0)=3 
    
    theta_hat
    f
    
    [t,y]=ode45(@sir_rhs,tspan,y0,[],pars);

    figure(1)
    plot(t,y(:,2))
    
    hold on
    data=[3;7.33011471138265;17.7299719159973;41.8670611177493;93.7018028080027;188.378209579639;317.714565799195;434.750912802697;498.452764354942;510.680101141760;492.520343356411;460.788974154104;424.603096448541;388.218142479646];
    
    data_times=[0:13]; % data points are at t=0, 1, ... , 13
    plot(data_times,data,'x')
    hold off
    
end

function ESS=error_sum_of_squares(input_pars)

    beta=input_pars(1);
    gamma=input_pars(2);
    N=763;
    
    pars=[beta,gamma,N];
    
    tspan=[0:13];	  % 14 days of data, including initial value
    y0=[760;3];       % take initial conditions to be S(0)=760, I(0)=3 
    
    [t,y]=ode45(@sir_rhs,tspan,y0,[],pars);

    ESS=0;
    data=[3;7.33011471138265;17.7299719159973;41.8670611177493;93.7018028080027;188.378209579639;317.714565799195;434.750912802697;498.452764354942;510.680101141760;492.520343356411;460.788974154104;424.603096448541;388.218142479646];
    
    % initial value plus 13 more values
    
    diff=data-y(:,2);  % calculate differences between data and predictions

    ESS=sum(diff.^2);  % square entries of diff and then sum
    % note: there are 14 terms in this sum, but one (the first) is zero
end

function f = sir_rhs(t,y,pars)
f=zeros(2,1);
f(1)=-pars(1)*y(1)*y(2)/pars(3);
f(2)=pars(1)*y(1)*y(2)/pars(3)-pars(2)*y(2);
end

%Order: MultiStart starts error_sum_of_squares with starting points
    %(x0=[input_pars(1), input_pars(2)], those values are put into pars
    %vector so that ode45 can be evaluated with those parameters, then the
    %model is compared with the data. This is done over and over again
    %until the optimal solution is found, which is called theta_hat, which
    %then is put into pars vector and used to plot optimal solution with
    %the data in MultiStart function. 