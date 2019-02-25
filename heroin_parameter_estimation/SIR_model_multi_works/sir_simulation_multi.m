%function sir_simulation_multi

    beta=1.0;
    gamma=1.0/10.0;    % five day infectious period	
    N=763;
    pars=[beta,gamma,N];
    
     
    tspan=linspace(0,13,14);	  
    y0=[760;3];       % take initial conditions to be S(0)=760, I(0)=3 
    
    
    [t,y]=ode45(@sir_rhs,tspan,y0,[],pars);
    
    
    plot(t,y(:,2));   % plot prevalence of infection (i.e. I(t)) over time

    hold on
    plot(t,y(:,1),'--');   % add susceptibles vs time to the plot
                           % as a dashed line
                           
    xlabel('Time')
    ylabel('State Variables')
    legend('I','S')
    
    hold off
    
%end

function f = sir_rhs(t,y,pars)
   f=zeros(2,1);
   f(1)=-pars(1)*y(1)*y(2)/pars(3);
   f(2)=pars(1)*y(1)*y(2)/pars(3)-pars(2)*y(2);
end
