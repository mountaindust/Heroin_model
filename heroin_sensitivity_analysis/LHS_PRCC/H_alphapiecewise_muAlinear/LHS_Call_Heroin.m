function s=LHS_Call_Heroin(xmin,xmean,xmax,xsd,nsample,distrib,logscale)
% s=latin_hs(xmean,xsd,nsample,nvar)
% LHS from normal distribution, no correlation
% method of Stein
% Stein, M. 1987. Large Sample Properties of Simulations Using Latin Hypercube Sampling. 
%                 Technometrics 29:143-151

if nsample==1
    s=xmean;
    return
end
if nargin<6
    logscale=1e5; % original value = 1e3
    distrib='unif';
end
if nargin<7
    logscale=1e7; % original value = 1e5
end

[sample,nvar]=size(xmean);
if distrib == 'norm'   % you only need to specify xmean & xsd
    rand('twister',sum(100*clock));
    ran=rand(nsample,nvar); % rand('twister',sum(100*clock))
    s=zeros(nsample,nvar);
    %method of Stein
    for j=1: nvar
        idx=randperm(nsample);
        P=(idx'-ran(:,j))/nsample;       % probability of the cdf
        %SQRT(2) * ERFINV(2*P-1)
        s(:,j) = xmean(j) + sqrt(2)*erfinv(2*P-1).* xsd(j); % this can be replaced by any inverse distribution function
        %s(:,j) = xmean(j) + ltqnorm(P).* xsd(j); % this can be replaced by any inverse distribution function
    end
end


if distrib == 'unif'  % you only need to specify xmin & xmax
    if xmin==0
        xmin=1e-300;
    end
    nvar=length(xmin);
    %rand('twister',sum(100*clock));
    ran=rand(nsample,nvar);
    s=zeros(nsample,nvar);
    for j=1: nvar
        idx=randperm(nsample);
        P =(idx'-ran(:,j))/nsample;
        xmax(j);
        xmin(j);
        xmax(j)/xmin(j);
        if (xmax(j)<1 & xmin(j)<1) || (xmax(j)>1 & xmin(j)>1)
            'SAME RANGE';
            if (xmax(j)/xmin(j))<logscale%1e3 %% It uses the log scale if the order of magnitude of [xmax-xmin] is bigger than 1e2,
                '<1e3: LINEAR SCALE';
                s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
            else
                '>=1e3: LOG SCALE';
                s(:,j) = log(xmin(j)) + P.*abs(abs(log(xmax(j)))-abs(log(xmin(j))));
                s(:,j) = exp(s(:,j));
            end
        else
            'e- to e+';
            if (xmax(j)/xmin(j))<logscale%1e3 %% It uses the log scale if the order of magnitude of [xmax-xmin] is bigger than 1e2,
                '<1e3: LINEAR SCALE';
                s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
            else
                '>=1e3: LOG SCALE';
                s(:,j) = log(xmin(j)) + P.*abs(log(xmax(j))-log(xmin(j)));
                s(:,j) = exp(s(:,j));
            end
        end
    end
end
%hist(s);
