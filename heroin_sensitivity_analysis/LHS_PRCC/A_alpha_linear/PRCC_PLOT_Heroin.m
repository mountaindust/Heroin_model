%% Plot the residual of the partial regression of X (input - LHS matrix) and Y (output)
%% at column s (time points saved). PCC Coefficients are calculated on these
%% var: labels of the parameters varied in the X (as legend)
%% The Title of the plot is the Pearson correlation coefficient of the
%% transformed data, that is  the PRCC calculated on the original data.
%% The p-value is also showed in the title
%% by Simeone Marino, June 5 2007 %%

function PRCC_PLOT_Heroin(X,Y,s,PRCC_var,y_var)

Y=Y(s,:);
% Y;
[a k]=size(X); % Define the size of LHS matrix
Xranked=rankingN_Heroin(X);
Yranked=ranking_Heroin(Y);
for i=1:k  % Loop for the whole submatrices, Zi
    c1=['LHStemp=Xranked;LHStemp(:,',num2str(i),')=[];Z',num2str(i),'=[ones(a,1) LHStemp];LHStemp=[];'];
    eval(c1);
end
for i=1:k
    c2=['[b',num2str(i),',bint',num2str(i),',r',num2str(i),']= regress(Yranked,Z',num2str(i),');'];
    c3=['[b',num2str(i),',bint',num2str(i),',rx',num2str(i),']= regress(Xranked(:,',num2str(i),'),Z',num2str(i),');'];
    eval(c2);
    eval(c3);
end
% for i=1:k
%     c4=['r',num2str(i)];
%     c5=['rx',num2str(i)];
%     [r p]=corr(eval(c4),eval(c5));
%     a=['[PRCC , p-value] = ' '[' num2str(r) ' , '  num2str(p) '].'];% ' Time point=' num2str(s-1)];
%     figure,plot((eval(c4)),(eval(c5)),'.'),title(a),...
%             legend(PRCC_var{i}),xlabel(PRCC_var{i}),ylabel(y_var);%eval(c6);
% 
% end

AA(1,:)={'Variable'          'PRCC'              'p-value'};

figure(11);
set(gcf, 'Position',  [1, 1, 1700, 800])

for i=1:16
    c4=['r',num2str(i)];
    c5=['rx',num2str(i)];
    [r p]=corr(eval(c4),eval(c5));
    a=['[PRCC , p-value] = ' '[' num2str(r) ' , '  num2str(p) '].'];% ' Time point=' num2str(s-1)];
    %figure,
    subplot(4,4,i)
    plot((eval(c4)),(eval(c5)),'.'),title(a),...
            legend(PRCC_var{i}),xlabel(PRCC_var{i}),ylabel(y_var);
AA(i+1,:)={PRCC_var{i} num2str(r) num2str(p)};
end


%figure(12);
%set(gcf, 'Position',  [1, 1, 1300, 515])

%for i=10:16
  %  c4=['r',num2str(i)];
  %  c5=['rx',num2str(i)];
  %  [r p]=corr(eval(c4),eval(c5));
  %  a=['[PRCC , p-value] = ' '[' num2str(r) ' , '  num2str(p) '].'];% ' Time point=' num2str(s-1)];
  %  %figure,
  %  subplot(3,3,i-9)
  %  plot((eval(c4)),(eval(c5)),'.'),title(a),...
  %          legend(PRCC_var{i}),xlabel(PRCC_var{i}),ylabel(y_var);
%AA(i+1,:)={PRCC_var{i} num2str(r) num2str(p)};
%end

% figure(4);
% 
% for i=19:20
%     c4=['r',num2str(i)];
%     c5=['rx',num2str(i)];
%     [r p]=corr(eval(c4),eval(c5));
%     a=['[PRCC , p-value] = ' '[' num2str(r) ' , '  num2str(p) '].'];% ' Time point=' num2str(s-1)];
%     %figure,
%     subplot(3,3,i-18)
%     plot((eval(c4)),(eval(c5)),'.'),title(a),...
%             legend(PRCC_var{i}),xlabel(PRCC_var{i}),ylabel(y_var);
% AA(i+1,:)={PRCC_var{i} num2str(r) num2str(p)};
% end

% figure(5);
% 

%Shows PRCC results in command window 
AA




