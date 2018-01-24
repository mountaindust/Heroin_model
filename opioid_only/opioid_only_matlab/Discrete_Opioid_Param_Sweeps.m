%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code solves a model of a basic opioid addiction epidemic using
% discrete dynamics
%
%
% Author: Nick Battista
% Date: August 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Discrete_Opioid_Param_Sweeps()


%
% Temporal information
%
Tstop = 1000;
t=0:1:Tstop;


%
% Dynamical Coupling Parameters
%

da = 0.01; % <--- Shows up as x-axis in VisIt
dz = 0.01; % <--- Shows up as y-axis in VisIt

alphaVec = 0:da:1;  % S->G : people who are prescribed prescription opioids
eps = 0.74;         % G->S : people who use their prescriptions and then go back to susceptible
beta = 0.006;       % S->H : people who get opioids from their relatives/friends/etc to abuse them
d = 0.00824;        %      : natural death rate
dSTAR = 0.00834 ;   %      : enhanced death rate for opioid abusers
gamma =(1-eps);     % G->H : percent of prescribed opioid class who get addicted to opioids
zetaVec = 0:dz:1;   % H->R : rate at which Opioid abusers start treatment
delta = 0.75;       % R->S : people who finish their treatment and then go back to susceptible class
mu = 0.293*(1-delta);     % R->H : rate at which users in treatment fall back into drug use
sigma = 0.707*(1-delta); % R->H : rate at which people in treatment fall back into use themselves.


S_sol = zeros(length(alphaVec),length(zetaVec));
R_sol = S_sol;
H_sol = S_sol;
G_sol = S_sol;

for j = 1:length(alphaVec)
    
    alpha = alphaVec(j);
    
    for k=1:length(zetaVec)

        zeta = zetaVec(k);

        %
        % initialize
        %
        S = zeros(1,Tstop+1); G = S; R = S; H = S;
        
        %
        % initial conditions
        %
        S(1) = 0.95; 
        G(1)= 0.05; 
        R(1)= 0.00; 
        H(1) = 0.0;

        %
        % Time Evolution Step
        %
        for i=1:Tstop 
           Lambda = d*(S(i)+G(i)+R(i)) + dSTAR*H(i); 
           S(i+1) = S(i) + Lambda + delta*R(i) - alpha*S(i) - beta*S(i)*H(i) + eps*G(i) - d*S(i);
           G(i+1) = G(i) + alpha*S(i) - gamma*G(i) - eps*G(i) - d*G(i);
           R(i+1) = R(i) + zeta*H(i) - mu*R(i)*H(i) - delta*R(i) - d*R(i) - sigma*R(i);
           H(i+1) = 1-S(i+1)-G(i+1)-R(i+1);
        end

        if ( (S(end)<0) || (G(end)<0) || (R(end)<0) || (H(end)<0) )
        
            S_sol(j,k)=-1;
            G_sol(j,k)=-1;
            R_sol(j,k)=-1;
            H_sol(j,k)=-1;
            
        elseif ( (S(end)>1) || (G(end)>1) || (R(end)>1) || (H(end)>1) )
        
            S_sol(j,k)=-1;
            G_sol(j,k)=-1;
            R_sol(j,k)=-1;
            H_sol(j,k)=-1; 
            
        elseif ( ( abs( S_sol(end)-S_sol(end-1) )>1e-3 ) || ( abs( G_sol(end)-G_sol(end-1) )>1e-3 ) || ( abs( R_sol(end)-R_sol(end-1) )>1e-3 ) ||( abs( H_sol(end)-H_sol(end-1) )>1e-3 ) )

            S_sol(j,k)=1;
            G_sol(j,k)=1;
            R_sol(j,k)=1;
            H_sol(j,k)=1; 
            
        else
            
            S_sol(j,k) = S(end);
            G_sol(j,k) = G(end);
            R_sol(j,k) = R(end);
            H_sol(j,k) = H(end);
            
        end
        
        clear S R H R;
        
    end
end

save_Solutions('Sols',S_sol,G_sol,H_sol,R_sol,da,dz);

function save_Solutions(str,S_sol,G_sol,H_sol,R_sol,da,dz)

    mkdir('Solutions');
    cd('Solutions');
    S_Name = [str '.00.vtk'];
    savevtk_scalar(S_sol, S_Name, str,da,dz);

    G_Name = [str '.01.vtk'];
    savevtk_scalar(G_sol, G_Name, str,da,dz);

    R_Name = [str '.02.vtk'];
    savevtk_scalar(R_sol, R_Name, str,da,dz);

    H_Name = [str '.03.vtk'];
    savevtk_scalar(H_sol, H_Name, str,da,dz);

    cd ..

%
% Plotting solutions
%
%plot_Phase_Planes(S,G,H,R);
%plot_Time_Evolutions(t,S,G,H,R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints scalar matrix to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_scalar(array, filename, colorMap,dx,dy)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [ny, nx, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', ny, nx, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING   1.000   1.000   1.000\n'); if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
    fprintf(fid, ['SPACING   ' num2str(dx) ' '   num2str(dy) '   1.000\n']);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['SCALARS ' colorMap ' double\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:nx
            for c=1:ny
                fprintf(fid, '%d ', array(c,b,a));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Phase_Planes(x1,x2,x3,x4)

figure(1)
plot(x1,x2,'r-','LineWidth',3); hold on;
plot(x1,x3,'b-','LineWidth',3); hold on;
plot(x1,x4,'k-','LineWidth',3); hold on;
xlabel('x1');
ylabel('x2,x3,x4');
legend('G vs. S','H vs. S','R vs. S');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Evolutions(t,x1,x2,x3,x4)

figure(2)
plot(t,x1,'m-','LineWidth',4); hold on;
plot(t,x2,'r-','LineWidth',4); hold on;
plot(t,x3,'b-','LineWidth',4); hold on;
plot(t,x4,'k-','LineWidth',4); hold on;
xlabel('time');
ylabel('populations');
legend('Susceptible', 'Prescribed', 'Opioid Abuse', 'Treatment');