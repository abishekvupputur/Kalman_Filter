%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstration of the Linear Kalman Filter
%
%   Author: C.C. de Visser, Delft University of Technology, 2013
%   email: c.c.devisser@tudelft.nl
%   Version: 1.0
%
% Extended To Perform Nonlinear EKF
%   By: Abishek Narasimhan  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc;
clear all;

% Use these commands to initialize the randomizer to a fixed (reproducable) state.
% rng('default'); % init randomizer (default, fixed)-> version 2014a,b
% RandStream.setDefaultStream(RandStream('mt19937ar','seed', 300));-> version 2013a,b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n               = 18;
dt              = 0.01;
N               = 12000;
epsilon         = 1e-5;
maxIterations   = 100;
%Files = ["da3211.mat" ; "de3211.mat" ; "dadoublet.mat" ; "drdoublet.mat"; "dr3211.mat"];
dataset = 5;
ignore_data = 15;
do_plots = 0;


printfigs = 0;
figpath = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for states and statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m       = 6; % number of inputs

% Initial estimate for covariance matrix
stdx_0  = [1*ones(3,1); 1*ones(3,1);1*ones(9,1);1*ones(3,1)]; %TODO
P_0     = diag(stdx_0.^2);

% System noise statistics:
Ew = 0; % bias
stdw = 0.001*ones(6,1); % noise variance
Q = diag(stdw.^2);

% Measurement noise statistics:
Ev = zeros(12,1); % bias
stdv = [0.1,0.1,0.1,10,10,10,deg2rad(0.1),deg2rad(0.1),deg2rad(0.1),0.1,deg2rad(0.1),deg2rad(0.1)]; % noise variance
R = diag(stdv.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate batch with measurement data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
T2
Ex_0    = [u(1), v(1), w(1), x(1), y(1), z(1), phi(1), theta(1), psi(1),0.1*ones(1,6), 0,7,1]'; % initial estimate of optimal value of x_k_1k_1 %TODO
x_0     = [u(1), v(1), w(1), x(1), y(1), z(1), phi(1), theta(1), psi(1),0.1*ones(1,6), 0,7,1]'; % initial state %TODO

% Real simulated state-variable and measurements data:
x = x_0;
X_k = zeros(n, N);
%Z_k = zeros(nm, N);
%U_k = zeros(m, N);
Z_k = Z_k';
U_k = U_k';
XX_k1k1 = zeros(n, N);
PP_k1k1 = zeros(N, n, n);
VVe = zeros(N, n, n);
STDx_cor = zeros(n, N);
z_pred = zeros(12, N);
IEKFitcount = zeros(N, 1);

x_k_1k_1 = Ex_0; % x(0|0)=E{x_0}
P_k_1k_1 = P_0; % P(0|0)=P(0)

time1 = toc;
load('Jacobian_H.mat');
load('Jacobian_F.mat');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% Extended Kalman Filter (EKF)
ti = 0; 
tf = dt;
n = 18;%length(x_k); % n: state dimension
i=0;
% Run the filter through all N samples
for k = 1:N
    % Prediction x(k+1|k) 
    [t, x_kk_1] = rk4(@kf_calc_f, x_k_1k_1,U_k(:,k), [ti tf]); 

    % z(k+1|k) (predicted output)
    z_kk_1 = kf_calc_h(0, x_kk_1, U_k(:,k)); %x_kk_1.^3; 
    z_pred(:,k) = z_kk_1;
    
    %Get B, G matrix
    G = Get_G_matrix(x_kk_1);
    B = Get_B_matrix(x_kk_1,U_k);
    
    % Calc Phi(k+1,k) and Gamma(k+1, k)
    Fx = kf_calc_Fx(0, x, U_k(:,k), JFx); % perturbation of f(x,u,t)
    % the continuous to discrete time transformation of Df(x,u,t) and G
    [dummy, Psi] = c2d(Fx, B, dt);   
    [Phi, Gamma] = c2d(Fx, G, dt);   
    
    % P(k+1|k) (prediction covariance matrix)
    P_kk_1 = Phi*P_k_1k_1*Phi' + Gamma*Q*Gamma'; 
    P_pred = diag(P_kk_1);
    stdx_pred = sqrt(diag(P_kk_1));

    
    % Correction
    Hx = kf_calc_Hx(0, x_kk_1, U_k(:,k), JHx); % perturbation of h(x,u,t)
    % Pz(k+1|k) (covariance matrix of innovation)
    Ve = (Hx*P_kk_1 * Hx' + R); 

    % K(k+1) (gain)
    K = P_kk_1 * Hx' / Ve;
    % Calculate optimal state x(k+1|k+1) 
    x_k_1k_1 = x_kk_1 + K * (Z_k(:,k) - z_kk_1); 

    
    P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K';  
    P_cor = diag(P_k_1k_1);
    stdx_cor = sqrt(diag(P_k_1k_1));

    % Next step
    ti = tf; 
    tf = tf + dt;
    
    % store results
    XX_k1k1(:,k) = x_k_1k_1;
    PP_k1k1(k,:,:) = P_k_1k_1;
    STDx_cor(:,k) = stdx_cor;
    %VVe(k,:,:) = Ve;
    i=i+1;
    
    if mod(i,1000)==0
        i
    end
    
end

time2 = toc;
%%
%heatmap(P_k_1k_1)
for i=2:N+1
   Z_k1k1(1:12,i)=kf_calc_h(0, XX_k1k1(:,i-1), zeros(6,1));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_plots
    doPlots
end
%LSE
%%Plot Bias Estimates
%plotBiasEst
%doPlots
%convergence
File_names = ["LSE_da3211.mat" ; "LSE_de3211.mat" ; "LSE_dadoublet.mat" ; "LSE_drdoublet.mat"; "LSE_dr3211.mat"];
save(File_names(dataset), 'Z_k1k1', 'XX_k1k1', 'U_k', 'da', 'de', 'dr', 'Tc1', 'Tc2', 'dt', 'dataset')
clear all