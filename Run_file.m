% run file for Kernel based detector for NYPA Model
% Andrew Omiloli 07/01/2026

clear
clc
close all
%% Needs from NYPA
% 1. Bus indices of:
%     a. the five solar plants [Franklin Solar, Brookside Solar, North Side Solar, Rich Road Solar Energy Center, Excelsior Energy Center]
%     b. all measurement nodes (randomizing now with some assumptions)
% 2. Average dynamic parameters for generators (Swing equation equivalent)
% 3. Typical load profiles at the demand points

%% Global parameters

%=== MATPOWER global indices ==============================================
% column indices for bus data
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

% column indices for branch data
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

% column indices for gen data
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%==========================================================================

%% Network

%=== Loading and upnpacking data ========================================
%mpc_raw = loadcase('case14');
mpc_raw = psse2mpc('NYS_5000_Zone-158-162-175_Disconnected_147786_147932_v33.raw');
mpc = ext2int(mpc_raw);

%========================================================================

% ==== Variables sizes ==================================================
% - refernce bus voltage is known (measured)

% removing redundant rows using unique branches (F_BUS, T_BUS)
branch = aggregate_branches(mpc.branch);

% removing redundant rows using unique BUS_I
bus = merge_duplicate_buses(mpc.bus);

% removing redundant rows using unique GEN_BUS
gen = aggregate_generators(mpc.gen);

% resulting aggregated mpc structure
mpc_agg        = mpc;
mpc_agg.bus    = bus;
mpc_agg.branch = branch;
mpc_agg.gen    = gen;
baseMVA        = mpc_agg.baseMVA;
gen            = mpc_agg.gen;
branch         = mpc.branch;
gen_buses      = mpc_agg.gen(:,1);
n_bus          = size(bus,1);
n_gen          = size(gen,1);

%% Dimension Variables
n_states = 2*n_gen;
p        = n_states + n_bus;
n_meas   = n_bus;
%% simulation time 
Ts = 1e-2; % sample time
%% Admittance Matrix (Ybus)
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% get conductance and subsceptance matrix
G_conduct = full(real(Ybus));   % Conductance matrix
B_suscep = full(imag(Ybus));   % Susceptance matrix

%% get voltage magnitude
Vm       = bus(:,VM);         % Voltage magnitude
theta_m  = bus(:,VA)*pi/180;  % Voltage angle 

%% Get Laplacian 
Lap = diag(Vm) * (diag(B_suscep * Vm) - B_suscep * diag(Vm));

gen_idx  = gen_buses;                       % generator buses
load_idx = setdiff(1:n_bus, gen_idx);     % remaining load buses

%% Extract loads not connected to generator
% Pd   = mpc.bus(:,3);
% Pl   = Pd(load_idx)./baseMVA;
L_gg = Lap(gen_idx, gen_idx);        % generator-generator block
L_gl = Lap(gen_idx, load_idx);       % generator-load
L_lg = L_gl.';                       % load-generator
L_ll = Lap(load_idx, load_idx);      % load-load

disp('Sanity Check Eigenvalues:')
disp(eig(full(L_gg - (L_gl*(L_ll\L_lg)))).')

%% Inertia Matrix
f_ref  = 60;
omega_r = 2*pi*f_ref;

% Generator inertia constant and damping
H_range     = [1 8];
D_range    = [1 6];
% H  = randi(H_range, n_gen, 1);      % example inertia constants
% Dg = randi(D_range, n_gen, 1);      % example damping

H = 8 * ones(n_gen,1);
Dg = 0.15 * ones(n_gen,1);

Mass_gen = 2*H/omega_r;
M = diag(Mass_gen);
Dg = diag(Dg);

%% System Matrices
%%%% state update matrix                        - 2*n_gen x 2*n_gen
A = [zeros(n_gen)                   eye(n_gen);
        -M\(L_gg-(L_gl*(L_ll\L_lg)))    -M\Dg];
%%%% input matrix                               - 2*n_gen x n_gen+n_bus
%        Pg           Pd
B = [zeros(n_gen) zeros(n_gen, n_bus- n_gen);
        inv(M)       -M\(L_gl/L_ll)];
%%%% output matrix in terms of x_bar = [delta; theta_l];
C = [zeros(n_gen) eye(n_gen);...
         -L_ll\L_lg     zeros(n_bus-n_gen, n_gen)];

D  =[zeros(n_gen)   zeros(n_gen,n_bus-n_gen);...
zeros(n_bus-n_gen,n_gen)   -L_ll^(-1)];


%% system dimensions
n_x = size(A, 1);        % number of states
n_u = size(B, 2);    % number of controls
n_y = size(C, 1);        % number of measured outputs

% checking for controllability and observability
C_mat = ctrb(A,B);
O_mat = obsv(A,C);

%% Control and observer parameters
P = -randi(10,n_x,1);  % desired poles 
K_ctrl = place(A,B,P);
L = place(A.',C.',4*P).'; % placing oberver poles further than control poles

B_perp = pinv(B);
K_r    = B_perp*A-K_ctrl;

%% Simulation parameters
t_final = 15; % final simulation time (>= 10 sec)
x_init  = randn(n_x,1); % initial system state


%% Reference generator
N      = null([A B]);
N_x    = N(1:n_x,:);
% N_u    = N(n_x+1:end,:);
lambda = rand(size(N,2),1);
x_r    = N_x*lambda;

%% Part B: Clean Residual Estimate 
%
%
%
%% (1) Extract clean training data (attack free)
data        = load('residual_nominal.mat');  
time_train  = data.residual_nominal.Time;
residual    = squeeze(data.residual_nominal.Data).';
 
%% (2) Define Kernel 
sigma      = 0.01;
kernel_fun = @(t1,t2) exp(-(t1 - t2).^2 / (2*sigma^2));
lambda     = 1e-3;  % regularization parameter

%% (3) extract training data 
n                      = round(0.8 * length(0:Ts:t_final));               % number of training samples
residual_extract       = residual(1:n,:);   % Extract measurement data
time_extract           = time_train(1:n,:); % Extract time data
[n_samples, n_signals] = size(residual_extract);

%% (4) Latent function (g_hat) initializations
g_hat_vec_all          = zeros(n_samples, n_signals);
alpha_hat_all          = cell(1, n_signals);
g_hat_fun_all          = cell(1, n_signals);

%% (5) Estimate residuls 
for i = 1:n_signals
    r_i = residual_extract(:, i);   % extract one signal
    
    [g_hat_vec, alpha_hat, g_hat_fun] = ...
    RKHS_residual_estimate(r_i, time_extract, kernel_fun, lambda);
    
    g_hat_vec_all(:, i) = g_hat_vec;
    alpha_hat_all{i}    = alpha_hat;
    g_hat_fun_all{i}    = g_hat_fun;
end

%% (6) Output Estimation Error
est_error = residual_extract - g_hat_vec_all;
overall_err = norm(est_error, 'fro');
fprintf('Frobenius norm error for clean residual RKHS estimate = %.6f\n', overall_err);

%% Part C: Learning Basis Functions
%
%
%
%% (1) Define functional PCA parameters
n_basis     = 20;   % number of basis elements      
gamma_min   = 1e-6;
gamma_max   = 10;
n_gamma     = 30;   % number of gammas
lambda_pca  = 1e-3; % regularization term


%% (2) Get optimum functional principal directions
[Z_star, gamma_star, eigvals_star, f_handles, K_train, K_inv] = ...
RKHS_fPCA(g_hat_vec_all, time_extract, kernel_fun, n_basis, gamma_min, gamma_max, n_gamma, lambda);

%% Evaluating functional bases, f_i(t) just for test purposes
% f1_vals = f_handles{1}(t_train);  % f_1(t)

%% KEFSD detection enable time (after window  is full)
Window            = 20;                 % window size
T_detect_RKHS     = (Window-1) * Ts ;   % enable time
tau_KEFSD         = 0.03;              % KEFSD threshold
tau_norm          = 0.01;               % 2 norm detetor threshold

%% detector parameters
params_detector.K_train          = K_train;
params_detector.Z_star           = Z_star;     
params_detector.tau_KEFSD        = tau_KEFSD;
params_detector.tau_norm         = tau_norm;
params_detector.lambda           = lambda;
params_detector.sigma            = sigma;
params_detector.time_extract     = time_extract;
params_detector.Window           = Window;



%% Test attack 
Attack = randn();
f       = 1; 