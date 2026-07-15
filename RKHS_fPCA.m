function [Z_star, gamma_star, eigvals_star, f_handles, K, K_inv] = ...
    RKHS_fPCA(Y, t, kernel_fun, r, gamma_min, gamma_max, n_gamma, lambda)
% RKHS-Constrained Functional PCA
%
% Inputs:
%   Y          : (n x m) matrix of latent residual samples
%                each column = g_i^{nom}(t_j)
%   t          : (n x 1) time samples
%   kernel_fun : function handle k(t1,t2)
%   r          : number of principal components
%   gamma_min  : minimum gamma
%   gamma_max  : maximum gamma
%   n_gamma    : number of gamma grid points
%   lambda     : regularization term
%
% Outputs:
%   Z_star        : (n x r) optimal principal directions
%   gamma_star    : selected gamma
%   eigvals_star  : top r eigenvalues
%   f_handles     : cell array of functional basis handles
%   K             : Gram Matrix
%   K_inv         : inverse Gram Matrix


%% Flatten time and get size of training data, Y
t      = t(:);
[n, m] = size(Y);

%% Build Gram Matrix K 
K = zeros(n,n);
for i = 1:n
    for j = 1:n
        K(i,j) = kernel_fun(t(i), t(j));
    end
end

%% Regularize K slightly for numerical stability
K        = K + lambda*eye(n);
K_inv    = K^(-1);

%% Sample covariance 
Sigma_m  = (1/m) * (Y * Y');

%% Gamma grid
gamma_grid = linspace(gamma_min, gamma_max, n_gamma);

%% Initializations
best_score = inf;
Z_star       = [];
gamma_star   = gamma_grid(1);
eigvals_star = [];

%%  Grid search over gamma 
for g = 1:n_gamma
    gamma = gamma_grid(g);

    % Penalized matrix
    A = Sigma_m - gamma * K_inv;

    % Ensure symmetry
    A = (A + A')/2;

    % Eigen-decomposition
    [V,D] = eig(A);

    % Sort eigenvalues descending
    [eigvals, idx] = sort(diag(D), 'descend');
    V = V(:, idx);

    % Take first r eigenvectors
    Z = V(:, 1:r);

    % RKHS norm measure: trace(Z' K^{-1} Z)
    smoothness = trace(Z' * K_inv * Z);

    if smoothness < best_score
        best_score   = smoothness;
        Z_star       = Z;
        gamma_star   = gamma;
        eigvals_star = eigvals(1:r);
    end
end

%% Construct functional basis 
f_handles = cell(r,1);

for i = 1:r
    coeff = K_inv * Z_star(:,i);

    f_handles{i} = @(t_new) evaluate_kernel_basis( ...
                            t_new, t, coeff, kernel_fun);
end

end


%% Evaluate functional principal component (helper function)
function f_val = evaluate_kernel_basis(t_new, t_train, coeff, kernel_fun)
t_new      = t_new(:);
n_new      = length(t_new);
n_train    = length(t_train);
f_val      = zeros(n_new,1);

for i = 1:n_new
    for j = 1:n_train
        f_val(i) = f_val(i) + coeff(j) * kernel_fun(t_train(j), t_new(i));
    end
end

end