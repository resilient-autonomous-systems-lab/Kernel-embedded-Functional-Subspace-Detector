function [g_hat_vec, alpha_hat, g_hat_fun] = RKHS_residual_estimate(r, t, kernel_fun, lambda)
% RKHS residual estimation using kernel ridge regression
%
% Inputs:
%   r           : (n x 1) residual samples
%   t           : (n x 1) time samples
%   kernel_fun  : function handle k = @(t1,t2) ...
%   lambda      : regularization parameter (scalar > 0)
%
% Outputs:
%   g_hat_vec   : (n x 1) estimated latent residual at sample points
%   alpha_hat   : (n x 1) RKHS coefficients
%   g_hat_fun   : function handle for evaluating g_hat at new time points


% Ensure column vectors
r = r(:);
t = t(:);

n = length(t);

%% Build Gram matrix K
K = zeros(n,n);
for i = 1:n
    for j = 1:n
        K(i,j) = kernel_fun(t(i), t(j));
    end
end

%% Solve kernel ridge regression
% alpha_hat = (K + lambda I)^(-1) r
alpha_hat = (K + lambda*eye(n)) \ r;

%% Estimated latent residual at sample points
g_hat_vec = K * alpha_hat;

%% Function handle for evaluating at new time points
g_hat_fun = @(t_new) evaluate_kernel_expansion(t_new, t, alpha_hat, kernel_fun);

end


%% Helper function for evaluating the RKHS expansion
function g_val = evaluate_kernel_expansion(t_new, t_train, alpha, kernel_fun)

    t_new = t_new(:);
    n_new = length(t_new);
    n_train = length(t_train);

    g_val = zeros(n_new,1);

    for i = 1:n_new
        for j = 1:n_train
            g_val(i) = g_val(i) + alpha(j) * kernel_fun(t_new(i), t_train(j));
        end
    end

end