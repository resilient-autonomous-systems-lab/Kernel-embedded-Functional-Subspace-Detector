function [Pg, Qg] = electrical_power(theta_i, G, B, n_gen, gen_buses, V, theta_fixed)
% Computes electrical active and reactive power of generators
%
% Inputs:
%   theta_i      : [n_gen x 1] rotor angles of generators 
%   G, B         : [n_bus x n_bus] network matrices 
%   n_gen        : (scalar) number of generators 
%   gen_buses    : [n_gen x 1] bus indices for generators 
%   V            : [n_bus x 1]bus voltages 
%   theta_fixed  : [n_bus x 1] reference bus angles 
%
% Outputs:
%   Pg           : [n_gen x 1] active power 
%   Qg           : [n_gen x 1] reactive power

% Preallocate outputs
Pg = zeros(n_gen,1);
Qg = zeros(n_gen,1);

% Loop over generators
for k = 1:n_gen
    i = gen_buses(k);  % bus index of the k-th generator
    
    % Active power
    temp_P = 0;
    temp_Q = 0;
    
    % Loop over all buses
    for j = 1:length(V)
        theta_diff = theta_i(k) - theta_fixed(j);
        temp_P = temp_P + V(j) * ( G(i,j) * cos(theta_diff) + B(i,j) * sin(theta_diff) );
        temp_Q = temp_Q + V(j) * ( G(i,j) * sin(theta_diff) - B(i,j) * cos(theta_diff) );
    end
    
    Pg(k) = V(i) * temp_P;
    Qg(k) = V(i) * temp_Q;
end

end