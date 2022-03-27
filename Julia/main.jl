using LinearAlgebra
using Plots
using UnPack
include("/Users/junbu/Documents/GitHub/negishimethod/Julia/model_error.jl")
# initializae parameters

alpha = 1/3;
k2y = 10;
c2y = 3/4;
gamma = 2.0;
i2y = 1-3/4;

# steady state
nss = 1.0;
kss = k2y^(1.0/(1.0-alpha));
yss = kss^(alpha)*nss^(1.0-alpha);
delta = k2y^(-1.0) * (1.0- c2y);
css = yss - delta*kss;
beta = 1.0 / ( (1.0-delta) + alpha * kss^(alpha-1.0));


#psi selection
psi = [0.2,0.5,0.7,0.9];
idx_psi = 3;
z = 1.0;
# Negishi problem
sim_length = 12500;
a = 0.01;
b = 0.99;
lambda_init = a+b/2.0;
c01_ss = lambda_init*css;
c02_ss = (1.0- lambda_init) * css;
# grid construction 
k_low = 0.8 * kss;
k_high = 1.2 * kss;
k = LinRange(k_low,k_high,sim_length);

# parametrize c1 using k
init = Vector{Float64}(undef,3);
init[1] = c01_ss;
init[2] = 0.0001;
init[3] = 0.0001;
# bound and update parameter
update_param = 0.6;
iter = 0;
tol_search = 2e-4;
ltbc_dif = Inf;

# Bisection search

while ltbc_dif > tol_search
iteration = 0;
ltbc_dif = model_error(init, lambda_init, tol_search, sim_length, k, update_param, alpha, nss, kss, psi, idx_psi);
if ltbc_dif > 0
    a = lambda_init;
elseif ltbc_dif < 0 
    b = lambda_init;
end
lambda_init = (a + b)/2
ltbc_dif = abs(ltbc_dif)
iteration += 1
if rem(iteration,10) == 0
print(ltbc_dif)
end
end

lambda_init