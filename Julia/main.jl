
# initializae parameters

alpha = 1/3;
k2y = 10;
c2y = 3/4;
gamma = 2;
i2y = 1-3/4;

# steady state
nss = 1;
kss = k2y^(1/(1-alpha));
yss = kss^(alpha)*nss^(1-alpha);
delta = k2y^(-1) * (1- c2y);
css = yss - delta*kss;
beta = 1 / ( (1-delta) + alpha * kss^(alpha-1));

#psi selection
psi = [0.2,0.5,0.7,0.9];
idx_psi = 3;

# Negishi problem
T = 2500;
a = 0.01;
b = 0.99;
lambda = a+b/2;
c01_ss = lambda*css;
c02_ss = (1- lambda) * css;
# grid construction 
k_low = 0.8 * kss;
k_high = 1.2 * kss;
kgrid = LinRange(k_low,k_high,T)

# parametrize c1 using k
coef = [c01_ss; 0.0001];
# bound and update parameter
update = 0.6;
tol = 1e-5;

# Bisection search
# initialize lambda

iter = 0;
tol_search = 1e-7;
dif = Inf;
err_dif = Inf;


