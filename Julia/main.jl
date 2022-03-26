using LinearAlgebra
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
idx_psi = 3.0;
z = 1.0;
# Negishi problem
T = 2500;
a = 0.01;
b = 0.99;
lambda = a+b/2.0;
c01_ss = lambda*css;
c02_ss = (1.0- lambda) * css;
# grid construction 
k_low = 0.8 * kss;
k_high = 1.2 * kss;
kgrid = LinRange(k_low,k_high,T)

# parametrize c1 using k
coef = Vector{Float64}(undef,3);
coef[1] = c01_ss;
coef[2] = 0.0001;
coef[3] = 0.0001;
# bound and update parameter
update = 0.6;
tol = 1e-5;

# Bisection search
# initialize lambda

iter = 0;
tol_search = 5e-5;
dif = Inf;
err_dif = Inf;


c_1 = zeros(T,1);
c_2 = zeros(T,1);
kp = zeros(T,1);
e = zeros(T,1);

while dif > tol

    c_1 = (coef[1] .+ coef[2]*kgrid .+ coef[3]*log.(kgrid));
    c_2 = ((lambda/(1.0-lambda)) .* c_1.^(gamma)).^(-1.0/gamma) ;
    kp = (1-delta).*kgrid - (c_1 + c_2) + kgrid.^(alpha);
    c_1p =(coef[1] .+ coef[2]*kgrid .+ coef[3]*log.(kgrid));

    # fitting
    e = (beta.* c_1p.^(-gamma).*(alpha*kp.^(alpha-1) .+ (1.0-delta))).^(-1.0/gamma);
    X = [ones(T,1) kgrid log.(kgrid)];
    L,U = lu(transpose(X) * X);
    invxx = inv(U) * inv(L);
    zeta = invxx * (transpose(X) * e) 

    dif = norm(zeta - coef);
    if mod(iter,20) == 0
        dif
    end
    coef = (1-update) * coef + update * zeta;
    iter = iter +1;
    
end
    ## other variables
    w = (1.0-alpha)*kgrid.^alpha * nss.^(-alpha);
    r = alpha*kgrid.^(1.0-alpha) * nss.^(1.0-alpha);
    denum = 1.0 ./ (1.0 .- delta .+ r)
    p = Vector{Float64}(undef,T)
    p[1] = denum[1];
    for z = 1:T-1
        p[z+1] = p[z]*denum[z+1]
    end
