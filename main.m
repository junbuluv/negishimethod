%% Negishi method
clear all
clc
%% initializae parameters

alpha = 1/3;
k2y = 10;
c2y = 3/4;
gamma = 2;
i2y = 1-3/4;

%% steady state
nss = 1;
kss = k2y^(1/(1-alpha));
yss = kss^(alpha)*nss^(1-alpha);
delta = k2y^(-1) * (1- c2y);
css = yss - delta*kss;
beta = 1 / ( (1-delta) + alpha * kss^(alpha-1));
param = struct("nss",nss,"kss",kss,"yss",yss,"css",css,"delta",delta, "beta", beta,"gamma", gamma,"alpha",alpha);


%% psi
psi = [0.2,0.5,0.7,0.9];
idx_psi = 1;

%% Negishi problem
T = 2500;



a = 0.01;
b = 0.99;
lambda = a+b/2;
c01_ss = lambda*css;
c02_ss = (1- lambda) * css;
% grid construction 
k_low = 0.8 * kss;
k_high = 1.2 * kss;
kgrid = linspace(k_low,k_high,T)';


%% parametrize c1 using k
coef = [c01_ss; 0.0001];
%% bound and update parameter
update = 0.6;
tol = 1e-5;

%% Bisection search
% initialize lambda

iter = 0;
tol_search = 1e-5;
dif = Inf;
while dif > tol_search 

[err,coef] = model_error(lambda, coef, kgrid, param, T, psi, update, tol, idx_psi);
[err_a,~] = model_error(a,coef,kgrid,param,T,psi,update,tol,idx_psi);
    if err == tol_search
       return
    elseif sign(err) == sign(err_a)
        a = lambda;
    else
        b = lambda;
    end

    lambda = (a+b)/2;
    dif = abs(err - err_a)
    iter = iter +1;
    
end





c1_sim = zeros(T,1);
c2_sim = zeros(T,1);
k_sim = zeros(T,1);
k_sim(1,1) = 0.8 * kss;
n1_sim = 1/2;
w_sim = zeros(T,1);
r_sim = zeros(T,1);
denum = zeros(T,1);
p = zeros(T,1);
ltbc = zeros(T,1);

for t = 1:T
    c1_sim(t) = coef(1) + coef(2)*k_sim(t);
    c2_sim(t) = ((lambda / (1-lambda)) * c1_sim(t)^(-gamma) )^(-1/gamma);
    k_sim(t+1) = (1-delta)*k_sim(t) - (c1_sim(t) + c2_sim(t)) + k_sim(t)^(alpha);
    w_sim(t) = (1-alpha)*k_sim(t)^(alpha)*(1/2)^(-alpha);
    r_sim(t) = alpha*k_sim(t)^(alpha-1)*(1/2)^(1-alpha);
    denum(t) = 1/(1-r_sim(t) + delta) ;
    p(1) = denum(1);
    for z = 1:T
    p(z+1) = p(z)*denum(t);
    end
    ltbc(t) = p(t)*(c1_sim(t) - w_sim(t)*(n1_sim));
end







