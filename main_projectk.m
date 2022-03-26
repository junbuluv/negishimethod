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
coef = [log(c01_ss^(-gamma)); 0.0001 ];
%% bound and update parameter
update = 0.6;
tol = 1e-5;

%% Bisection search
% initialize lambda

iter = 0;
tol_search = 1e-7;
dif = Inf;
err_dif = Inf;


while dif > tol_search

[err,coef] = project_k(lambda, coef, kgrid, param, T, psi, update, tol, idx_psi);
[err_b,~] = project_k(b,coef,kgrid,param,T,psi,update,tol,idx_psi);
    if abs(err) < tol_search
       return
    elseif sign(err) == sign(err_b)
        b = lambda;
    else
        a = lambda;
    end
    lambda = (a+b)/2;
    dif = abs(err - err_b)
    iter = iter +1;
    
end



sim_length = 200;
c1_sim = zeros(sim_length,1);
c2_sim = zeros(sim_length,1);
k_sim = kgrid;
k_sim(1,1) = 0.8* kss;
kp_sim  = zeros(sim_length,1);
for t = 1:sim_length
    kp_sim(t) = exp(coef(1) + coef(2)*log(k_sim(t)));
    c2_sim(t) = (1/(1+((1-lambda)/lambda)^(-1/gamma)))*((1-delta)*k_sim(t) + k_sim(t).^(alpha)-kp_sim(t));
    c1_sim(t) = ((1-delta)*k_sim(t) + k_sim(t).^(alpha)-kp_sim(t)) - c2_sim(t);
end


hold on
plot(c1_sim)
plot(c2_sim)
hold off


