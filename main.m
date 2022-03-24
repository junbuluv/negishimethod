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

%% psi
psi = [0.2,0.5,0.7,0.9];
idx_psi = 1;

%% Negishi problem
% initialize lambda


lambda = 0.5;
T = 10000;

c01_ss = lambda*css;
c02_ss = (1- lambda) * css;
% grid construction 
k_low = 0.8 * kss;
k_high = 1.2 * kss;
kgrid = linspace(k_low,k_high,10000)';

% initialize 
c_1 = zeros(T,1);
c_2 = zeros(T,1);
kp = zeros(T,1);
c_1p = zeros(T,1);
e = zeros(T-1,1);

%% parametrize c1 using k

coef = [exp(c01_ss); 0.0001];



%% bound and update parameter

update = 0.6;
tol = 1e-7;
iter_limit = T;
dif = Inf;
iteration = 0;
%% loop

while dif > tol

t = 1:1:T;
c_1(t) = coef(1) + coef(2).*log(kgrid(t));
c_2(t) = ((lambda / (1-lambda)) .* c_1(t).^(-gamma) ).^(-1/gamma);
kp(t) = (1-delta).*kgrid(t) - (c_1(t) + c_2(t)) + kgrid(t).^(alpha);

c_1p = c_1(2:end);

t_p = 1:1:T-1;
e(t_p) = (beta.* c_1p(t_p).^(-gamma).*(alpha*kp(t_p).^(alpha-1) + (1-delta))).^(-1/gamma);


X = [ones(T-1,1), log(kgrid(1:end-1,1))];

%zeta = nlinfit(X,e,'fit',coef);
[L,U] = lu(X'*X);
invXX = inv(U) * inv(L);
zeta = invXX * X' * e;

dif = norm(coef - zeta);

if rem(iteration,100) == 0
    dif
end

coef = update * coef + (1-update)* zeta;
iteration = iteration +1;
end


%% simulate
c1_sim = zeros(T,1);
c2_sim = zeros(T,1);
k_sim = zeros(T,1);
k_sim(1,1) = 0.8 * kss;
n1_sim = 1/2;
n2_sim = 1/2;
w_sim = zeros(T,1);
r_sim = zeros(T,1);
denum = zeros(T,1);
p = zeros(T,1);

for t = 1:T
    c1_sim(t) = exp(coef(1) + coef(2)*log(k_sim(t)));
    c2_sim(t) = ((lambda / (1-lambda)) * c1_sim(t)^(-gamma) )^(-1/gamma);
    k_sim(t+1) = (1-delta)*k_sim(t) - (c1_sim(t) + c2_sim(t)) + k_sim(t)^(alpha);
    w_sim(t) = (1-alpha)*k_sim(t)^(alpha)*(1/2)^(-alpha);
    r_sim(t) = alpha*k_sim(t)^(alpha-1)*(1/2)^(1-alpha);
    denum(t) = (1-r_sim(t) + delta)^(-1) ;
    p(1) = denum(1);
    p(t+1) = p(t)*denum(t);
    ltbc(t) = p(t)*(c1_sim(t) - w_sim(t)*(n1_sim));
end

ltbc = p(t)*(c1_sim(t) - w_sim(t)*(n1_sim));
p_dif = kgrid(1)*psi(idx_psi) - sum(ltbc);

