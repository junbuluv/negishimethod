function [p_dif,coef] = model_error_revision(lambda, coef, kgrid, param, T, psi, update, tol,idx_psi)

alpha = param.alpha;
beta = param.beta;
gamma = param.gamma;
delta = param.delta;
nss = param.nss;
dif = Inf;
iteration = 0;

% initialize 

while dif > tol

c_1 = exp(coef(1) + coef(2).* log(kgrid)) .^(-1/gamma);
c_2 = ((lambda/(1-lambda)).*c_1.^(-gamma)).^(-1/gamma);
kp = (1-delta).*kgrid - (c_1 + c_2) + kgrid.^(alpha);

c_1p = exp(coef(1) + coef(2).* log(kp)).^(-1/gamma);

e = (beta.* c_1p.^(-gamma).*(alpha.*kp.^(alpha-1) + (1-delta)));


X = [ones(T,1), log(kgrid)];

zeta = nlinfit(X,e,'object',coef);
%[L,U] = lu(X'*X);
%invXX = inv(U) * inv(L);
%zeta = invXX * X' * e;

dif = norm(coef - zeta);

if rem(iteration,100) == 0
    dif;
end

coef = update * coef + (1-update)* zeta;
iteration = iteration +1;
end


%% simulate
c1_sim = zeros(T,1);
c2_sim = zeros(T,1);
k_sim = kgrid;
n1_sim = 1/2;
w_sim = zeros(T,1);
r_sim = zeros(T,1);
denum = zeros(T,1);
p = zeros(T,1);
ltbc = zeros(T,1);

for t = 1:T
    c1_sim(t) = exp(coef(1) + coef(2)*log(k_sim(t)))^(-1/gamma);
    c2_sim(t) = ((lambda / (1-lambda)) * c1_sim(t)^(-gamma))^(-1/gamma);
    k_sim(t+1) = (1-delta)*k_sim(t) - (c1_sim(t) + c2_sim(t)) + k_sim(t)^(alpha);
    w_sim(t) = (1-alpha)*k_sim(t)^(alpha)*(nss)^(-alpha);
    r_sim(t) = alpha*k_sim(t)^(alpha-1)*(nss)^(1-alpha);
    denum(t) = 1/(1 + r_sim(t) - delta) ;
    p(1) = denum(1);
    for z = 1:T-1
    p(z+1) = p(z)*denum(z+1);
    end
    ltbc(t) = p(t)*(c1_sim(t) - w_sim(t)*(n1_sim));
end

p_dif = kgrid(1)*psi(idx_psi) - sum(ltbc) - p(T) * k_sim(T+1);

