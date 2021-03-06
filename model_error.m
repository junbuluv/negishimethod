function [p_dif,coef] = model_error(lambda, coef, kgrid, param, T, psi, update, tol,idx_psi)

alpha = param.alpha;
beta = param.beta;
gamma = param.gamma;
delta = param.delta;
nss = param.nss;
kss = param.kss;
dif = Inf;
iteration = 0;

% initialize 
c_1 = zeros(T,1);
c_2 = zeros(T,1);
kp = zeros(T,1);
e = zeros(T,1);

while dif > tol

t = 1:1:T;
c_1(t) = (coef(1) + coef(2).*kgrid(t));
c_2(t) = ((lambda / (1-lambda)) .* c_1(t).^(-gamma) ).^(-1/gamma);
kp(t) = (1-delta).*kgrid(t) - (c_1(t) + c_2(t)) + kgrid(t).^(alpha);

c_1p = (coef(1) + coef(2).*kp(t));


e(t) = (beta.* c_1p(t).^(-gamma).*(alpha*kp(t).^(alpha-1) + (1-delta))).^(-1/gamma);


X = [ones(T,1), kgrid(1:end,1)];

%zeta = nlinfit(X,e,'fit',coef);
[L,U] = lu(X'*X);
invXX = inv(U) * inv(L);
zeta = invXX * X' * e;

dif = norm(coef - zeta);

if rem(iteration,100) == 0
    dif;
end

coef = update * coef + (1-update)* zeta;
iteration = iteration +1;
end


%% simulate
k_sim = zeros(T,1);
k_sim(1,1) = 0.8*kss;
n1_sim = 1/2;
p = zeros(T,1);
c1_sim = zeros(T,1);
c2_sim = zeros(T,1);
w_sim = zeros(T,1);
r_sim = zeros(T,1);
denum = zeros(T,1);
ltbc = zeros(T,1);

for t = 1: T
    c1_sim(t) = coef(1) + coef(2)*k_sim(t);
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
p_dif = kgrid(1)*psi(idx_psi) - sum(ltbc);


