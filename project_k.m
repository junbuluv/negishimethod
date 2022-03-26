function [p_dif,coef] = project_k(lambda, coef, kgrid, param, T, psi, update, tol,idx_psi)

alpha = param.alpha;
beta = param.beta;
gamma = param.gamma;
delta = param.delta;
kss = param.kss;
dif = Inf;
iteration = 0;

% initialize 
c_1 = zeros(T,1);
c_2 = zeros(T,1);
kp = zeros(T,1);
e = zeros(T-1,1);

while dif > tol


kp = exp(coef(1) + coef(2)*log(kgrid));

c_2 = (1/(1+((1-lambda)/lambda)^(-1/gamma)))*((1-delta)*kgrid + kgrid.^(alpha)-kp);

c_1 = ((1-delta)*kgrid + kgrid.^(alpha)-kp) - c_2;


e = ((1/(alpha*beta)).*(c_1(1:end-1)./cp_1).^(-gamma) - (1-delta)/alpha).^(1/(alpha-1)) ;


X = [ones(T-1,1), log(kgrid(1:end-1,1))];

zeta = nlinfit(X,e,'fit',coef);
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
kp_sim = zeros(T,1);
n1_sim = 1/2;
w_sim = zeros(T,1);
r_sim = zeros(T,1);
denum = zeros(T,1);
p = zeros(T,1);
ltbc = zeros(T,1);

for t = 1:T
    kp_sim(t) = exp(coef(1) + coef(2)*log(k_sim(t)));
    c2_sim(t) = (1/(1+((1-lambda)/lambda)^(-1/gamma)))*((1-delta)*k_sim(t) + k_sim(t).^(alpha)-kp_sim(t));
    c1_sim(t) = ((1-delta)*k_sim(t) + k_sim(t).^(alpha)-kp_sim(t)) - c2_sim(t);
    w_sim(t) = (1-alpha)*k_sim(t)^(alpha)*(1)^(-alpha);
    r_sim(t) = alpha*k_sim(t)^(alpha-1)*(1)^(1-alpha);
    denum(t) = 1/(1 + r_sim(t) - delta) ;
    p(1) = denum(1);
    for z = 1:T
    p(z+1) = p(z)*denum(t);
    end
    ltbc(t) = p(t)*(c1_sim(t) - w_sim(t)*(n1_sim));
end

p_dif = kgrid(1)*psi(idx_psi) - sum(ltbc);


