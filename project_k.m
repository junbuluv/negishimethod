function [p_dif,coef] = model_error(lambda, coef, kgrid, param, T, psi, update, tol,idx_psi)

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

t = 1:1:T;

kp(t) = coef(1) + coef(2).* kgrid(t);

c_2(t) = ( 1/ (1+ ((1-lambda)/lambda)^(-1/gamma) )).*((1-delta).*kgrid(t) + kgrid(t).^(alpha)-kp(t));

c_1(t) = ((1-delta).*kgrid(t) + kgrid(t).^(alpha)-kp(t)) - c_2(t);


t_p = 1:1:T-1;
e(t_p) = (1-delta).*kgrid(t_p) + kgrid(t_p).^(alpha) - (beta.*(c_2(t_p+1).^(-gamma) .* ((1-delta) + alpha.*kgrid(t_p+1)))).^(-1/gamma) .*...
    (1 + ((1-lambda)/lambda).^(-1/gamma)) ;


X = [ones(T-1,1), kgrid(1:end-1,1)];

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
c1_sim = zeros(T,1);
c2_sim = zeros(T,1);
k_sim = zeros(T,1);
k_sim(1,1) = 0.8 * kss;
kp_sim = zeros(T,1);
n1_sim = 1/2;
w_sim = zeros(T,1);
r_sim = zeros(T,1);
denum = zeros(T,1);
p = zeros(T,1);
ltbc = zeros(T,1);

for t = 1:T
    kp_sim(t) = coef(1) + coef(2).* k_sim(t);

    c1_sim(t) = ( 1/ (1+ ((1-lambda)/lambda)^(-1/gamma) )).*((1-delta).*k_sim(t) + k_sim(t).^(alpha)-kp_sim(t));
    c2_sim(t) = ((1-delta).*k_sim(t) + k_sim(t).^(alpha)-kp_sim(t)) - c2_sim(t);
    w_sim(t) = (1-alpha)*k_sim(t)^(alpha)*(1/2)^(-alpha);
    r_sim(t) = alpha*k_sim(t)^(alpha-1)*(1/2)^(1-alpha);
    denum(t) = 1/(1-r_sim(t) + delta) ;
    p(1) = denum(1);
    for z = 1:T
    p(z+1) = p(z)*denum(t);
    end
    ltbc(t) = p(t)*(c1_sim(t) - w_sim(t)*(n1_sim));
end

p_dif = kgrid(1)*psi(idx_psi) - sum(ltbc);


