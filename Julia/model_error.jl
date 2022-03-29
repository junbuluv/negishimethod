

function model_error(coef, lambda, tol, T, kgrid, update, alpha, nss, kss, psi, idx_psi)
    dif = Inf;
    iter = 0;
    c_1 = zeros(T,1);
    c_2 = zeros(T,1);
    kp = zeros(T,1);
    e = zeros(T,1);
   while dif > tol
    
        c_1 = (coef[1] .+ coef[2]*kgrid);
        c_2 = ((lambda/(1.0-lambda)) .* c_1.^(gamma)).^(-1/gamma) ;
        kp = (1.0-delta).*kgrid - (c_1 + c_2) + kgrid.^(alpha);
        c_1p =(coef[1] .+ coef[2]*kp);
     
        # fitting
        e = (beta.* c_1p.^(-gamma).*(alpha*kp.^(alpha-1) .+ (1.0-delta))).^(-1.0/gamma);
        X = [ones(T,1) kgrid];
        zeta = X\e;
        dif = norm(zeta - coef);
        if mod(iter,20) == 0
            dif
        end
        coef = (1.0-update) * coef + update * zeta;
        iter = iter +1.0;
        
    end
    
    k_sim = zeros(T+1);
    k_sim[1] = 0.8 * kss;
    c1_sim = zeros(T);
    c2_sim = zeros(T);
    n1_sim = 1/2;
    w = zeros(T);
    r = zeros(T);
    idx_n = zeros(T);
    idx_z = zeros(T-1);
    p = Vector{Float64}(undef,T)
    denum = zeros(T);
    p = zeros(T);
    ltbc = zeros(T);
    for t in eachindex(idx_n)
        ## other variables
        c1_sim[t] = coef[1] + coef[2] * k_sim[t];
        c2_sim[t] = ((lambda/ (1-lambda)) * c1_sim[t]^(-gamma))^(-1/gamma);
        k_sim[t+1] = (1-delta)*k_sim[t] - (c1_sim[t] + c2_sim[t]) + k_sim[t]^(alpha) *nss^(1-alpha);
        w[t] = (1.0-alpha)*k_sim[t]^alpha * nss^(-alpha);
        r[t] = alpha*k_sim[t]^(alpha-1.0) * nss^(1.0-alpha);
        denum[t] = 1.0 / (1.0 - delta + r[t])
        p[1] = denum[1];
        for z in eachindex(idx_z)
            p[z+1] = p[z]*denum[z+1];
        end
    
        ltbc[t] = p[t]*(c1_sim[t] - w[t] * (n1_sim));
    end
    
    val = kss*psi[idx_psi] - sum(ltbc)
    
    return val
       
    end
