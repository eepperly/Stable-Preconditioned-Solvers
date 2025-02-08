function [y,stats] = lsqrir(Afun,Atfun,Pfun,Ptfun,b,tol,iterschedule,varargin)
    summary = [];
    if ~isempty(varargin)
        summary = varargin{1};
    end
    if length(varargin) > 1 && ~isempty(varargin{2})
        verbose = varargin{2};
    else
        verbose = false;
    end

    y = [];

    stats = [];
    matvec = @(y) Afun(Pfun(y));
    adjvec = @(y) Ptfun(Atfun(y));

    oldberr = [];
    if length(iterschedule) == 1 && tol > 0
        assert(tol > 0)
        stagnation_tol = 0.9;
        stagnation_check_frequency = iterschedule;
        iterschedule = 100*ones(50,1);
    else
        stagnation_tol = Inf;
        stagnation_check_frequency = Inf;
    end
    found_solution = false;

    if tol > 0
        n = length(b);
        x = randn(n,1);
        for i = 1:ceil(log(n))
            x = Atfun(Afun(x)); x = x / norm(x);
        end
        Anorm = norm(Atfun(x));
    end

    beta = norm(b); u = b / beta;
    for i = 1:length(iterschedule)
        v = adjvec(u); alpha = norm(v); v = v / alpha;
        if isempty(y); y = zeros(size(v)); end
        w = v;
        phibar = beta; rhobar = alpha;
        resnorm = beta * alpha;
        if ~isempty(summary) && i == 1; stats(end+1,:) = summary(y); end 
        numiters = iterschedule(i);
        for iter = 1:numiters
            u = matvec(v) - alpha*u; beta = norm(u); u = u / beta;
            v = adjvec(u) - beta*v; alpha = norm(v); v = v / alpha;
            rho = sqrt(rhobar^2 + beta^2);
            c = rhobar / rho;
            s = beta / rho;
            theta = s * alpha;
            rhobar = - c * alpha;
            phi = c * phibar;
            phibar = s * phibar;
            y = y + (phi/rho) * w;
            w = v - (theta/rho) * w;
            if ~isempty(summary); stats(end+1,:) = summary(y); end %#ok<AGROW> 
            if verbose
                fprintf('%d\t%e\n',iter,abs(phibar * alpha * c))
            end

            if mod(iter,stagnation_check_frequency) == 0 && tol > 0
                x = Pfun(y);
                berr = norm(b - Afun(x)) / (norm(x) * Anorm);

                if berr <= tol
                    found_solution = true;
                    break; 
                end

                if ~isempty(oldberr) && berr > stagnation_tol * oldberr
                    if verbose
                        fprintf("Restarting at %d/%d: %e/%e\n",iter,i,berr,oldberr)
                    end
                    oldberr = berr;
                    break;
                end

                oldberr = berr;
            end
        end

        if found_solution; break; end    
        u = b - matvec(y); beta = norm(u); u = u / beta;
    end
end