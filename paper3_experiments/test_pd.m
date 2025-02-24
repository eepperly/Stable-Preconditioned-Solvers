addpath("../code")
addpath("../utils")

colors

rng(12349301)

n = 1e3;
cond_A = 1e14;
steps = 50;
refines = 3;
U = haarorth(n);
A = U * diag(logspace(0,-log10(cond_A),n)) * U';
G = randn(4*n,n);
P = U * diag(logspace(0,-log10(cond_A)/2,n)) * (G' * G) ...
    * diag(logspace(0,-log10(cond_A)/2,n)) * U';
R = chol(P);
x = randn(n,1);
b = A*x;
normA = norm(A);
normx = norm(x);
Avpa = vpa(A,24);
bvpa = vpa(b,24);

summary = @(xx) norm(bvpa-Avpa*vpa(xx,24)) / (normA*normx);
[~,~,pcgstats] = mycg(@(x) A*x,@(x) R\(R'\x),b,0,(refines+1)*steps,summary);

[xx,~,pcgirstats] = mycg(@(x) A*x,@(x) R\(R'\x),b,0,steps,summary);
for refine = 1:refines
    summary = @(y) norm(bvpa-Avpa*vpa(y+xx,24)) / (normA*normx);
    [dx,~,newstats] = mycg(@(x) A*x,@(x) R\(R'\x),b-A*xx,0,steps,summary);
    xx = xx + dx;
    pcgirstats(end+1:end+(length(newstats)-1)) = newstats(2:end);
end

pmatrix = R' \ (A / R); pmatrix = (pmatrix+pmatrix')/2;
evals2 = eig(pmatrix);

load("../data/COMET_MC_SAMPLE.mat")

n = 10000; k = 850;
S = datasample(1:size(Xtr,1),n,"Replace",false);
X = Xtr(S,:);
X = (X - mean(X)) ./ std(X);
b = Ytr(S).';
A = exp(-pdist2(X,X).^2/(2*size(X,2)));
F = rpcholesky(@(i) A(:,i),ones(size(A,1),1),k);
[U,s,~] = svd(F,"econ","vector");
lamb = 1e-10;
d = 1 ./ (s.^2 + lamb);
A = A + lamb*eye(n);
% x = A\b;
x = randn(n,1);
b = A*x;
normA = norm(A);
normx = norm(x);
refines = 1;
steps = 150;
pfun = @(y) apply_pre(y,U,d,lamb);

summary = @(y) norm(b-A*y) / (normA*normx);
[x0,~,pcgstats2] = mycg(@(y) A*y,pfun,b,0,(refines+1)*steps,summary);
[xx,~,pcgirstats2] = mycg(@(x) A*x,pfun,b,0,steps,summary);
for refine = 1:refines
    summary = @(y) norm(b-A*(y+xx)) / (normA*normx);
    [dx,~,newstats] = mycg(@(x) A*x,pfun,b-A*xx,0,steps,summary);
    xx = xx + dx;
    pcgirstats2(end+1:end+(length(newstats)-1)) = newstats(2:end);
end

% cond(A)
% evals2 = eig(apply_pre(A,U,d,lamb));

%% Make figures

close all
figure('Position', [100, 100, 1100, 400])  
subplot(1,2,1)
semilogy(0:(length(pcgstats)-1),pcgstats,"Color",blue,"LineWidth",3); hold on
semilogy(0:(length(pcgirstats)-1),pcgirstats,"--","Color",orange,"LineWidth",3)
yline(norm(b-A*(A\b)) / (normA*normx),":","Color",black,"LineWidth",3)
axis([-Inf Inf 1e-17 1e0])
xlabel("Iteration $i$")
ylabel("Residual $\|\mbox{\boldmath $b$}-\mbox{\boldmath $A$}\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $A$}\| \|\mbox{\boldmath $x$}\|$")
legend({"PCG","PCG-IR","Direct"})

subplot(1,2,2)
semilogy(0:(length(pcgstats2)-1),pcgstats2,"Color",blue,"LineWidth",3); hold on
semilogy(0:(length(pcgirstats2)-1),pcgirstats2,"--","Color",orange,"LineWidth",3)
yline(norm(b-A*(A\b)) / (normA*normx),":","Color",black,"LineWidth",3)
axis([-Inf Inf 1e-17 1e0])
xlabel("Iteration $i$")

exportgraphics(gcf,'../figs/pd.png')
saveas(gcf,'../figs/pd.fig')

function z = apply_pre(y,U,d,lamb)
    Uy = U'*y;
    z = U * (d .* Uy) + (y - U*Uy) / lamb;
end