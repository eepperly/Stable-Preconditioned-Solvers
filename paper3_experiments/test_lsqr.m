addpath('../code')
addpath('../utils')
colors

rng(134200343)

n = 1000;
U = haarorth(n);
V = haarorth(n);
s = logspace(-10,0,n);
A = U*diag(s)*V';
Pinv = V / diag(s .* linspace(1,4,n));
x = randn(n,1); b = A*x;
cond(A*Pinv)
normA = norm(A);
normx = norm(x);
summary = @(y) [norm(b-A*(Pinv*y))/(normA*normx) norm(x-Pinv*y)/normx norm(b-A*(Pinv*y))/(normA*norm(Pinv*y))];
niter = 100;

[~,~,lsqr_res] = mylsqr(@(y) A*(Pinv*y), @(y) Pinv'*(A'*y),b,0,niter,summary,[],[]);
[~,lsqrir_res] = lsqrir(@(y) A*y, @(y) A'*y, @(y) Pinv*y, @(y) Pinv'*y, b, 0, [50,50], summary);
[~,gmres_res] = mygmres(@(y) A*(Pinv*y),b,niter,summary);

close all
figure
semilogy(0:niter,lsqr_res(:,1),"Color",blue,"LineWidth",3); hold on
semilogy(0:(length(lsqrir_res(:,1))-1),lsqrir_res(:,1),"--","Color",orange,"LineWidth",3)
semilogy(0:niter,gmres_res(:,1),"-.","Color",pink,"LineWidth",3)
yline(norm(b-A*(A\b))/(normA*normx),":","Color",black,"LineWidth",3)
axis([0 niter 1e-17 1e0])
legend({"PLSQR","PLSQR-IR","PGMRES","Direct"},"Location","best")
xlabel("Iteration $i$")
ylabel("Residual $\|\mbox{\boldmath $b$}-\mbox{\boldmath $A$}\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $A$}\| \|\mbox{\boldmath $x$}\|$")
exportgraphics(gcf,"../figs/lsqr_backward.png")
saveas(gcf,"../figs/lsqr_backward.fig")

figure
semilogy(0:niter,lsqr_res(:,2),"Color",blue,"LineWidth",3); hold on
semilogy(0:(length(lsqrir_res(:,2))-1),lsqrir_res(:,2),"--","Color",orange,"LineWidth",3)
semilogy(0:niter,gmres_res(:,2),"-.","Color",pink,"LineWidth",3)
yline(norm(x-(A\b))/(normx),":","Color",black,"LineWidth",3)
axis([0 niter -Inf Inf])
xlabel("Iteration $i$")
ylabel("Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $x$}\|$")
exportgraphics(gcf,"../figs/lsqr_forward.png")
saveas(gcf,"../figs/lsqr_forward.fig")

figure
semilogy(0:niter,lsqr_res(:,3),"Color",blue,"LineWidth",3); hold on
semilogy(0:(length(lsqrir_res(:,3))-1),lsqrir_res(:,3),"--","Color",orange,"LineWidth",3)
semilogy(0:niter,gmres_res(:,3),"-.","Color",pink,"LineWidth",3)
yline(norm(b-A*(A\b))/(normA*norm(A\b)),":","Color",black,"LineWidth",3)
xlabel("Iteration $i$")
ylabel("Backward error $\|\mbox{\boldmath $b$}-\mbox{\boldmath $A$}\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $A$}\| \|\mbox{\boldmath $x$}_i\|$")
exportgraphics(gcf,"../figs/lsqr_backward_2.png")
saveas(gcf,"../figs/lsqr_backward_2.fig")