addpath("../code/")
addpath("../utils/")

colors

n = 1000;
U = haarorth(n);
V = haarorth(n);
s = logspace(-10,0,n);
A = U*diag(s)*V';
Pinv = diag(s .* linspace(1,4,n)) \ U';
x = randn(n,1); b = A*x;
normA = norm(A);
normx = norm(x);

summary = @(y) [norm(b-A*y)/(normA*normx) norm(x-y)/normx];
[~,~,lsqrstats] = mylsqr(@(y) Pinv*(A*y), @(y) A'*(Pinv'*y), Pinv*b, 0, 100, summary);
[~,~,cgstats] = mycg(@(y) A'*(Pinv'*(Pinv*(A*y))),@(y) y,A'*(Pinv'*(Pinv*b)),0,100,summary);
summary = @(y) [norm(b-A*(A'*(Pinv'*y)))/(normA*normx) norm(x-A'*(Pinv'*y))/normx];
[~,~,cgadjstats] = mycg(@(y) Pinv*(A*(A'*(Pinv'*y))),@(y) y,Pinv*b,0,100,summary);
summary = @(y) [norm(b-A*(A'*y))/(normA*normx) norm(x-A'*y)/normx];
[~,~,pcgstats] = mycg(@(y) A*(A'*y),@(y) Pinv'*(Pinv*y),b,0,100,summary);

close all
for i = 1:2
    figure(i)
    semilogy(0:100,lsqrstats(:,i),"Color",orange,"LineWidth",3); hold on
    semilogy(0:100,cgstats(:,i),"--","Color",blue,"LineWidth",3)
    semilogy(0:100,cgadjstats(:,i),"-.","Color",pink,"LineWidth",3)
    % semilogy(0:100,pcgstats(:,i),"-.","Color",purple,"LineWidth",3)
    xlabel("Iteration $i$")
    if i == 1
        yline(norm(b-A*(A\b))/(normA*normx),":","Color",black,"LineWidth",3)
        ylabel("Residual $\|\mbox{\boldmath $b$}-\mbox{\boldmath $A$}\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $A$}\| \|\mbox{\boldmath $x$}\|$")
        exportgraphics(gcf,"../figs/left_lsqr_backward.png")
        saveas(gcf,"../figs/left_lsqr_backward.fig")
    else
        yline(norm(x-(A\b))/(normx),":","Color",black,"LineWidth",3)
        legend({"LSQR","CG on (3.1)","CG on (3.2)","Direct"},"Location","northeast")
        ylabel("Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $x$}\|$")
        exportgraphics(gcf,"../figs/left_lsqr_forward.png")
        saveas(gcf,"../figs/left_lsqr_forward.fig")
    end
end