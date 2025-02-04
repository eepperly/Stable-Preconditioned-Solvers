rng(1323901)
n = 1e3;

%% Matrix 1
U = haarorth(n);
s = linspace(1,2,n);
V = haarorth(n);
W = U*diag(s)*V';
A1 = W*diag(logspace(0,2,n))/W;

%% Matrix 2
A2 = U * diag(1 ./ linspace(1,4,n));
% A2 = U * diag((1+rand(1,n)) .* exp(1i*pi/4*rand(1,n))) * U';

%% Matrix 3
A3 = W * diag(1 ./ linspace(1,4,n)) / W;

b = randn(n,1);
nsteps = [200 200 200];

As = {A1,A2,A3};

close all
figure('Position', [100, 100, 1800, 400])  

for A_idx = 1:3
    A = As{A_idx};
    steps = nsteps(A_idx);
    mycond = cond(A);
    cgrate = (sqrt(mycond) - 1)/(sqrt(mycond)+1);
    lsqrrate = (mycond - 1)/(mycond+1);
    normA = norm(A);
    normx = norm(x);
    
    % Create subplot
    subplot(1,3,A_idx)
    
    [~,~,stats] = mylsqr(@(y) A*y,@(y) A'*y,b,0,steps,@(y) norm(b-A*y)/normA/normx);
    semilogy(0:(length(stats)-1),stats,"Color",orange,"LineWidth",3); hold on
    semilogy(0:steps,lsqrrate .^ (0:steps),"--","Color",black,"LineWidth",3)
    [~,stats] = mygmres(@(y) A*y,b,steps,@(y) norm(b-A*y)/normA/normx);
    semilogy(0:steps,stats,"-.","Color",pink,"LineWidth",3)
    semilogy(0:steps,cgrate .^ (0:steps),":","Color",black,"LineWidth",3)
    axis([0 steps 1e-17 1e0])
    
    xlabel("Iteration $i$")
    if A_idx == 1
        legend({"LSQR","$\varrho_{\rm LSQR}^i$","GMRES","$\varrho_{\rm CG}^i$"},"Location","best")
        ylabel("Residual $\|\mbox{\boldmath $b$}-\mbox{\boldmath $A$}\mbox{\boldmath $x$}_i\| / \|\mbox{\boldmath $A$}\| \|\mbox{\boldmath $x$}\|$")
    elseif A_idx == 2
        pos = get(gca, 'Position');  % Get position of the current subplot
        inset_pos = [pos(1) + 0.4*pos(3), pos(2) + 0.3*pos(4), 0.5*pos(3), 0.5*pos(4)];  % Position inset inside the subplot
        
        ax_inset = axes('Position', inset_pos);
        
        evals = eig(A);
        
        % Plot scatter plot in the inset
        scatter(ax_inset, real(evals), imag(evals), 10, 'filled', "MarkerFaceColor", purple, 'MarkerFaceAlpha', 0.4)
        title(ax_inset, 'Eigenvalues')
        xlabel("Re")
        ylabel("Im")
        % axis(ax_inset, 'off')
        set(ax_inset, 'Box', 'off')
        ax_inset.XAxisLocation = 'origin';
        ax_inset.YAxisLocation = 'origin';
        % axis(ax_inset, 'equal')
        axis([-0.7 0.7 -0.7 0.7])
        % set(ax_inset, 'Box', 'on')

        xticks(ax_inset, [])
        yticks(ax_inset, [])
        % xticklabels(ax_inset, {})
        % yticklabels(ax_inset, {'-1', '1'})
        
        % Add arrows for the x-axis and y-axis
        annotation('arrow', [inset_pos(1), inset_pos(1) + inset_pos(3)], [inset_pos(2)+0.5*inset_pos(4), inset_pos(2)+0.5*inset_pos(4)], 'Color', 'black', 'LineWidth', 1.5);  % X-axis arrow
        annotation('arrow', [inset_pos(1)+0.5*inset_pos(3), inset_pos(1)+0.5*inset_pos(3)], [inset_pos(2), inset_pos(2) + inset_pos(4)], 'Color', 'black', 'LineWidth', 1.5);  % Y-axis arrow
    end
end

exportgraphics(gcf,'../figs/gmres_lsqr.png')
saveas(gcf,'../figs/gmres_lsqr.fig')

drawnow

