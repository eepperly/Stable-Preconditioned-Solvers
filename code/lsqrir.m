function [x,stats] = lsqrir(matvec,adjvec,b,tol,iter,varargin)
    summary = [];
    if ~isempty(varargin); summary = varargin{1}; end
    stats = [];
    x = zeros(size(adjvec(zeros(size(b,1),0)),1));
    for i = 1:length(iter)
        if i > 1 && ~isempty(stats)
            stats(end,:) = [];
        end
        [x,~,newstats] = mylsqr(matvec,adjvec,b,tol,iter(i),summary,x);
        stats = [stats;newstats];
    end
end