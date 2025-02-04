function [x,stats] = mygmres(matvec,b,it,varargin)
summary = []; stats = [];
if ~isempty(varargin); summary = varargin{1}; end
beta = norm(b);
H = zeros(it+1,it);
V = zeros(length(b),it+1);
V(:,1) = b / beta;
if ~isempty(summary); stats(end+1,:) = summary(zeros(size(b))); end 
for j=1:it
    w = matvec(V(:,j));
    for i = 1:j
        H(i,j) = w'*V(:,i);
        w = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    V(:,j+1) = w / H(j+1,j);
    if ~isempty(summary)
        e1 = zeros(j+1,1); e1(1) = 1;
        stats(end+1,:) = summary(V(:,1:j) * (H(1:j+1,1:j) \ (beta*e1))); %#ok<AGROW>
    end 
end
e1 = zeros(it+1,1); e1(1) = 1;
x = V(:,1:it) * (H \ (beta*e1));
end