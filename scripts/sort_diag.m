function [row_ind, col_ind] = sort_diag(x, n, dist_metric)

if nargin < 3
    dist_metric = 'correlation';
end

% sort columns based on distance
L = linkage(x', 'single', dist_metric);
col_ind = optimalleaforder(L, pdist(x',dist_metric));
x = x(:, col_ind);

% sort rows sequentially by decreasing value in each column
row_ind = [];
for i = 1:size(x, 2)
    [~,sort_ind] = sort(x(:,i), 'descend');
    row_ind = [row_ind; setdiff(sort_ind(1:n), row_ind, 'stable')];
end
row_ind = flipud(row_ind);  % since the heatmap function inverts the rows


