function P = multi_evaluate_ops(pindex, points, s)

% Generate a multi-variate orthogonal polynomial at points
[rows, cols] = size(pindex);

% Number of polynomails
for k = 1 : rows
    p{k} = evaluate_ops(s(k), max(pindex(k,:) + 1), points(:,k) ); % get maximum order for each direction!
end


for i = 1 : cols
   [no_quad_pts, dims ] = size(points);
   temp = ones(1,no_quad_pts);
    for k = 1 : rows
        P(i,:) = p{k}(pindex(k,i)+1, :) .* temp;
        temp = P(i,:);
   end
    
end
end
