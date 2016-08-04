clear; clc; close all; 

% Testing solving least square problems with kronecker products!
% Main idea is to solve a multivariate poly problem without ever
% constructing the full tensor matrix.

fun = @(x) exp(x(1) + x(2) );

% parameter 1
s1 = parameter('Legendre', -1, 1);
order1 = 3;
p1 = gaussian_quadrature(s1, order1);
A1 = evaluate_ops(s1,order1, p1);

% parameter 2
s2 = parameter('Legendre', -1, 1);
order2 = 4;
p2 = gaussian_quadrature(s2, order2);
A2 = evaluate_ops(s2, order2, p2);

% Test.
tensor_pts = gaussian_quadrature([s1, s2], [order1, order2]);
b = funceval(fun, tensor_pts);
pindex = index_set('tensor', [order1-1, order2-1]);
Pactual = multi_evaluate_ops(pindex, tensor_pts, [s1, s2]);
Pkron = kron(A1,A2);

% Solution
x = Pkron \ b
x2 = kron(inv(A1' * A1) * A1', inv(A2' * A2) * A2') * b 

% Using Zha's method!
[Q1,R1,P1] = qr(A1);
[Q2,R2,P2] = qr(A2);
[rows_R1, cols_R1] = size(R1);
[rows_R2, cols_R2] = size(R2);
B = reshape(b, [order2, order1]);
K = Q2' * B * Q1;
K11 = K(1:rows_R2, 1:cols_R1);
X = P2 * inv(R2) * K11 * inv(R1') * P1

% Plots
plot(Pkron', 'r-', 'LineWidth', 2); hold on;
plot(Pactual', 'b--');
