function Kron_test
clear all; close all; clc;

% Testing solving least square problems with kronecker products!
% Main idea is to solve a multivariate poly problem without ever
% constructing the full tensor matrix.

fun = @(x) x(1) + x(2) ;

% parameter 1
s1 = parameter('Legendre', -1, 1);
order1 = 3;
[p1, w1] = gaussian_quadrature(s1, order1);
A1 = evaluate_ops(s1,order1, p1)';

% parameter 2
s2 = parameter('Legendre', -1, 1);
order2 = 4;
[p2, w2] = gaussian_quadrature(s2, order2);
A2 = evaluate_ops(s2, order2, p2)';

% Test.
tensor_pts = gaussian_quadrature([s1, s2], [order1, order2]);
b = funceval(fun, tensor_pts);

% Using conjugate gradient solver!
A{1} = A1;
A{2} = A2;
maxit = 60;
tol = 1e-8
x1 = cgs(@afun, b, tol, maxit)

    function y = afun(x)
        y = kronmult(A, x);
    end
end

