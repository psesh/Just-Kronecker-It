function X = kronlsq(A1, A2, b)
%--------------------------------------------------------------------------
% Kronecker product least squares without every computing the full
% Kronecker product! 
%
% A1: m-by-p matrix, where m>p
% A2: n-by-q matrix, where n>q
% b: mn-by-1 vector
%
% Pranay Seshadri
% University of Cambridge
% May 19th, 2017
%--------------------------------------------------------------------------
[m, p] = size(A1);
[n, q] = size(A2);
B = reshape(b, [n, m]);
[Q1, R1] = qr(A1, 0);
[Q2, R2] = qr(A2, 0);
K = Q2' * B * Q1;
X = inv(R2) * K * inv(R1') ; 
end
