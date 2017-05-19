clear all; close all; clc;
%% Testing Kronlsq
A1 = rand(15,10);
A2 = rand(13,8);
b = rand(195,1);
A = kron(A1, A2);
x = A \ b
X = kronlsq(A1, A2, b)
