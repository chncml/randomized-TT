function [G, time] = prandomizedTTeven(A, TTr, OV, q)
%% input
% target tensor a; the desired TT-rank TTr; oversampling paramter OV
%% output
% all the core tensors stored in a cell G; running time time
time = tic;
M = size(A);N = length(M);G = cell(1, N);
L = TTr(1) + OV;A = reshape(A, M(1), []);
Q = randn(M(1), L);
for k = 1 : q
    B = A' * Q;[Q, ~] = qr(B, 0);B = A * Q;[Q, ~] = qr(B, 0);
end
Q = Q(:, 1 : TTr(1));G{1} = Q;
for i = 2 : N - 1
    A = Q' * A;A = reshape(A, TTr(i-1) * M(i), []);
    L = TTr(i) + OV;
    Q = randn(TTr(i-1) * M(i), L);
    for k = 1 : q
        B = A' * Q;[Q, ~] = qr(B, 0);B = A * Q;[Q, ~] = qr(B, 0);
    end
    Q = Q(:, 1 : TTr(i));G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
end
G{N} = Q' * A;time =  toc(time);
end