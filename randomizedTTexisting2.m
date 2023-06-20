function [G, time] = randomizedTTexisting2(A, TTr, OV, q)
%% T. Shi, M. Ruth and A. Townsend, Parallel algorithms for computing the tensor-train decomposition (2023)
% PSTT
time = tic;
M = size(A);N = length(M);G = cell(1, N);
% the first core
L = TTr(1) + OV;Atemp = reshape(A, M(1), []);
randG = randn(prod(M(2 : N)), L);
B = Atemp * randG;
for k=1 : q
    [Q, ~] = qr(B, 0);B = Atemp' *Q;[Q, ~] = qr(B, 0);B = Atemp * Q;
end
[Q, ~] = qr(B, 0);G{1} = Q(:, 1 : TTr(1));
% 2 to N-1 cores
Qtemp = cell(1, N - 1);
Qtemp{1} = G{1};
for i = 2 : (N - 1)
    Atemp = reshape(A, prod(M(1 : i)), []);
    L = TTr(i) + OV;randG = randn(prod(M(i + 1 : N)), L);
    B = Atemp * randG;
    for k=1 : q
        [Q, ~] = qr(B, 0);B = Atemp' *Q;[Q, ~] = qr(B, 0);B = Atemp * Q;
    end
    [Q, ~] = qr(B, 0);Q = Q(:, 1 :  TTr(i));
    Qtemp{i} = Q;
end
for i = 1 : (N - 2)
    W = Qtemp{i}' * reshape(Qtemp{i+1}, prod(M(1 : i)), []);
    G{i + 1} = reshape(W, [TTr(i), M(i + 1), TTr(i + 1)]);
end
% the final core
G{N} = Q' * Atemp;time =  toc(time);
end