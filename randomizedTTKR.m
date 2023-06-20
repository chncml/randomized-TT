function [G, time] = randomizedTTKR(A, TTr, OV, q)
time = tic;
M = size(A);N = length(M);G = cell(1, N);
% the first core
L = TTr(1) + OV;A = reshape(A, M(1), []);
randG = randn(M(2), L);
for k = 3 : N
    temp = randn(M(k), L);randG = kr(randG, temp);
end
B = A * randG;
for k=1 : q
    [Q, ~] = qr(B, 0);B = A' *Q;[Q, ~] = qr(B, 0);B = A * Q;
end
[Q, ~] = qr(B, 0);Q = Q(:, 1 : TTr(1));G{1} = Q;
% 2 to N-1 cores
for i = 2 : (N - 1)
    A = Q' * A;A = reshape(A, TTr(i-1) * M(i), []);
    L = TTr(i) + OV;randG = randn(M(i + 1), L);
    for k = i + 2 : N
        temp = randn(M(k), L);randG = kr(randG, temp);
    end
    B = A * randG;
    for k=1 : q
        [Q, ~] = qr(B, 0);B = A' *Q;[Q, ~] = qr(B, 0);B = A * Q;
    end
    [Q, ~] = qr(B, 0);
    Q = Q(:, 1 : TTr(i));G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
end
% the final core
G{N} = Q' * A;time =  toc(time);
end