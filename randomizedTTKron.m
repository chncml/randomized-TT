function [G, time] = randomizedTTKron(A, TTr, OV, q)
time = tic;
M = size(A);N = length(M);G = cell(1, N);

% the first core
L1 = TTr(1) + OV;l1 = ceil(L1^(1/(N-1)));
randG = randn(M(2), l1);
for kk = 3 : N
    temp = randn(M(kk), l1);randG = kron(randG, temp);
end
if l1^(N-1)  > L1
    randG = randG(:, 1 : L1);
end
A = reshape(A, M(1), []);B = A * randG;
for l = 1 : q
    [Q, ~] = qr(B, 0);B = A' *Q;[Q, ~] = qr(B, 0);B = A * Q;
end
[Q, ~] = qr(B, 0);Q = Q(:, 1 : TTr(1));G{1} = Q;
% 2 to N-1 cores
for i = 2 : (N-1)
    A = Q' * A;A = reshape(A, TTr(i-1) * M(i), []);
    L1 = TTr(i) + OV;l1 = ceil(L1^(1/(N-i)));randG = randn(M(i + 1), l1);
    for kk = (i + 2) : N
        temp = randn(M(kk), l1);randG = kron(randG, temp);
    end
    if l1^(N-i)  > L1
        randG = randG(:, 1 : L1);
    end
    B = A * randG;
    for l = 1 : q
        [Q, ~] = qr(B, 0);B = A' *Q;[Q, ~] = qr(B, 0);B = A * Q;
    end
    [Q, ~] = qr(B, 0);Q = Q(:, 1 : TTr(i));G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
end
% the final core
G{N} = Q' * A;
time =  toc(time);
end