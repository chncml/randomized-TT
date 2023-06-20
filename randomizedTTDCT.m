function [G, time] = randomizedTTDCT(A, TTr, OV, q)
%%
time = tic;
M = size(A);N = length(M);G = cell(1, N);
% the first core
L = TTr(1) + OV;A = reshape(A, M(1), []);
K = prod(M(2 : N));ind = randi(2, [1, K]) * 2 - 3;
A1 = bsxfun(@times, A, ind);A1 = dct(A1')';
ind1 = sort(randsample(K, L));B = A1(:, ind1);
for k=1 : q
    [Q, ~] = qr(B, 0);B = A' *Q;[Q, ~] = qr(B, 0);B = A * Q;
end
[Q, ~] = qr(B, 0);Q = Q(:, 1 : TTr(1));G{1} = Q;
% 2 to N-1 cores
for i = 2 : (N - 1)
    A = Q' * A;A = reshape(A, TTr(i-1) * M(i), []);
    L = TTr(i) + OV;K = prod(M(i + 1 : N));
    if K <= L
        [Q, ~] = qr(A, 0);
    else
        ind = randi(2, [1, K]) * 2 - 3;A1 = bsxfun(@times, A, ind);A1 = (dct(A1'))';
        ind1 = sort(randsample(K, L));B = A1(:, ind1);
        for k=1 : q
            [Q, ~] = qr(B, 0);B = A' *Q;[Q, ~] = qr(B, 0);B = A * Q;
        end
        [Q, ~] = qr(B, 0);
    end
    Q = Q(:, 1 : TTr(i));
    G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
end
% the final core
G{N} = Q' * A;time =  toc(time);
end