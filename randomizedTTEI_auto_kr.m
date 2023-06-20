function [G, TTr, time] = randomizedTTEI_auto_kr(A, relerr, b, q)
%%
time = tic;
M = size(A);N = length(M);G = cell(1, N);TTr = zeros(1, N - 1);
relerr = relerr/sqrt(N-1);
normA = norm(A(:))^2;
relerr = relerr^2;relerr = relerr * normA;
M = size(A);N = length(M);G = cell(1, N);TTr = zeros(1, N - 1);
% the first core
[Q, B] = randQB_EI_auto_kr(A, relerr, b, q);
TTr(1) = size(Q, 2);G{1} = Q;
% 2 to N-1 cores
for i = 2 : (N - 1)
    A = reshape(B, [TTr(i-1) * M(i), M(i + 1 : N)]);
    [Q, B] = randQB_EI_auto_kr(A, relerr, b, q);
    TTr(i) = size(Q, 2);
    G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
end
% the final core
G{N} = B;time =  toc(time);
end