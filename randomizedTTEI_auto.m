function [G, TTr, time] = randomizedTTEI_auto(A, relerr, b, q)
%%
time = tic;
M = size(A);N = length(M);G = cell(1, N);TTr = zeros(1, N - 1);
relerr = relerr/sqrt(N-1);
normA = norm(A(:))^2;
relerr = relerr^2;relerr = relerr * normA;
% the first core
A = reshape(A, M(1), []);
[Q, B, k] = randQB_EI_auto(A, relerr, b, q);
TTr(1) = k;G{1} = Q;
% 2 to N-1 cores
for i = 2 : (N - 1)
    A = reshape(B, TTr(i-1) * M(i), []);
    [Q, B, k] = randQB_EI_auto(A, relerr, b, q);
    TTr(i) = k;
    G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
end
% the final core
G{N} = B;time =  toc(time);
end