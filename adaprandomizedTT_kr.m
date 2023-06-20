function [G, TTr, time] = adaprandomizedTT_kr(A, relerr, OV)
%%
time = tic;
M = size(A);N = length(M);G = cell(1, N);TTr = zeros(1, N - 1);
relerr = relerr/sqrt(N - 1);
normA = norm(A(:));
relerr = relerr * normA;
% the first core
Bmat = reshape(A, M(1), []);
[Q, mu, ~] = adap_range_finder_kr(A, relerr, OV);
TTr(1) = mu;G{1} = Q;
Bup = Q' * Bmat;
% 2 to N-1 cores
for i = 2 : (N - 1)
    Bmat = reshape(Bup, TTr(i-1) * M(i), []);
    Aten = reshape(Bmat, [TTr(i-1) * M(i), M(i+1 : N)]);
    [Q, mu, ~] = adap_range_finder_kr(Aten, relerr, OV);
    TTr(i) = mu;
    G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
    Bup = Q' * Bmat;
end
% the final core
G{N} = Bup;time =  toc(time);
end