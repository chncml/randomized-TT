function [G, TTr, time] = adaprandomizedTT(A, relerr, OV)
%%
time = tic;
M = size(A);N = length(M);G = cell(1, N);TTr = zeros(1, N - 1);
relerr = relerr/sqrt(N - 1);
normA = norm(A(:));
relerr = relerr * normA;
% the first core
Amat = reshape(A, M(1), []);
[Q, mu, ~] = adap_range_finder(Amat, relerr, OV);
TTr(1) = mu;
G{1} = Q;
Aup = Q' * Amat;
% 2 to N-1 cores
for i = 2 : (N - 1)
    Amat = reshape(Aup, TTr(i-1) * M(i), []);
    [Q, mu, ~] = adap_range_finder(Amat, relerr, OV);
    TTr(i) = mu;
    G{i} = reshape(Q, [TTr(i-1), M(i), TTr(i)]);
    Aup = Q' * Amat;
end
% the final core
G{N} = Aup;time =  toc(time);
end