function test_different_model
clc;clear;
% nn = [40, 40, 40, 40, 40];
% index = linspace(0, 1, 40);
% [I1, I2, I3, I4, I5] = ndgrid(index);
% c = rand(5, 1);
% A = sin((I1.^2 + I2.^2 + I3.^2 + I4.^2 + I5.^2).^(1/2));

% load('dataYaleB.mat');
% A = X;normA = norm(A(:));
% load('coil100A.mat');
% load('coil100B.mat');
% A(:, :, :, 1 : 50, :) = A1;A(:, :, :, 51 : 100, :) = A2;
% size(A)
% A = reshape(A, [1024, 1152, 300]);normA = norm(A(:));
% kk = [10 : 10 : 150];
load('washington.mat');
A = washington_DC_mall_double_data;
normA = norm(A(:));
kk = [10 : 20 : 150];

T1 = [];T2 = [];T7 = [];T5 = [];T6 = [];
ERR1 = [];ERR2 = [];ERR7 = [];ERR5 = [];ERR6 = [];
T11 = [];T21 = [];T71 = [];T51 = [];
ERR11 = [];ERR21 = [];ERR71 = [];ERR51 = [];
for k = 1 : length(kk)
    k
    r = kk(k);TTr = [7 * r, r];multirank = [7 * r, 2 * r, r];
    for sample = 1 : 10
        sample
        
        
        
        [Q, t1] = randomizedTTGaussian(A, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err1 = norm(A2 - tensor(A))/normA;
        
        [Q, t11] = randomizedTTGaussian(A, TTr, 10, 0);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err11 = norm(A2 - tensor(A))/normA;
          
        [Q, t4] = randomizedTTKR(A, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err4 = norm(A2 - tensor(A))/normA;
        
        [Q, t41] = randomizedTTKR(A, TTr, 10, 0);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err41 = norm(A2 - tensor(A))/normA;
        
        [Q, t5] = randomizedTTKron(A, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err5 = norm(A2 - tensor(A))/normA;
        
        [Q, t51] = randomizedTTKron(A, TTr, 10, 0);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err51 = norm(A2 - tensor(A))/normA;
        
        [Q, t6] = prandomizedTTeven(A, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
%         A2 = ttt(A2,tensor(Q{4}),[4],[1]);
%         A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err6 = norm(A2 - tensor(A))/normA;
        
        tic;[U, S] = mlsvd(A, multirank, 0);t7 = toc
        C7 = lmlragen(U, S);err7 = norm(A(:) - C7(:))/normA;
    
        t71 = tic;[U, S] = mlsvd(A, multirank, [1, 2, 3], 'LargeScale', true);t71 = toc(t71)
        C71 = lmlragen(U, S);err71 = norm(A(:) - C71(:))/normA;
        
        T1(k, sample) = t1;T4(k, sample) = t4;T7(k, sample) = t7;
        T5(k, sample) = t5;
        T11(k, sample) = t11;T41(k, sample) = t41;T71(k, sample) = t71;
        T51(k, sample) = t51;
        T6(k, sample) = t6;
        ERR1(k, sample) = err1;ERR4(k, sample) = err4;ERR7(k, sample) = err7;
        ERR5(k, sample) = err5;
        ERR11(k, sample) = err11;ERR41(k, sample) = err41;ERR71(k, sample) = err71;
        ERR51(k, sample) = err51;
        ERR6(k, sample) = err6;
        
    end
end
sample = 10;
figure(1)
subplot(121)
plot(kk, sum(ERR11, 2)/sample, '+', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR1, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR41, 2)/sample, 'o', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR4, 2)/sample, '--o', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR51, 2)/sample, 'x', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR5, 2)/sample, '-x', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR6, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR71, 2)/sample, 'd', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(ERR7, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
subplot(122)
plot(kk, sum(T11, 2)/sample, '+', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T1, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T41, 2)/sample, 'o', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T4, 2)/sample, '--o', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T51, 2)/sample, 'x', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T5, 2)/sample, '-x', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T6, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T71, 2)/sample, 'd', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T7, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
figure(2)
subplot(121)
semilogy(kk, sum(ERR11, 2)/sample, '+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR1, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR41, 2)/sample, 'o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR4, 2)/sample, '--o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR51, 2)/sample, 'x', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR5, 2)/sample, '-x', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR6, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR71, 2)/sample, 'd', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR7, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
subplot(122)
semilogy(kk, sum(T11, 2)/sample, '+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T1, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T41, 2)/sample, 'o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T4, 2)/sample, '--o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T51, 2)/sample, 'x', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T5, 2)/sample, '-x', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T6, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T71, 2)/sample, 'd', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T7, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
end