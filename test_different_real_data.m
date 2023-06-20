function test_different_real_data
%% a test tensor from yaleB
clc;clear;
load('dataYaleB.mat')
A = X(1 : 451, 1 : 602, 1 : 540);
A = reshape(A, [41, 42, 43, 44, 45]);
nn = [41, 42, 43, 44, 45];d = 5;
% load('washington.mat');
% A = washington_DC_mall_double_data;
% nn = [35 37 38 39 39];A = A(:);A = A(1 : prod(nn));
% A = reshape(A, nn);d = 5;
% load('salesman_qcif.mat')
% A = double(X/255);A = A(:, :, :, 1 : 440);
% A = reshape(A, [30, 32, 32, 33, 33]);
% nn = [30, 32, 32, 33, 33];d = 5;
% load('coil100A.mat');
% load('coil100B.mat');
% A(:, :, :, 1 : 50, :) = A1;A(:, :, :, 51 : 100, :) = A2;
% A = reshape(A, [64, 64, 50, 72, 24]);nn = [64, 64, 50, 72, 24];d = 5;
% load('salesman_qcif.mat')
% A = double(X/255);A = A(:, :, :, 1 : 440);
% A = reshape(A, [30, 32, 32, 33, 33]);
% nn = [30, 32, 32, 33, 33];d = 5;
B = A;A = tensor(A);normA = norm(A);

AT1 = [];AT2 = [];AT3 = [];AT4 = [];AT5 = [];
AERR1 = [];AERR2 = [];AERR3 = [];AERR4 = [];AERR5 = [];R = [];
kk = [1 : 0.5 : 5];
for P = 1 : length(kk)
    P
    k = kk(P);
    for sample = 1 : 10
        sample
        
        t5 = tic;Q = tt_tensor(B, 10^(-k));t5 = toc(t5);
        R = rank(Q);
        err5 = norm(full(Q) - B(:))/normA;
        TTr = R(2 : end-1);
        
        [Q, t1] = randomizedTTGaussian(B, TTr, 10, 0);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err1 = norm(A2 - A)/normA;
        
        [Q, t11] = randomizedTTGaussian(B, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err11 = norm(A2 - A)/normA;
        
        [Q, t2] = randomizedTTKR(B, TTr, 10, 0);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err2 = norm(A2 - A)/normA;
        
        [Q, t21] = randomizedTTKR(B, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err21 = norm(A2 - A)/normA;
        
        [Q, t3] = randomizedTTKron(B, TTr, 10, 0);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err3 = norm(A2 - A)/normA;
        
        [Q, t31] = randomizedTTKron(B, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err31 = norm(A2 - A)/normA;
        
        [Q, t4] = prandomizedTTeven(B, TTr, 10, 1);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err4 = norm(A2 - A)/normA;
        
        [Q, t41] = prandomizedTTeven(B, TTr, 10, 2);
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        err41 = norm(A2 - A)/normA;
        
        T1(P, sample) = t1;T2(P, sample) = t2;T3(P, sample) = t3;T4(P, sample) = t4;
        T11(P, sample) = t11;T21(P, sample) = t21;T31(P, sample) = t31;T41(P, sample) = t41;
        T5(P, sample) = t5;
        ERR1(P, sample) = err1;ERR2(P, sample) = err2;ERR3(P, sample) = err3;ERR4(P, sample) = err4;
        ERR11(P, sample) = err11;ERR21(P, sample) = err21;ERR31(P, sample) = err31;ERR41(P, sample) = err41;
        ERR5(P, sample) = err5;
    end
    TTr
end

sample = 10;
figure(1)
subplot(121)
semilogy(kk, sum(ERR5, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR1, 2)/sample, ':+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR11, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR2, 2)/sample, ':*', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR21, 2)/sample, '-*', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR3, 2)/sample, ':o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR31, 2)/sample, '-o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR4, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR41, 2)/sample, '-d', 'LineWidth', 2, 'MarkerSize', 10);hold on
set(gca, 'fontsize', 20);
subplot(122)
semilogy(kk, sum(T5, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T1, 2)/sample, ':+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T11, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T2, 2)/sample, ':*', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T21, 2)/sample, '-*', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T3, 2)/sample, ':o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T31, 2)/sample, '-o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T4, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(T41, 2)/sample, '-d', 'LineWidth', 2, 'MarkerSize', 10);hold on
set(gca, 'fontsize', 20);
figure(2)
subplot(121)
semilogy(kk, sum(ERR5, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR1, 2)/sample, ':+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR11, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR2, 2)/sample, ':*', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR21, 2)/sample, '-*', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR3, 2)/sample, ':o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR31, 2)/sample, '-o', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR4, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
semilogy(kk, sum(ERR41, 2)/sample, '-d', 'LineWidth', 2, 'MarkerSize', 10);hold on
set(gca, 'fontsize', 20);
subplot(122)
plot(kk, sum(T5, 2)/sample, '--s', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T1, 2)/sample, ':+', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T11, 2)/sample, '-+', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T2, 2)/sample, ':*', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T21, 2)/sample, '-*', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T3, 2)/sample, ':o', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T31, 2)/sample, '-o', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T4, 2)/sample, ':d', 'LineWidth', 2, 'MarkerSize', 10);hold on
plot(kk, sum(T41, 2)/sample, '-d', 'LineWidth', 2, 'MarkerSize', 10);hold on
set(gca, 'fontsize', 20);
end