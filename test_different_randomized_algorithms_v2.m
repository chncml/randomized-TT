function test_different_randomized_algorithms_v2
%% existing randomized algorithms
clc;clear;nn = [50, 50, 50, 50, 50];
AT1 = [];AT2 = [];AT3 = [];AT4 = [];AT5 = [];
AT6 = [];AT7 = [];AT8 = [];AT9 = [];AT10 = [];
AERR1 = [];AERR2 = [];AERR3 = [];AERR4 = [];AERR5 = [];
AERR6 = [];AERR7 = [];AERR8 = [];AERR9 = [];AERR10 = [];
P = 20;
tic
G = cell(1,5);G{1} = randn(50, P);
G{2} = randn(P, 50, P);G{3} = randn(P, 50, P);
G{4} = randn(P, 50, P);G{5} = randn(P, 50);
B = ttt(tensor(G{1}), tensor(G{2}), [2], [1]);
B = ttt(B, tensor(G{3}), [3], [1]);
B = ttt(B, tensor(G{4}), [4], [1]);
B = ttt(B, tensor(G{5}), [5], [1]);
B = double(B);
toc
beta = [-30 : 5 : 30];TTr = [P, P, P, P];
for k = 1 : length(beta)
    k
    for sample = 1 : 10
        sample
        A = noisegen(B(:),beta(k));
        A = reshape(A, nn);normA = norm(A(:));
        A1 = A;A2 = A;A3 = A;A4 = A;A5 = A;
        A6 = A;A7 = A;A8 = A;A9 = A;A10 = A;
        [Q, T1] = randomizedTTGaussian(A, TTr, 10, 0);
        AT1(k, sample) = T1;
        a1 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a1 = ttt(a1,tensor(Q{3}),[3],[1]);
        a1 = ttt(a1,tensor(Q{4}),[4],[1]);
        a1 = ttt(a1,tensor(Q{5}),[5],[1]);
        AERR1(k, sample) = norm(a1 - tensor(A))/normA;
        
        [Q, T2] = randomizedTTGaussian(A, TTr, 10, 1);
        AT2(k, sample) = T2;
        a2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a2 = ttt(a2,tensor(Q{3}),[3],[1]);
        a2 = ttt(a2,tensor(Q{4}),[4],[1]);
        a2 = ttt(a2,tensor(Q{5}),[5],[1]);
        AERR2(k, sample) = norm(a2 - tensor(A))/normA;
        
        [Q, T3] = randomizedTTKR(A, TTr, 10, 0);
        AT3(k, sample) = T3;
        a3 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a3 = ttt(a3,tensor(Q{3}),[3],[1]);
        a3 = ttt(a3,tensor(Q{4}),[4],[1]);
        a3 = ttt(a3,tensor(Q{5}),[5],[1]);
        AERR3(k, sample) = norm(a3 - tensor(A))/normA;
        
        [Q, T4] = randomizedTTKR(A, TTr, 10, 1);
        AT4(k, sample) = T4;
        a4 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a4 = ttt(a4,tensor(Q{3}),[3],[1]);
        a4 = ttt(a4,tensor(Q{4}),[4],[1]);
        a4 = ttt(a4,tensor(Q{5}),[5],[1]);
        AERR4(k, sample) = norm(a4 - tensor(A))/normA;
        
        [Q, T5] = prandomizedTTeven(A, TTr, 10, 1);
        AT5(k, sample) = T5;
        a5 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a5 = ttt(a5,tensor(Q{3}),[3],[1]);
        a5 = ttt(a5,tensor(Q{4}),[4],[1]);
        a5 = ttt(a5,tensor(Q{5}),[5],[1]);
        AERR5(k, sample) = norm(a5 - tensor(A))/normA;
        
        [Q, T6] = prandomizedTTeven(A, TTr, 10, 2);
        AT6(k, sample) = T6;
        a6 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a6 = ttt(a6,tensor(Q{3}),[3],[1]);
        a6 = ttt(a6,tensor(Q{4}),[4],[1]);
        a6 = ttt(a6,tensor(Q{5}),[5],[1]);
        AERR6(k, sample) = norm(a6 - tensor(A))/normA;
        
%         [Q, T7] = randomizedTTexisting1(A, TTr, 10, 0);AT7(k, sample) = T7;
%         a7 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
%         a7 = ttt(a7,tensor(Q{3}),[3],[1]);
%         a7 = ttt(a7,tensor(Q{4}),[4],[1]);
%         a7 = ttt(a7,tensor(Q{5}),[5],[1]);
%         AERR7(k, sample) = norm(a7 - tensor(A))/normA;
%         
        [Q, t8] = TTSVD1(A, TTr);AT8(k, sample) = t8;
        A2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        A2 = ttt(A2,tensor(Q{3}),[3],[1]);
        A2 = ttt(A2,tensor(Q{4}),[4],[1]);
        A2 = ttt(A2,tensor(Q{5}),[5],[1]);
        AERR8(k, sample) = norm(A2 - tensor(A))/normA;
        
        [Q, T9] = randomizedTTexisting2(A, TTr, 10, 0);
        AT9(k, sample) = T9;
        a9 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a9 = ttt(a9,tensor(Q{3}),[3],[1]);
        a9 = ttt(a9,tensor(Q{4}),[4],[1]);
        a9 = ttt(a9,tensor(Q{5}),[5],[1]);
        AERR9(k, sample) = norm(a9 - tensor(A))/normA;
        
        [Q, T10] = randomizedTTexisting2(A, TTr, 10, 1);
        AT10(k, sample) = T10;
        a10 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a10 = ttt(a10,tensor(Q{3}),[3],[1]);
        a10 = ttt(a10,tensor(Q{4}),[4],[1]);
        a10 = ttt(a10,tensor(Q{5}),[5],[1]);
        AERR10(k, sample) = norm(a10 - tensor(A))/normA;
    end
end
sample = 10;
figure(1);
subplot(1,2,1);
plot(beta, 1 - sum(AERR8,2)/sample,'s','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR9,2)/sample,'-*','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR10,2)/sample,':*','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR1,2)/sample,'-+','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR2,2)/sample,':+','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR3,2)/sample,'-o','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR4,2)/sample,':o','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR5,2)/sample,'-d','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR6,2)/sample,':d','LineWidth',2,'MarkerSize',10);hold on;

xlabel('P (number of indices)');
ylabel('RLNE');
set(gca, 'fontsize', 20);
subplot(1,2,2);
plot(beta, sum(AT8,2)/sample,'s','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT9,2)/sample,'-*','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT10,2)/sample,':*','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT1,2)/sample,'-+','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT2,2)/sample,':+','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT3,2)/sample,'-o','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT4,2)/sample,':o','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT5,2)/sample,'-d','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, sum(AT6,2)/sample,':d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('CPU time');
set(gca, 'fontsize', 20);
figure(2);
subplot(1,2,1);
plot(beta, 1 - sum(AERR8,2)/sample,'s','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR9,2)/sample,'-*','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR10,2)/sample,':*','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR1,2)/sample,'-+','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR2,2)/sample,':+','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR3,2)/sample,'-o','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR4,2)/sample,':o','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR5,2)/sample,'-d','LineWidth',2,'MarkerSize',10);hold on;
plot(beta, 1 - sum(AERR6,2)/sample,':d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('RLNE');
set(gca, 'fontsize', 20);
subplot(1,2,2);
semilogy(beta, sum(AT8,2)/sample,'s','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT9,2)/sample,'-*','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT10,2)/sample,':*','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT1,2)/sample,'-+','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT2,2)/sample,':+','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT3,2)/sample,'-o','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT4,2)/sample,':o','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT5,2)/sample,'-d','LineWidth',2,'MarkerSize',10);hold on;
semilogy(beta, sum(AT6,2)/sample,':d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('CPU time');
set(gca, 'fontsize', 20);
end