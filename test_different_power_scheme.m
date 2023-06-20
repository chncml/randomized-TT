function test_different_power_scheme
clc;clear;nn = [50, 50, 50, 50, 50];
AT1 = [];AT2 = [];AT3 = [];AT4 = [];
AT5 = [];AT6 = [];AT7 = [];AT8 = [];AT9 = [];AT10 = [];
AERR1 = [];AERR2 = [];AERR3 = [];AERR4 = [];
AERR5 = [];AERR6 = [];AERR7 = [];AERR8 = [];AERR9 = [];AERR10 = [];
P = 20;
tic
G = cell(1,5);G{1} = randn(50, P);
G{2} = randn(P, 50, P);G{3} = randn(P, 50, P);
G{4} = randn(P, 50, P);G{5} = randn(P, 50);
B = ttt(tensor(G{1}), tensor(G{2}), [2], [1]);
B = ttt(B, tensor(G{3}), [3], [1]);
B = ttt(B, tensor(G{4}), [4], [1]);
B = ttt(B, tensor(G{5}), [5], [1]);
gamma = 1e-4;
A = double(B);normA = norm(A(:));
A = A + gamma * normA/sqrt(sum(nn.^2)) * randn(nn);
toc
kk = [2 : 2 : P];
normA = norm(A(:));
A1 = A;A2 = A;A3 = A;A4 = A;A5 = A;A6 = A;A7 = A;A8 = A;A9 = A;A10 = A;
for k = 1:length(kk)
    k
    TTr = 2 * [k, k, k, k];
    for sample = 1 : 10
        sample
        [Q, T1] = randomizedTTGaussian(A, TTr, 10, 0);AT1(k, sample) = T1;
        a1 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a1 = ttt(a1,tensor(Q{3}),[3],[1]);
        a1 = ttt(a1,tensor(Q{4}),[4],[1]);
        a1 = ttt(a1,tensor(Q{5}),[5],[1]);
        AERR1(k, sample) = norm(a1 - tensor(A))/normA;
        
        [Q, T2] = randomizedTTGaussian(A, TTr, 10, 1);AT2(k, sample) = T2;
        a2 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a2 = ttt(a2,tensor(Q{3}),[3],[1]);
        a2 = ttt(a2,tensor(Q{4}),[4],[1]);
        a2 = ttt(a2,tensor(Q{5}),[5],[1]);
        AERR2(k, sample) = norm(a2 - tensor(A))/normA;
        
        [Q, T3] = randomizedTTGaussian(A, TTr, 10, 2);AT3(k, sample) = T3;
        a3 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a3 = ttt(a3,tensor(Q{3}),[3],[1]);
        a3 = ttt(a3,tensor(Q{4}),[4],[1]);
        a3 = ttt(a3,tensor(Q{5}),[5],[1]);
        AERR3(k, sample) = norm(a3 - tensor(A))/normA;
        
        [Q, T4] = randomizedTTGaussian(A, TTr, 10, 3);AT4(k, sample) = T4;
        a4 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a4 = ttt(a4,tensor(Q{3}),[3],[1]);
        a4 = ttt(a4,tensor(Q{4}),[4],[1]);
        a4 = ttt(a4,tensor(Q{5}),[5],[1]);
        AERR4(k, sample) = norm(a4 - tensor(A))/normA;
        
        [Q, T5] = randomizedTTKR(A, TTr, 10, 0);AT5(k, sample) = T5;
        a5 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a5 = ttt(a5,tensor(Q{3}),[3],[1]);
        a5 = ttt(a5,tensor(Q{4}),[4],[1]);
        a5 = ttt(a5,tensor(Q{5}),[5],[1]);
        AERR5(k, sample) = norm(a5 - tensor(A))/normA;
        
        [Q, T6] = randomizedTTKR(A, TTr, 10, 1);AT6(k, sample) = T6;
        a6 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a6 = ttt(a6,tensor(Q{3}),[3],[1]);
        a6 = ttt(a6,tensor(Q{4}),[4],[1]);
        a6 = ttt(a6,tensor(Q{5}),[5],[1]);
        AERR6(k, sample) = norm(a6 - tensor(A))/normA;
        
        [Q, T7] = randomizedTTKR(A, TTr, 10, 2);AT7(k, sample) = T7;
        a7 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a7 = ttt(a7,tensor(Q{3}),[3],[1]);
        a7 = ttt(a7,tensor(Q{4}),[4],[1]);
        a7 = ttt(a7,tensor(Q{5}),[5],[1]);
        AERR7(k, sample) = norm(a7 - tensor(A))/normA;
        
        [Q, T8] = randomizedTTKR(A, TTr, 10, 3);AT8(k, sample) = T8;
        a8 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a8 = ttt(a8,tensor(Q{3}),[3],[1]);
        a8 = ttt(a8,tensor(Q{4}),[4],[1]);
        a8 = ttt(a8,tensor(Q{5}),[5],[1]);
        AERR8(k, sample) = norm(a8 - tensor(A))/normA;
        
        [Q, T9] = prandomizedTTeven(A, TTr, 10, 1);AT9(k, sample) = T9;
        a9 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a9 = ttt(a9,tensor(Q{3}),[3],[1]);
        a9 = ttt(a9,tensor(Q{4}),[4],[1]);
        a9 = ttt(a9,tensor(Q{5}),[5],[1]);
        AERR9(k, sample) = norm(a9 - tensor(A))/normA;
        
        [Q, T10] = prandomizedTTeven(A, TTr, 10, 2);AT10(k, sample) = T10;
        a10 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a10 = ttt(a10,tensor(Q{3}),[3],[1]);
        a10 = ttt(a10,tensor(Q{4}),[4],[1]);
        a10 = ttt(a10,tensor(Q{5}),[5],[1]);
        AERR10(k, sample) = norm(a10 - tensor(A))/normA;
        
        [Q, T11] = prandomizedTTeven(A, TTr, 10, 3);AT11(k, sample) = T11;
        a11 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
        a11 = ttt(a11,tensor(Q{3}),[3],[1]);
        a11 = ttt(a11,tensor(Q{4}),[4],[1]);
        a11 = ttt(a11,tensor(Q{5}),[5],[1]);
        AERR11(k, sample) = norm(a11 - tensor(A))/normA;
        
%         [Q, T12] = prandomizedTTeven(A, TTr, 10, 4);AT12(k, sample) = T12;
%         a12 = ttt(tensor(Q{1}),tensor(Q{2}),[2],[1]);
%         a12 = ttt(a12,tensor(Q{3}),[3],[1]);
%         a12 = ttt(a12,tensor(Q{4}),[4],[1]);
%         a12 = ttt(a12,tensor(Q{5}),[5],[1]);
%         AERR12(k, sample) = norm(a12 - tensor(A))/normA;
    end
end
sample = 10;
figure(1);
subplot(1,2,1);
semilogy(kk,sum(AERR1,2)/sample,'-+','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR2,2)/sample,'-o','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR3,2)/sample,'-x','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR4,2)/sample,'-d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('RLNE');
set(gca, 'fontsize', 20);
subplot(1,2,2);
plot(kk,sum(AT1,2)/sample,'-+','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT2,2)/sample,'-o','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT3,2)/sample,'-x','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT4,2)/sample,'-d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('CPU time');
set(gca, 'fontsize', 20);
figure(2);
subplot(1,2,1);
semilogy(kk,sum(AERR5,2)/sample,'--+','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR6,2)/sample,'--o','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR7,2)/sample,'--x','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR8,2)/sample,'--d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('RLNE');
set(gca, 'fontsize', 20);
subplot(1,2,2);
plot(kk,sum(AT5,2)/sample,'--+','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT6,2)/sample,'--o','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT7,2)/sample,'--x','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT8,2)/sample,'--d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('CPU time');
set(gca, 'fontsize', 20);
figure(3);
subplot(1,2,1);
semilogy(kk,sum(AERR9,2)/sample,':+','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR10,2)/sample,':o','LineWidth',2,'MarkerSize',10);hold on;
semilogy(kk,sum(AERR11,2)/sample,':x','LineWidth',2,'MarkerSize',10);hold on;
% semilogy(kk,sum(AERR12,2)/sample,':d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('RLNE');
set(gca, 'fontsize', 20);
subplot(1,2,2);
plot(kk,sum(AT9,2)/sample,':+','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT10,2)/sample,':o','LineWidth',2,'MarkerSize',10);hold on;
plot(kk,sum(AT11,2)/sample,':x','LineWidth',2,'MarkerSize',10);hold on;
% plot(kk,sum(AT12,2)/sample,':d','LineWidth',2,'MarkerSize',10);hold on;
xlabel('P (number of indices)');
ylabel('CPU time');
set(gca, 'fontsize', 20);
end