clc;
clear all;
close all;

%% Input
N = input("Enter the value of N: ");
x = input("Enter the input sequence x[n]: ");
n = 0:1:N-1;

%% FFT
X_fft = fft(x);

%% Manual IDFT
x_manual = zeros(N,1);
for nn = 0:N-1
    for k = 0:N-1
        x_manual(nn+1) = x_manual(nn+1) + (X_fft(k+1)*exp(j*2*pi*nn*k/N));
    end
end

%% IDFT
x_manual = abs(x_manual/N);
x_ifft = ifft(X_fft);

%% Plots
subplot(3,1,1);
stem(n,x);

subplot(3,1,2);
stem(n,x_manual);

subplot(3,1,3);
stem(n,x_ifft);
