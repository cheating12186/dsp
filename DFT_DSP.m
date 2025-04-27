clc;
clear all;
close all;

% Input
N = input("Enter the value of N: ");
x = input("Enter the input sequence x[n]: ");
n = 0:1:N-1;

%% Manual Computation of DFT
X = zeros(N, 1);  % Initialize X array for DFT

for k = 0:N-1
    for nn = 0:N-1
        X(k+1) = X(k+1) + (x(nn+1) * exp(-1j * 2 * pi * nn * k / N));
    end
end

% Get magnitude of X (DFT)
X_mag = abs(X);

%% FFT Command
X_fft = fft(x);  % Apply fft directly to the input sequence x[n]

%% Manual IDFT
x_manual = zeros(N, 1);
for nn = 0:N-1
    for k = 0:N-1
        x_manual(nn+1) = x_manual(nn+1) + (X(k+1) * exp(1j * 2 * pi * nn * k / N));
    end
end

% Scale the manual IDFT result by 1/N
x_manual = x_manual / N;

%% IFFT
x_ifft = ifft(X_fft);

%% Plots
subplot(3, 2, [1,2]);
stem(n, abs(x));   % Plot the input sequence
grid on;
xlabel("n");
ylabel("x[n]");
title("Sequence x[n]");

subplot(3, 2, 3);
stem(n, X_mag);   % Plot the magnitude of the manually computed DFT
grid on;
xlabel("Frequency");
ylabel("X(k)");
title("Manual Computation of DFT");

subplot(3, 2, 5);
stem(n, abs(X_fft));   % Plot the magnitude of the FFT
grid on;
xlabel("Frequency");
ylabel("X(k)");
title("FFT(x[n])");

subplot(3, 2, 4);
stem(n, real(x_manual));   % Plot the real part of the manual IDFT
grid on;
xlabel("n");
ylabel("x[n]");
title("Manual IDFT(x[n])");

subplot(3, 2, 6);
stem(n, real(x_ifft));   % Plot the real part of the IFFT result
grid on;
xlabel("n");
ylabel("x[n]");
title("IFFT(x[n])");

