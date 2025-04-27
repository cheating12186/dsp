clear all;
close all;
clc;

% Define the range of n
n = -50:50;  % Time index
omega = -4*pi:0.01:4*pi;  % Frequency range for DTFT
X = zeros(size(omega));   % Initialize X(omega)

% Assuming X(omega) is known (you could compute it using DTFT of any signal)
% Here we'll use a dummy X(omega) for demonstration purposes (e.g., a simple sine wave)
X = 1 ./ (1 + (omega/pi).^2);  % Example DTFT (Lorentzian function)

% Inverse DTFT computation
x_rec = zeros(size(n));  % Initialize the reconstructed signal

% Discretize the inverse DTFT integral
for nn = 1:length(n)
    integral_sum = 0;
    for k = 1:length(omega)
        integral_sum = integral_sum + X(k) * exp(1j * omega(k) * n(nn));  % Numerical sum for IDTFT
    end
    x_rec(nn) = real(integral_sum) / (2 * pi);  % Normalize and take the real part
end

% Plot the results
figure;

% Plot the original signal
subplot(3,2,1:2);
stem(n, x_rec, 'filled');
xlabel('n');
ylabel('x[n]');
title('Reconstructed Signal from IDTFT');

% Plot the magnitude of X(omega)
subplot(3,2,3:4);
plot(omega/pi, abs(X)); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude of X(ω)");

% Plot the phase of X(omega)
subplot(3,2,5:6);
plot(omega/pi, angle(X)); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase of X(ω)");
