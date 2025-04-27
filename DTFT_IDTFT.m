% To find the DTFT and IDTFT of Discrete Signals
clear all;
close all;
clc;

n = -50:50;

%% Unit Impulse Function
figure("Name", "Unit Impulse DTFT");
for ii = 1:length(n)
    if n(ii) == 0
        x(ii) = 1;
    else
        x(ii) = 0;
    end
end
subplot(3,2,1:2);
stem(n,x);
xlabel('n');
ylabel('x[n]');
title("Discrete Impulse Sequence");
omega = -4*pi:0.01:4*pi;
Xtemp = [];
X = [];
for a = 1:length(omega)
    for b = 1:length(n)
        Xtemp = [Xtemp x(b)*exp(-j*omega(a)*n(b))];
    end
    Xtemp = sum(Xtemp);
    X = [X Xtemp];
    Xtemp = [];
end
magx = abs(X);
subplot(3,2,3:4);
plot(omega/pi, magx); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude Plot");
angx = angle(X);
subplot(3,2,5:6);
plot(omega/pi, angx); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase Plot");

figure("Name", "Unit Step DTFT");
%% Unit Step Sequence
uT = zeros(size(n));
uT(n >= 0) = 1; 
subplot(3,2,1:2);
stem(n,uT);
xlabel('n');
ylabel('x[n]');
title("Discrete Unit Sequence");
omega = -4*pi:0.01:4*pi;
Xtemp = [];
X = [];
for a = 1:length(omega)
    for b = 1:length(n)
        Xtemp = [Xtemp uT(b)*exp(-j*omega(a)*n(b))];
    end
    Xtemp = sum(Xtemp);
    X = [X Xtemp];
    Xtemp = [];
end
magx = abs(X);
subplot(3,2,3:4);
plot(omega/pi, magx); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude Plot");
angx = angle(X);
subplot(3,2,5:6);
plot(omega/pi, angx); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase Plot");

figure("Name", "Exponential DTFT");
%% Exponential Sequence
eT = zeros(size(n));
eT(n >= 0) = (0.5).^n(n >= 0);
subplot(3,2,1:2);
stem(n,eT);
xlabel('n');
ylabel('x[n]');
title("Discrete Exponential Sequence");
omega = -4*pi:0.01:4*pi;
Xtemp = [];
X = [];
for a = 1:length(omega)
    for b = 1:length(n)
        Xtemp = [Xtemp eT(b)*exp(-j*omega(a)*n(b))];
    end
    Xtemp = sum(Xtemp);
    X = [X Xtemp];
    Xtemp = [];
end
magx = abs(X);
subplot(3,2,3:4);
plot(omega/pi, magx); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude Plot");
angx = angle(X);
subplot(3,2,5:6);
plot(omega/pi, angx); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase Plot");

figure("Name", "Sine DTFT");
%% Sine Sequence
Sin = zeros(size(n)); M = 1; N = 10;
Sin = 1*sin(2*pi*M/N*n);
subplot(3,2,1:2);
stem(n,Sin);
xlabel('n');
ylabel('x[n]');
title("Discrete Sine Sequence");
omega = -4*pi:0.01:4*pi;
Xtemp = [];
X = [];
for a = 1:length(omega)
    for b = 1:length(n)
        Xtemp = [Xtemp Sin(b)*exp(-j*omega(a)*n(b))];
    end
    Xtemp = sum(Xtemp);
    X = [X Xtemp];
    Xtemp = [];
end
magx = abs(X);
subplot(3,2,3:4);
plot(omega/pi, magx); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude Plot");
angx = angle(X);
subplot(3,2,5:6);
plot(omega/pi, angx); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase Plot");

figure("Name", "Rectangle DTFT");
%% Rect Sequence
Ret = zeros(size(n)); 
Ret(n > -5 & n < 5) = 1;
subplot(3,2,1:2);
stem(n,Ret);
xlabel('n');
ylabel('x[n]');
title("Discrete Rectangular Sequence");
omega = -4*pi:0.01:4*pi;
Xtemp = [];
X = [];
for a = 1:length(omega)
    for b = 1:length(n)
        Xtemp = [Xtemp Ret(b)*exp(-j*omega(a)*n(b))];
    end
    Xtemp = sum(Xtemp);
    X = [X Xtemp];
    Xtemp = [];
end
magx = abs(X);
subplot(3,2,3:4);
plot(omega/pi, magx); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude Plot");
angx = angle(X);
subplot(3,2,5:6);
plot(omega/pi, angx); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase Plot");

figure("Name", "Sinc DTFT");
%% Sinc Sequence
Sinc = zeros(size(n)); 
Sinc(n == 0) = 1;  % Define Sinc(0) = 1 to avoid division by zero at n=0
Sinc(n ~= 0) = sin(pi*n(n ~= 0)/5) ./ (pi*n(n ~= 0)/5); % Correct Sinc definition
subplot(3,2,1:2);
stem(n,Sinc);
xlabel('n');
ylabel('x[n]');
title("Discrete Sinc Sequence");
omega = -4*pi:0.01:4*pi;
Xtemp = [];
X = [];
for a = 1:length(omega)
    for b = 1:length(n)
        Xtemp = [Xtemp Sinc(b)*exp(-j*omega(a)*n(b))];
    end
    Xtemp = sum(Xtemp);
    X = [X Xtemp];
    Xtemp = [];
end
magx = abs(X);
subplot(3,2,3:4);
plot(omega/pi, magx); grid on;
xlabel("frequency in pi");
ylabel("|X(ω)|");
title("Magnitude Plot");
angx = angle(X);
subplot(3,2,5:6);
plot(omega/pi, angx); grid on;
xlabel("frequency in pi");
ylabel("∠X(ω)");
title("Phase Plot");

%% Practical 6:
% Prove linearity property
% Time shifting
% Time reversal
% Parsevals 