clc;
close all;
clear all;
f1 = 1000; f2 = 1900;
fmax = max(f1,f2);
T = 1/min(f1,f2);
t = 0:0.01*T:2*T;
x = cos(2*pi*f1*t)+cos(2*pi*f2*t);
subplot(2,3,[1,2,3]);
plot(t,x);
xlabel('t');
ylabel('x(t)');
title("Continuous Timed Signal");
grid on;

%% over sampling condition
fs1 = 10*fmax;
n1 = 0:1/fs1:2*T;
x1 = cos(2*pi*f1*n1)+cos(2*pi*f2*n1);
subplot(2,3,4);
stem(n1,x1);
hold on;
plot(n1,x1,"red",LineStyle="-");
hold off;
title("fs > 2fm signal");
grid on;

%% nyquist condition
fs2 = 2*fmax;
n2 = 0:1/fs2:2*T;
x2 = cos(2*pi*f1*n2)+cos(2*pi*f2*n2);
subplot(2,3,5);
stem(n2,x2);
hold on;
plot(n2,x2,"red",LineStyle="-");
hold off;
title("fs = 2fm signal");
grid on;

%% under sampled
fs3 = fmax;
n3 = 0:1/fs3:2*T;
x3 = cos(2*pi*f1*n3)+cos(2*pi*f2*n3);
subplot(2,3,6);
stem(n3,x3);
hold on;
plot(n3,x3,"red",LineStyle="-");
hold off;
title("fs = 2fm signal");
grid on;