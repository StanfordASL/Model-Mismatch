clear all; close all; clc;

%% Get approximation

%choose function here (in vectorized operation form)
fnc = @(x) 1./x;

%number of points to use to compute best Chebyshev approximation
n = 100;

%order of approximation
m = 3;

%interval of approximation
intvl = [3,10];
bma=0.5*(intvl(2)-intvl(1));
bpa=0.5*(intvl(2)+intvl(1));

c = get_cheby_approx(fnc,intvl,n,m);
c(1) = 0.5*c(1);

%% Check approximation
k = 100;

x = linspace(intvl(1),intvl(2),k+1);
% x = linspace(-pi/4,pi/4,k+1);
y = (x-bpa)/bma;

[T,~] = compute_cheby(k,m,y);

f_val = fnc(x);
f_app = (c*T);

%% Plot
figure()
subplot(2,1,1)
plot(x,f_val,'r-','linewidth',2);
hold on
plot(x,f_app,'b-','linewidth',2);
grid on
legend('actual','approx');

subplot(2,1,2)
plot(x,f_val-f_app,'r-','linewidth',2);
grid on
title('error');





