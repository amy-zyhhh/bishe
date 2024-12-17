clear all; close all; clc;

syms k__0 k__1 k__2 lambda q

A = [-k__0 * lambda + 2 * k__0 + k__1 + k__2 - k__0 / lambda -(k__1 * lambda) - exp(-i * q) * k__2 / lambda; -exp(i * q) * lambda * k__2 - (k__1 / lambda) -k__0 * lambda + 2 * k__0 + k__1 + k__2 - k__0 / lambda;];

k__0 = 1;
k__1 = 3;
k__2 = 2;
subs(A);
lambda = solve(det(A)==0, lambda);
q = linspace(-pi, pi, 100);
lambda = double(subs(lambda));

figure(2)
plot(q, abs(lambda))
axis([-pi pi 0 1])
title('abs(lambda)')
ylabel('abs(lambda)')
xlabel('q')

figure(3)
plot(q, imag(lambda))
axis([-pi pi 0 1])
title('arg(lambda)')
ylabel('imag(lambda)')
xlabel('q')