clear; close all; clc;

syms k__0 k__1 k__2 lambda q dn

% %% 1.1.2 square-lattice eq32
% A = [-k__0 .* lambda + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda -(k__1 .* lambda) - exp(-i .* q) .* k__2 ./ lambda; -lambda .* exp(i .* q) .* k__2 - (k__1 ./ lambda) -k__0 .* lambda + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda;];
% 
% q = linspace(-pi, pi, 100);
% k__0 = 1;
% k__1 = 2;
% k__2 = 3;
% 
% lambda = double(subs(solve(det(A) == 0, lambda)));
% 
% figure(1)
% plot(q, abs(lambda(1,:)))
% title("abs(lambda)")
% xlabel("q")
% ylabel("abs(lambda)")
% axis([-pi pi 0 1])
% 
% figure(2)
% plot(q, angle(lambda(1,:)))
% title("arg(lambda)")
% xlabel("q")
% ylabel("arg(lambda)")

% %% 1.2.1 square lattice eq53
% A = [-k__0 .* lambda + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda -(k__2 .* lambda .^ 2) - exp(-i .* q) .* k__1 ./ (lambda .^ 2); -(lambda .^ 2) .* exp(i .* q) .* k__1 - (k__2 ./ lambda .^ 2) -k__0 .* lambda + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda;];
% 
% q = linspace(-pi, pi, 100);
% k__0 = 1;
% k__1 = 3;
% k__2 = 2; % 这里k__2.k__3和别的是相反的
% 
% lambda = double(subs(solve(det(A) == 0, lambda)));
% 
% figure(1)
% plot(q, abs(lambda))
% title("abs(lambda)")
% xlabel("q")
% ylabel("abs(lambda)")
% axis([-pi pi 0 1])
% 
% figure(2)
% plot(q, angle(lambda))
% title("arg(lambda)")
% xlabel("q")
% ylabel("arg(lambda)")

% %% 1.3.1 square lattice eq72
% A = [-k__0 .* lambda + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda -(k__1 .* lambda .^ 3) - exp(-i .* q) .* k__2 ./ (lambda .^ 3); -(lambda .^ 3) .* exp(q .* i) .* k__2 - (k__1 ./ lambda .^ 3) -k__0 .* lambda + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda;];
% 
% q = linspace(-pi, pi, 100);
% k__0 = 1;
% k__1 = 2;
% k__2 = 3;
% 
% lambda = double(subs(solve(det(A) == 0, lambda)));
% 
% figure(1)
% plot(q, abs(lambda))
% title("abs(lambda)")
% xlabel("q")
% ylabel("abs(lambda)")
% axis([-pi pi 0 1])
% 
% figure(2)
% plot(q, angle(lambda))
% title("arg(lambda)")
% xlabel("q")
% ylabel("arg(lambda)")

% %% 1.3.1拓展 square lattice eq91
% % A = [-lambda .* k__0 + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda -(lambda .^ (-dn)) .* exp(-i .* q) .* k__2 - (lambda .^ dn .* k__1); -(lambda .^ dn) .* exp(q .* i) .* k__2 - (k__1 .* lambda .^ (-dn)) -lambda .* k__0 + 2 .* k__0 + k__1 + k__2 - k__0 ./ lambda;];
% syms k__0 k__1 k__2 lambda q dn
% 
% dn = 1;
% k__0 = 1;
% k__1 = 2;
% k__2 = 3;
% 
% lambda = solve(k__0 .^ 2 .* lambda .^ 4 .* exp(q .* i) .* (lambda .^ dn) .^ 2 - 4 .* k__0 .^ 2 .* lambda .^ 3 .* exp(q .* i) .* (lambda .^ dn) .^ 2 - 2 .* k__0 .* k__1 .* lambda .^ 3 .* exp(q .* i) .* (lambda .^ dn) .^ 2 - 2 .* k__0 .* k__2 .* lambda .^ 3 .* exp(q .* i) .* (lambda .^ dn) .^ 2 - k__1 .* k__2 .* exp(q .* i) .^ 2 .* (lambda .^ dn) .^ 4 .* lambda .^ 2 + 6 .* k__0 .^ 2 .* lambda .^ 2 .* exp(q .* i) .* (lambda .^ dn) .^ 2 + 4 .* k__0 .* k__1 .* lambda .^ 2 .* exp(q .* i) .* (lambda .^ dn) .^ 2 + 4 .* k__0 .* k__2 .* lambda .^ 2 .* exp(q .* i) .* (lambda .^ dn) .^ 2 + 2 .* k__1 .* k__2 .* lambda .^ 2 .* exp(q .* i) .* (lambda .^ dn) .^ 2 - k__1 .* k__2 .* lambda .^ 2 - 4 .* k__0 .^ 2 .* lambda .* exp(q .* i) .* (lambda .^ dn) .^ 2 - 2 .* k__0 .* k__1 .* lambda .* exp(q .* i) .* (lambda .^ dn) .^ 2 - 2 .* k__0 .* k__2 .* lambda .* exp(q .* i) .* (lambda .^ dn) .^ 2 + exp(q .* i) .* (lambda .^ dn) .^ 2 .* k__0 .^ 2 == 0, lambda);
% 
% q = linspace(-pi, pi, 1000);
% 
% lambda = double(subs(lambda));
% 
% figure(1)
% plot(q, abs(lambda))
% title("abs(lambda)")
% xlabel("q")
% ylabel("abs(lambda)")
% axis([-pi pi 0 1])
% 
% figure(2)
% plot(q, angle(lambda))
% title("arg(lambda)")
% xlabel("q")
% ylabel("arg(lambda)")
% 
% figure(3)
% plot3(real(lambda(3,:)), imag(lambda(3,:)) , q)
% hold on
% plot3(real(lambda(4,:)), imag(lambda(4,:)) , q)
% hold on
% title("Re(lambda)-Im(lambda)-q")
% xlabel("Re(lambda)")
% ylabel("Im(lambda)")
% zlabel("q")

% %% 1.5 square lattice eq188
% A = [(2 .* k__0) - k__0 .* lambda .^ (0.1e1 ./ 0.5e1) .* exp(-0.2e1 ./ 0.5e1.*i .* q) - k__0 .* lambda .^ (-0.1e1 ./ 0.5e1) .* exp(0.2e1 ./ 0.5e1.*i .* q) + k__1 + k__2 -k__1 - k__2 .* lambda .^ (-0.3e1 ./ 0.5e1) .* exp(0.1e1 ./ 0.5e1.*i .* q); -k__2 .* lambda .^ (0.3e1 ./ 0.5e1) .* exp(-0.1e1 ./ 0.5e1.*i .* q) - k__1 (2 .* k__0) - k__0 .* lambda .^ (0.1e1 ./ 0.5e1) .* exp(-0.2e1 ./ 0.5e1.*i .* q) - k__0 .* lambda .^ (-0.1e1 ./ 0.5e1) .* exp(0.2e1 ./ 0.5e1.*i .* q) + k__1 + k__2;];
% 
% k__0 = 1;
% k__1 = 2;
% k__2 = 2;
% 
% eq = lambda .^ 6 .* exp(-q .* i ./ 5) .^ 2 .* k__1 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 - lambda .^ 5 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 4 .* k__0 .^ 2 .* exp(-q .* i ./ 5) + (4 .* k__0 .^ 2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 3 .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__1 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 3 .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 3 .* exp(-q .* i ./ 5)) .* lambda .^ 4 + (-6 .* k__0 .^ 2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 .* exp(-q .* i ./ 5) - 4 .* k__0 .* k__1 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 .* exp(-q .* i ./ 5) - 4 .* k__0 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 .* exp(-q .* i ./ 5) - 2 .* exp(-q .* i ./ 5) .* k__1 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2) .* lambda .^ 3 + (4 .* k__0 .^ 2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__1 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .* exp(-q .* i ./ 5)) .* lambda .^ 2 - lambda .* k__0 .^ 2 .* exp(-q .* i ./ 5) + k__1 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 == 0;
% lambda = solve(eq, lambda);
% 
% q = linspace(-pi, pi, 100);
% lambda = double(subs(lambda));
% 
% figure(1)
% plot(q, abs(lambda))
% title("abs(lambda)")
% xlabel("q")
% ylabel("abs(lambda)")
% axis([-pi pi 0 1])
% 
% figure(2)
% plot(q, angle(lambda))
% title("arg(lambda)")
% xlabel("q")
% ylabel("arg(lambda)")

%% 1.6 square lattice eq220
k__0 = 1;
k__1 = 2;
k__2 = 3;

% % 不能直接用matlab求解用
% eq = lambda .^ 6 .* exp(-q .* i ./ 5) .^ 2 .* k__1 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 - lambda .^ 5 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 4 .* k__0 .^ 2 .* exp(-q .* i ./ 5) + (4 .* k__0 .^ 2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 3 .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__1 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 3 .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 3 .* exp(-q .* i ./ 5)) .* lambda .^ 4 + (-6 .* k__0 .^ 2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 .* exp(-q .* i ./ 5) - 4 .* k__0 .* k__1 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 .* exp(-q .* i ./ 5) - 4 .* k__0 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 .* exp(-q .* i ./ 5) - 2 .* exp(-q .* i ./ 5) .* k__1 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2) .* lambda .^ 3 + (4 .* k__0 .^ 2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__1 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .* exp(-q .* i ./ 5) + 2 .* k__0 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .* exp(-q .* i ./ 5)) .* lambda .^ 2 - lambda .* k__0 .^ 2 .* exp(-q .* i ./ 5) + k__1 .* k__2 .* exp(-(0.2e1 ./ 0.5e1) .* q .* i) .^ 2 == 0;
% lambda = solve(eq, lambda);
% 
% q = linspace(-pi, pi, 100);
% lambda = double(subs(lambda));

% 可以直接用matlab求解
A = [2 * k__0 - k__0 / lambda - k__0 * lambda + k__1 + k__2 -k__1 - k__2 * lambda / exp(q * i); -k__1 - k__2 / lambda * exp(q * i) 2 * k__0 - k__0 / lambda - k__0 * lambda + k__1 + k__2;];

lambda = solve(det(A) == 0, lambda);
q = linspace(-5*pi, 5*pi, 101);
lambda = double(subs(lambda));

figure(1)
plot(q, abs(lambda))
title("abs(lambda)")
xlabel("q")
ylabel("abs(lambda)")
axis([-pi pi 0 1])

figure(2)
for i = 1:max(size(lambda(:,1)))
    plot(q, angle(lambda(i,:)))
    hold on
end
title("arg(lambda)")
xlabel("q")
ylabel("arg(lambda)")

figure(3)
for i = 1:max(size(lambda(:,1)))
    plot(real(lambda(i,:)), imag(lambda(i,:)))
    hold on
end
title("Im(lambda)-Re(lambda)")
xlabel("Re(lambda)")
ylabel("Im(lambda)")

figure(4)
plot(q, real(lambda), '-', q, imag(lambda), '--')
title("q-Re(lambda)[Im(lambda)]")
xlabel("q")
ylabel("Re(lambda)[Im(lambda)]")

