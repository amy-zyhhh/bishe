% clear; close all; clc;

flag = 0;
a = 1;
b = 0;
c = 2;
d = 2;
dr = 0;
ds = 1;
k__0 = 1;
k__1 = 6;
k__2 = 6;
syms lambda q
% % 用det(A)求解
% detA = ((lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 4 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 ^ 2 - 4 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 3 * exp(1i * q / (a * d - b * c) * b) * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 ^ 2 - 2 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 3 * exp(1i * q / (a * d - b * c) * b) * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 * k__1 - 2 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 3 * exp(1i * q / (a * d - b * c) * b) * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 * k__2 - (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 * exp(1i * q / (a * d - b * c) * b) ^ 2 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 4 * k__1 * k__2 + 6 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 * exp(1i * q / (a * d - b * c) * b) ^ 2 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 ^ 2 + 4 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 * exp(1i * q / (a * d - b * c) * b) ^ 2 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 * k__1 + 4 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 * exp(1i * q / (a * d - b * c) * b) ^ 2 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 * k__2 + 2 * (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 * exp(1i * q / (a * d - b * c) * b) ^ 2 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__1 * k__2 - (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 * exp(1i * q / (a * d - b * c) * b) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 4 * k__1 * k__2 - 4 * lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) ^ 3 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 ^ 2 - 2 * lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) ^ 3 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 * k__1 - 2 * lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) ^ 3 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 * k__2 + exp(1i * q / (a * d - b * c) * b) ^ 4 * (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 * exp(1i * q / (a * d - b * c) * a) ^ 2 * k__0 ^ 2) / (lambda ^ (0.1e1 / (a * d - b * c) * d)) ^ 2 / exp(1i * q / (a * d - b * c) * b) ^ 2 / (lambda ^ (0.1e1 / (a * d - b * c) * c)) ^ 2 / exp(1i * q / (a * d - b * c) * a) ^ 2;
% lambda = solve(detA == 0, lambda);
% q = linspace(-pi, pi, 101);
% lambda = double(subs(lambda));

% 直接用A求解
A = [(2 * k__0) - k__0 / lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) - k__0 * lambda ^ (0.1e1 / (a * d - b * c) * d) / exp(1i * q / (a * d - b * c) * b) + k__1 + k__2 -k__1 * lambda ^ (0.1e1 / (a * d - b * c) * c * ds) / lambda ^ (0.1e1 / (a * d - b * c) * d * dr) / lambda ^ (0.1e1 / (a * d - b * c) * c) / exp(1i * q / (a * d - b * c) * a * ds) * exp(1i * q / (a * d - b * c) * b * dr) * exp(1i * q / (a * d - b * c) * a) - k__2 * lambda ^ (0.1e1 / (a * d - b * c) * c * ds) / lambda ^ (0.1e1 / (a * d - b * c) * d * dr) * lambda ^ (0.1e1 / (a * d - b * c) * c) / exp(1i * q / (a * d - b * c) * a * ds) * exp(1i * q / (a * d - b * c) * b * dr) / exp(1i * q / (a * d - b * c) * a); -k__1 / lambda ^ (0.1e1 / (a * d - b * c) * c * ds) * lambda ^ (0.1e1 / (a * d - b * c) * d * dr) * lambda ^ (0.1e1 / (a * d - b * c) * c) * exp(1i * q / (a * d - b * c) * a * ds) / exp(1i * q / (a * d - b * c) * b * dr) / exp(1i * q / (a * d - b * c) * a) - k__2 / lambda ^ (0.1e1 / (a * d - b * c) * c * ds) * lambda ^ (0.1e1 / (a * d - b * c) * d * dr) / lambda ^ (0.1e1 / (a * d - b * c) * c) * exp(1i * q / (a * d - b * c) * a * ds) / exp(1i * q / (a * d - b * c) * b * dr) * exp(1i * q / (a * d - b * c) * a) (2 * k__0) - k__0 / lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) - k__0 * lambda ^ (0.1e1 / (a * d - b * c) * d) / exp(1i * q / (a * d - b * c) * b) + k__1 + k__2;];
lambda = solve(det(subs(A)) == 0, lambda);
q = linspace(-5*pi, 5*pi, 101);
lambda = double(subs(lambda));

% A = [(2 * k__0) - k__0 / lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) - k__0 * lambda ^ (0.1e1 / (a * d - b * c) * d) / exp(1i * q / (a * d - b * c) * b) + k__1 + k__2 -k__1 * lambda ^ (0.1e1 / (a * d - b * c) * c * ds) / lambda ^ (0.1e1 / (a * d - b * c) * d * dr) / lambda ^ (0.1e1 / (a * d - b * c) * c) / exp(1i * q / (a * d - b * c) * a * ds) * exp(1i * q / (a * d - b * c) * b * dr) * exp(1i * q / (a * d - b * c) * a) - k__2 * lambda ^ (0.1e1 / (a * d - b * c) * c * ds) / lambda ^ (0.1e1 / (a * d - b * c) * d * dr) * lambda ^ (0.1e1 / (a * d - b * c) * c) / exp(1i * q / (a * d - b * c) * a * ds) * exp(1i * q / (a * d - b * c) * b * dr) / exp(1i * q / (a * d - b * c) * a); -k__1 / lambda ^ (0.1e1 / (a * d - b * c) * c * ds) * lambda ^ (0.1e1 / (a * d - b * c) * d * dr) * lambda ^ (0.1e1 / (a * d - b * c) * c) * exp(1i * q / (a * d - b * c) * a * ds) / exp(1i * q / (a * d - b * c) * b * dr) / exp(1i * q / (a * d - b * c) * a) - k__2 / lambda ^ (0.1e1 / (a * d - b * c) * c * ds) * lambda ^ (0.1e1 / (a * d - b * c) * d * dr) / lambda ^ (0.1e1 / (a * d - b * c) * c) * exp(1i * q / (a * d - b * c) * a * ds) / exp(1i * q / (a * d - b * c) * b * dr) * exp(1i * q / (a * d - b * c) * a) (2 * k__0) - k__0 / lambda ^ (0.1e1 / (a * d - b * c) * d) * exp(1i * q / (a * d - b * c) * b) - k__0 * lambda ^ (0.1e1 / (a * d - b * c) * d) / exp(1i * q / (a * d - b * c) * b) + k__1 + k__2;];
% for i = 0:100
%     q = -pi + i*2*pi/100;
%     q_result(i+1) = -pi + i*2*pi/100;
%     lambda_result(i+1)=solve(det(subs(A)) == 0, lambda, ReturnConditions=true);
% end
% q = q_result;
% lambda = lambda_result;
pm = sprintf('k_0 = %d, k_1 = %d, k_2 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d), (dr, ds) = (%d, %d)', k__0, k__1, k__2, a, b, c, d, dr, ds);

figure%('Visible','off')
plot(q, abs(lambda))
title("|\lambda|-q", pm)
xlabel("q")
ylabel("|\lambda|")
axis([-pi pi 0 1])
if flag == 1
    saveas(gcf, strcat('PBC ', pm, ' 1.png'));
end

figure%('Visible','off')
for i = 1:max(size(lambda(:,1)))
    plot(q, angle(lambda(i,:)))
    hold on
end
title("arg(\lambda)-q", pm)
xlabel("q")
ylabel("arg(\lambda)")
if flag == 1
    saveas(gcf, strcat('PBC ', pm, ' 2.png'));
end

figure%('Visible','off')
for i = 1:max(size(lambda(:,1)))
    plot(real(lambda(i,:)), imag(lambda(i,:)))
    hold on
end
title("Im(\lambda)-Re(\lambda)", pm)
xlabel("Re(\lambda)")
ylabel("Im(\lambda)")
axis([-1 1 -1 1])
if flag == 1
    saveas(gcf, strcat('PBC ', pm, ' 3.png'));
end
