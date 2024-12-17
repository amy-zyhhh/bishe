clear;close all;clc;

% eq25,eq27绘制衰减谱
syms theta k__1 k__2 k__0
q = linspace(-pi,pi,100);
% 25[1]
lambda__1 = (-sqrt(-0.4e1 .* (cos(theta) .* (k__1 + k__2) / 0.2e1 + k__0) .* abs(cos(theta)) .* sqrt(k__2 .^ 2 + k__1 .^ 2 + 0.2e1 .* k__1 .* k__2 .* cos(q)) + 0.4e1 .* cos(theta) .* ((k__1 .* k__2 .* cos(q) + k__1 .^ 2 + k__1 .* k__2 + k__2 .^ 2) .* cos(theta) / 0.2e1 + k__0 .* (k__1 + k__2))) - sqrt((exp(q .* i) .* k__2 + k__1) .* (exp(-i .* q) .* k__2 + k__1)) .* abs(cos(theta)) + (cos(theta) .* (k__1 + k__2)) + (0.2e1 .* k__0)) / k__0 / 2;
% 22[1]
lambda__2 = (-sqrt(0.4e1 .* (cos(theta) .* (k__1 + k__2) / 0.2e1 + k__0) .* abs(cos(theta)) .* sqrt(k__2 .^ 2 + k__1 .^ 2 + 0.2e1 .* k__1 .* k__2 .* cos(q)) + 0.4e1 .* cos(theta) .* ((k__1 .* k__2 .* cos(q) + k__1 .^ 2 + k__1 .* k__2 + k__2 .^ 2) .* cos(theta) / 0.2e1 + k__0 .* (k__1 + k__2))) + sqrt((exp(q .* i) .* k__2 + k__1) .* (exp(-i .* q) .* k__2 + k__1)) .* abs(cos(theta)) + (cos(theta) .* (k__1 + k__2)) + (0.2e1 .* k__0)) / k__0 / 2;


for k = 14:20
    k__0 = 1;
    k__1 = 2;
    k__2 = 2;
    theta = k * pi/20;
    plot__1 = double(subs(lambda__1));
    plot__2 = double(subs(lambda__2));
    % figure(k)
    % plot(q, real(plot__1), q, real(plot__2),q, imag(plot__1), q, imag(plot__2))
    plot(q, plot__1, q, plot__2)
    title('衰减谱随theta的变化关系')
    title('三角和指数位移下的衰减谱')
    hold on
    xlabel('q')
    ylabel('lambda')
    hold on

    % k__0 = 1; k__1 = 2; k__2 = 2; theta = pi/20;
    % 23[1] cos
    plot__t1 = (-sqrt(0.4e1 .* (cos(theta) .* (k__1 + k__2) / 0.2e1 + k__0) .* (cos(q) .* k__2 + k__1) .* abs(cos(theta)) .* sign(cos(q) .* k__2 + k__1) + 0.4e1 .* ((cos(q) .^ 2 .* k__2 .^ 2 / 0.2e1 + cos(q) .* k__1 .* k__2 + k__1 .^ 2 + k__1 .* k__2 + k__2 .^ 2 / 0.2e1) .* cos(theta) / 0.2e1 + k__0 .* (k__1 + k__2)) .* cos(theta)) + sign(cos(q) .* k__2 + k__1) .* (cos(q) .* k__2 + k__1) .* abs(cos(theta)) + cos(theta) .* (k__1 + k__2) + 0.2e1 .* k__0) / k__0 / 0.2e1;
    % 27[1] cos
    plot__t2 = (-sqrt(-0.4e1 .* (cos(theta) .* (k__1 + k__2) / 0.2e1 + k__0) .* (cos(q) .* k__2 + k__1) .* abs(cos(theta)) .* sign(cos(q) .* k__2 + k__1) + 0.4e1 .* ((cos(q) .^ 2 .* k__2 .^ 2 / 0.2e1 + cos(q) .* k__1 .* k__2 + k__1 .^ 2 + k__1 .* k__2 + k__2 .^ 2 / 0.2e1) .* cos(theta) / 0.2e1 + k__0 .* (k__1 + k__2)) .* cos(theta)) + (-cos(q) .* k__2 - k__1) .* abs(cos(theta)) .* sign(cos(q) .* k__2 + k__1) + cos(theta) .* (k__1 + k__2) + 0.2e1 .* k__0) / k__0 / 0.2e1;
    % 30[1] sin
    plot__t3 = (-sqrt(-0.4e1 .* abs(cos(theta)) .* (cos(theta) .* (k__1 + k__2) / 0.2e1 + k__0) .* sqrt(k__1 .^ 2 - sin(q) .^ 2 .* k__2 .^ 2) + 0.4e1 .* ((cos(q) .^ 2 .* k__2 .^ 2 / 0.2e1 + k__1 .* (k__1 + k__2)) .* cos(theta) / 0.2e1 + k__0 .* (k__1 + k__2)) .* cos(theta)) - sqrt(k__1 .^ 2 - sin(q) .^ 2 .* k__2 .^ 2) .* abs(cos(theta)) + cos(theta) .* (k__1 + k__2) + 0.2e1 .* k__0) / k__0 / 0.2e1;
    % 33[1] sin
    plot__t4 = (-sqrt(0.4e1 .* abs(cos(theta)) .* (cos(theta) .* (k__1 + k__2) / 0.2e1 + k__0) .* sqrt(k__1 .^ 2 - sin(q) .^ 2 .* k__2 .^ 2) + 0.4e1 .* ((cos(q) .^ 2 .* k__2 .^ 2 / 0.2e1 + k__1 .* (k__1 + k__2)) .* cos(theta) / 0.2e1 + k__0 .* (k__1 + k__2)) .* cos(theta)) + sqrt(k__1 .^ 2 - sin(q) .^ 2 .* k__2 .^ 2) .* abs(cos(theta)) + cos(theta) .* (k__1 + k__2) + 0.2e1 .* k__0) / k__0 / 0.2e1;

    plot__t1 = double(subs(plot__t1));
    plot__t2 = double(subs(plot__t2));
    plot__t3 = double(subs(plot__t3));
    plot__t4 = double(subs(plot__t4));

    plot(q, plot__t1)
    hold on
    plot(q, plot__t2)
    hold on
    plot(q, plot__t3)
    hold on
    plot(q, plot__t4)
    legend('e1', 'e2', 'c1', 'c2', 's1', 's2')

end


% clear;clc;
% theta = linspace(0,pi,1000);
% y = -0.2e1 .* sqrt(0.2e1) .* cos(theta / 0.2e1) .* sqrt(cos(theta)) + 0.1e1 + 0.2e1 .* cos(theta);
% plot(theta, real(y), theta, imag(y));
% title('能带闭合位置')
% xlabel('lambda')
% ylabel('theta')
% legend('实部', '虚部')