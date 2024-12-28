clear; close all; clc;
tic

N = 20;
a = 1;
b = 2;
c = 1;
d = 3;

k__0 = 1;
k__1 = 6;
syms k__2
p = 11;
k__2min = 1;
k__2max = 11;
k__b = 2*k__0 + k__1 + k__2;

syms lambda

e_1 = [a; b];
e_2 = [c; d];
syms n m
p1 = solve(n*e_1 + m*e_2 == [1; 0], {n, m});
p2 = solve(n*e_1 + m*e_2 == [-1; 0], {n, m});
p3 = solve(n*e_1 + m*e_2 == [0; 1], {n, m});
p4 = solve(n*e_1 + m*e_2 == [0; -1], {n, m});

A = sym('a', N);
A(:,:) = 0;

for i = 1:N
    if i == 1 || i == N
        A(N + 1 - i, N + 1 - i) = k__b;
    else
        A(N + 1 - i, N + 1 - i) = 2*k__0 + k__1 + k__2;
    end
    if mod(i,2) == 1
        if N + 1 - i - p1.m > 0 && N + 1 - i - p1.m <=N
        A(N + 1 - i, N + 1 - i - p1.m) = A(N + 1 - i, N + 1 - i - p1.m) - k__0*lambda^(p1.n);
        end
        if N + 1 - i - p2.m > 0 && N + 1 - i - p2.m <=N
        A(N + 1 - i, N + 1 - i - p2.m) = A(N + 1 - i, N + 1 - i - p2.m) - k__0*lambda^(p2.n);
        end
        if N + 1 - i - p3.m > 0 && N + 1 - i - p3.m <=N
        A(N + 1 - i, N + 1 - i - p3.m) = A(N + 1 - i, N + 1 - i - p3.m) - k__2*lambda^(p3.n);
        end
        if N + 1 - i - p4.m > 0 && N + 1 - i - p4.m <=N
        A(N + 1 - i, N + 1 - i - p4.m) = A(N + 1 - i, N + 1 - i - p4.m) - k__1*lambda^(p4.n);
        end
    else
        if N + 1 - i - p1.m > 0 && N + 1 - i - p1.m <=N
        A(N + 1 - i, N + 1 - i - p1.m) = A(N + 1 - i, N + 1 - i - p1.m) - k__0*lambda^(p1.n);
        end
        if N + 1 - i - p2.m > 0 && N + 1 - i - p2.m <=N
        A(N + 1 - i, N + 1 - i - p2.m) = A(N + 1 - i, N + 1 - i - p2.m) - k__0*lambda^(p2.n);
        end
        if N + 1 - i - p3.m > 0 && N + 1 - i - p3.m <=N
        A(N + 1 - i, N + 1 - i - p3.m) = A(N + 1 - i, N + 1 - i - p3.m) - k__1*lambda^(p3.n);
        end
        if N + 1 - i - p4.m > 0 && N + 1 - i - p4.m <=N
        A(N + 1 - i, N + 1 - i - p4.m) = A(N + 1 - i, N + 1 - i - p4.m) - k__2*lambda^(p4.n);
        end
    end
end

for i = 1:p
    k__2 = k__2min + (k__2max - k__2min)*(i - 1)/p;
    OBC_lambda(:, i) = sort(vpa(solve(det(subs(A)) == 0,lambda))); % 每个k__2对应的lambda存一列
    OBC_k__2(i) = k__2;
end

figure
pm1 = sprintf('k_0 = %d, k_1 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, a, b, c, d);
for i = 1:max(size(OBC_lambda(:,1)))
    plot(OBC_k__2, OBC_lambda(i,:));
    hold on
end
title("|\lambda|-k_2", pm1)
xlabel("k_2")
ylabel("|\lambda|")

OBC_lambda2 = zeros(1,N);
for i = 1:p
    k__2 = OBC_k__2(i);
    num = 1;
    for j = 1:max(size(OBC_lambda(:,i)))
        lambda = OBC_lambda(j,i);
        B = double(subs(A));
        temp = null(B, 1e-1); % 对误差好像很敏感啊
        lambda_num(j) = max(size(temp(1,:)));
        % OBC_d的列表示特征向量，页表示刚度
        OBC_d(:,num:num + lambda_num(j) - 1, i) = temp;
        num = num + lambda_num(j);
        temp = max(size(OBC_lambda2(:,i))); % 这里temp变了
        if temp == 1 % 为了使初始值不影响赋值
            temp = temp - 1;
        end
        for k = 1:lambda_num(j) % 重根计重数
            OBC_lambda2(k + temp,i) = lambda;
        end
    end
end

i = 16; % 特征向量的个数
j = 6; % k__2的取值
lambda = OBC_lambda(i, j);
figure
pm2 = sprintf('k_0 = %d, k_1 = %d, lambda = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, lambda, a, b, c, d);
scatter(1:N, OBC_d(:,i,j))
title("d-m",pm2)
xlabel("m")
ylabel("\lambda")

i = 1; % 这里是k__2位置处的云图
[X,Y] = meshgrid(1:N,vpa(OBC_lambda2(1:N,i)));
Y = double(Y);
% OBC_d归一化
for j = 1:N
        OBC_d2(j,:) = OBC_d(:, j, i)./max(OBC_d(:, j, i));
end

figure
pm3 = sprintf('k_0 = %d, k_1 = %d, k_2 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, k__2, a, b, c, d);
pcolor(X,Y,OBC_d(:,1:N)');
% shading interp
colorbar
title("d-\lambda-m", pm3)
xlabel("m")
ylabel("\lambda")

figure
pm4 = sprintf('k_0 = %d, k_1 = %d, k_2 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, k__2, a, b, c, d);
pcolor(X,Y,OBC_d(:,N + 1:2*N)');
% shading interp
colorbar
title("d-\lambda-m",pm4)
xlabel("m")
ylabel("\lambda")

figure
pm3 = sprintf('k_0 = %d, k_1 = %d, k_2 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, k__2, a, b, c, d);
pcolor(X,Y,abs(OBC_d(:,1:N)'));
% shading interp
colorbar
title("d-|\lambda|-m", pm3)
xlabel("m")
ylabel("\lambda")

figure
pm3 = sprintf('k_0 = %d, k_1 = %d, k_2 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, k__2, a, b, c, d);
pcolor(X,Y,abs(OBC_d(:,1 + N:2*N)'));
% shading interp
colorbar
title("d-|\lambda|-m", pm3)
xlabel("m")
ylabel("\lambda")