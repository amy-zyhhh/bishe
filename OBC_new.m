clear; close all; clc;
tic

N = 20;
a = 1;
b = 0;
c = -1;
d = 1;

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

%%%%%%%%%%%%%%%%%
i = 1;
k__2 = OBC_k__2(i);
for j = 1:max(size(OBC_lambda(:,i)))
    lambda = OBC_lambda(j,i);
    B = double(subs(A));
    OBC_d(:,j) = null(B, 1e-5);
end

j = 10;
lambda = OBC_lambda(j,i);
figure
pm2 = sprintf('k_0 = %d, k_1 = %d, lambda = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, lambda, a, b, c, d);
scatter(1:N, OBC_d(:,i))
title("d-m",pm2)
xlabel("m")
ylabel("\lambda")

[X,Y] = meshgrid(1:N,vpa(OBC_lambda(1:N,i)));
Y = double(Y);

% OBC_d归一化，分离对称和反对称
u = 1;
v = 1;
for i = 1:2*N
    if OBC_d(1, i)*OBC_d(N, i) > 0
        OBC_d2(u,:) = OBC_d(:, i)./max(OBC_d(:, i));
        u = u + 1;
    else
        OBC_d1(v,:) = OBC_d(:, i)./max(OBC_d(:, i));
        v = v + 1;
    end
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