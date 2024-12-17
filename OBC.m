clear; close all; clc;

tic
syms a b c d r s dr ds
syms k__0 k__1 k__2 k__b
syms lambda q

a = 1;
b = 0;
c = 0;
d = 1;

%%%%%%%%%%%%%%%%%%%
v = 1;
flag = 1;
for a = 1:3
    for b = 0:3
        for c = 0:3
            for d = 1:3
                if abs(a*d-b*c) == 1
%%%%%%%%%%%%%%%%%%%

N = 10;

k__0 = 1;
k__1 = 6;
k__2min = 1;
k__2max = 11;
p = 10;

e_1 = [a; b];
e_2 = [c; d];
syms n m
p1 = solve(n*e_1 + m*e_2 == [1; 0], {n, m});
p2 = solve(n*e_1 + m*e_2 == [-1; 0], {n, m});
p3 = solve(n*e_1 + m*e_2 == [0; 1], {n, m});
p4 = solve(n*e_1 + m*e_2 == [0; -1], {n, m});

for i = 1:N
    if i ~= 1 && i ~= N
        A(N + 1 - i, N + 1 - i) = 2*k__0 + k__1 + k__2;
    else
        A(N + 1 - i, N + 1 - i) = k__b;
    end
    if mod(i, 2) == 1
        if N + 1- i - p1.m <= N && N + 1- i - p1.m > 0
    A(N + 1 - i, N + 1- i - p1.m) = A(N + 1 - i, N + 1- i - p1.m) - k__0*lambda^(p1.n);
        end
        if N + 1- i - p2.m <= N && N + 1- i - p2.m > 0
    A(N + 1 - i, N + 1- i - p2.m) = A(N + 1 - i, N + 1- i - p2.m) - k__0*lambda^(p2.n);
        end
        if N + 1- i - p3.m <= N && N + 1- i - p3.m > 0
    A(N + 1 - i, N + 1- i - p3.m) = A(N + 1 - i, N + 1- i - p3.m) - k__2*lambda^(p3.n);
        end
        if N + 1- i - p4.m <= N && N + 1- i - p4.m > 0
    A(N + 1 - i, N + 1- i - p4.m) = A(N + 1 - i, N + 1- i - p4.m) - k__1*lambda^(p4.n);
        end
    else
        if N + 1- i - p1.m <= N && N + 1- i - p1.m > 0
    A(N + 1 - i, N + 1- i - p1.m) = A(N + 1 - i, N + 1- i - p1.m) - k__0*lambda^(p1.n);
        end
        if N + 1- i - p2.m <= N && N + 1- i - p2.m > 0
    A(N + 1 - i, N + 1- i - p2.m) = A(N + 1 - i, N + 1- i - p2.m) - k__0*lambda^(p2.n);
        end
        if N + 1- i - p3.m <= N && N + 1- i - p3.m > 0
    A(N + 1 - i, N + 1- i - p3.m) = A(N + 1 - i, N + 1- i - p3.m) - k__1*lambda^(p3.n);
        end
        if N + 1- i - p4.m <= N && N + 1- i - p4.m > 0
    A(N + 1 - i, N + 1- i - p4.m) = A(N + 1 - i, N + 1- i - p4.m) - k__2*lambda^(p4.n);
        end
    end
end

for i = 0:p
    k__2 = k__2min + (k__2max - k__2min)*i/p;
    k__b = 2*k__0 + k__1 +k__2;
    OBC_lambda(:,i + 1) = double(solve(det(subs(A)) == 0, lambda));
    OBC_k_2(i + 1) = k__2;
    i
    toc
end

pm = sprintf('k_0 = %d, k_1 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d)', k__0, k__1, a, b, c, d);
figure(Visible="off")
for i = 1:max(size(OBC_lambda(:,1)))
    plot(OBC_k_2, OBC_lambda(i,:))
    hold on
end
title("|\lambda|-k_2", pm)
xlabel("k_2")
ylabel("|\lambda|")
if flag == 1
    saveas(gcf, strcat('OBC_', pm, '.png'));
end


%%%%%%%%%%%%%%%
OBC_result(v).k_2 = OBC_k_2;
OBC_result(v).lambda = OBC_lambda;
v = v + 1;
OBC_k_2 = [];
OBC_lambda = [];
                end
            end
        end
    end
end