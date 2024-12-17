clear; close all; clc;

a = 1;
b = 0;
c = 2;
d = 2;
dr = 0;
ds = 1;

N = 40;
p = 10; % OBC中计算点数
k__2max = 11;
k__2min = 1;
flag = 0; % 1保存，0不保存

%% PBC
k__0 = 1;
k__1 = 6;
k__2 = 3;
[PBC_q1, PBC_lambda1] = PBC_fun(a, b, c, d, dr, ds, k__0, k__1, k__2, flag);
k__0 = 1;
k__1 = 6;
k__2 = 6;
[PBC_q2, PBC_lambda2] = PBC_fun(a, b, c, d, dr, ds, k__0, k__1, k__2, flag);
k__0 = 1;
k__1 = 6;
k__2 = 9;
[PBC_q3, PBC_lambda3] = PBC_fun(a, b, c, d, dr, ds, k__0, k__1, k__2, flag);

%% OBC1
k__0 = 1;
k__1 = 6;
if mod(c, 2) == 0
    d_temp = d/2;
    c_temp = c/2;
else
    c_temp = c;
    d_temp = d;
end
[OBC1] = OBC1_fun(c_temp, d_temp, dr, ds, k__0, k__1, N, k__2max, k__2min, p, flag);
