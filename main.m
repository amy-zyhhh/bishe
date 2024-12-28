% clear; close all; clc;

a = 1;
b = 0;
c = 2;
d = 2;
dr = 0;
ds = 1;
k__0 = 1;
k__1 = 6;
k__2 = 6;
N = 12;
p = 10; % OBC中计算点数
k__2max = 11;
k__2min = 1;
flag = 0; % 1保存，0不保存

% for c = 0:5
%     for d = 2:2:10
%         flag1 = 0; % 用于判断向量是否穿过1点，0表示不穿过
%         for i = 0:c
%             for j = 0:2:d
%                 if i*d == c*j && norm([i, j] - [0, 0]) ~= 0 && norm([i, j] - [c, d]) ~= 0
%                     flag1 = 1;
%                 else
%                     if (i == 0 && j == 0) || (i == c && j == d)
%                         ;
%                     else
%                         flag1 = 0;
%                     end
%                 end
%             end
%         end
%         if flag1 == 0
%             [t, e_1] = find_e_1(c, d, 0);
%             if t == 0 % 防止e_2竖直向上
%                 e_1(1,1:2) = [1, 0];
%                 t = 2;
%             end
%             for u = 1:t-1
%                 a = e_1(u, 1);
%                 b = e_1(u, 2);

%% PBC
[PBC_q, PBC_lambda] = PBC_fun(a, b, c, d, dr, ds, k__0, k__1, k__2, flag);

%% OBC1
if mod(c, 2) == 0
    d_temp = d/2;
    c_temp = c/2;
else
    c_temp = c;
    d_temp = d;
end
[OBC1] = OBC1_fun(c_temp, d_temp, dr, ds, k__0, k__1, N, k__2max, k__2min, p, flag);

%% OBC2
[OBC2_k_2, OBC2_lambda] = OBC2_fun(a, b, c, d, dr, ds, k__0, k__1, N, k__2max, k__2min, p, flag);

%             end
%         end
%     end
% end