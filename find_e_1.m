function [t, e_1] = find_e_1(c, d, flag3) % flag为1，表示1点和2点都考虑，否则只考虑1点
    h_min = sqrt(c^2 + (d)^2); % 防止取不到最小值
    t = 0; % 防止出现多个e_1
    e_1 = [0, 0];
    for i = 0:c
        for j = 0:d % 不能取d/2，因为(c,d)可能没有跨过2节点
            if (i == 0 && j == 0) || (i == c && j == d)
                ;
            else
                if flag3 == 0 && mod(j, 2) == 1
                    ;
                else
                    h = sqrt(j ^ 2 + i ^ 2) * sqrt(1 - (conj(c) * i + conj(d) * j) ^ 2 / (abs(c) ^ 2 + abs(d) ^ 2) / (abs(i) ^ 2 + abs(j) ^ 2));
                    if h_min - h > 10^(-5) % 直接比较有误差影响
                        h_min = h;
                        e_1(:,1:2) = 0;
                        t = 1;
                        e_1(t,1:2) = [i; j];
                        t = t + 1;
                    elseif abs(h_min - h) < 10^(-5)
                        e_1(t,1:2) = [i; j];
                        t = t + 1;
                    end
                end
            end
        end
    end

end