function [OBC1] = OBC1_fun(c, d, dr, ds, k__0, k__1, N, k__2max, k__2min, p, flag)

    tic;
    v = 1; % 对不同e_1的计数
    [t_temp, e_1_temp] = find_e_1(c, d, 1);
    if t_temp == 0 % 防止e_2竖直向上
        e_1_temp(1,1:2) = [1, 0];
        t_temp = 2;
    end

    for u = 1:t_temp-1
        a = e_1_temp(u, 1);
        b = e_1_temp(u, 2);
        syms lambda k__b k__2
        A = sym('temp', [N, N]);
        A(:,:) = 0;
        
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
        
        pm = sprintf('k_0 = %d, k_1 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d), [(dr, ds) = (%d, %d)]', k__0, k__1, a, b, c, d, dr, ds);
        figure('Visible','off')
        for i = 1:max(size(OBC_lambda(:,1)))
            if isempty(OBC_lambda) == 1
                ;
            else
                plot(OBC_k_2, OBC_lambda(i,:))
                hold on
            end
        end
        title("|\lambda|-k_2", pm)
        xlabel("k_2")
        ylabel("|\lambda|")
        if flag == 1
            saveas(gcf, strcat('OBC1 ', pm, '.png'));
        end
        OBC1(v).k_2 = OBC_k_2;
        OBC1(v).lambda = OBC_lambda;
        v = v + 1;
        OBC_k_2 = [];
        OBC_lambda = [];
        
    end
end
