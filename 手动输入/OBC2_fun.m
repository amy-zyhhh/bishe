function [OBC2_k_2, OBC2_lambda] = OBC2_fun(a, b, c, d, dr, ds, k__0, k__1, N, k__2max, k__2min, p, flag)
    
    tic;
    
    syms lambda k__b k__2 q
    A = sym('temp', [N, N]);
    A(:,:) = 0;
    %%%%%%%%%%%%
    
    e_1 = [a; b];
    e_2 = [c; d];
    
    syms n m
    p1 = solve(n*e_1 + m*e_2 == [1; 0], {n, m});
    p2 = solve(n*e_1 + m*e_2 == [-1; 0], {n, m});
    p3 = solve(n*e_1 + m*e_2 == [0 - dr; 1 - ds], {n, m});
    p4 = solve(n*e_1 + m*e_2 == [0 - dr; -1 - ds], {n, m});
    p5 = solve(n*e_1 + m*e_2 == [0 + dr; 1 + ds], {n, m});
    p6 = solve(n*e_1 + m*e_2 == [0 + dr; -1 + ds], {n, m});
    
    for i = 1:N
        if i ~= 1 && i ~= N
            A(N + 1 - i, N + 1 - i) = 2*k__0 + k__1 + k__2;
        else
            A(N + 1 - i, N + 1 - i) = k__b;
        end
        if mod(i, 2) == 1
            if N + 1 - i - 2*p1.m <= N && N + 1 - i - 2*p1.m > 0
                A(N + 1 - i, N + 1 - i - 2*p1.m) = A(N + 1 - i, N + 1 - i - 2*p1.m) - k__0*lambda^(p1.n)*exp(1i*q*p1.m);
            end
            if N + 1 - i - 2*p2.m <= N && N + 1 - i - 2*p2.m > 0
                A(N + 1 - i, N + 1 - i - 2*p2.m) = A(N + 1 - i, N + 1 - i - 2*p2.m) - k__0*lambda^(p2.n)*exp(1i*q*p2.m);
            end
            if N + 1 - i - (2*p3.m + 1) <= N && N + 1 - i - (2*p3.m + 1) > 0
                A(N + 1 - i, N + 1 - i - (2*p3.m + 1)) = A(N + 1 - i, N + 1 - i - (2*p3.m + 1)) - k__2*lambda^(p3.n)*exp(1i*q*p3.m);
            end
            if N + 1 - i - (2*p4.m + 1) <= N && N + 1 - i - (2*p4.m + 1) > 0
                A(N + 1 - i, N + 1 - i - (2*p4.m + 1)) = A(N + 1 - i, N + 1 - i - (2*p4.m + 1)) - k__1*lambda^(p4.n)*exp(1i*q*p4.m);
            end
        else
            if N + 1 - i - 2*p1.m <= N && N + 1 - i - 2*p1.m > 0
                A(N + 1 - i, N + 1 - i - 2*p1.m) = A(N + 1 - i, N + 1 - i - 2*p1.m) - k__0*lambda^(p1.n)*exp(1i*q*p1.m);
            end
            if N + 1 - i - 2*p2.m <= N && N + 1 - i - 2*p2.m > 0
                A(N + 1 - i, N + 1 - i - 2*p2.m) = A(N + 1 - i, N + 1 - i - 2*p2.m) - k__0*lambda^(p2.n)*exp(1i*q*p2.m);
            end
            if N + 1 - i - (2*p5.m - 1) <= N && N + 1 - i - (2*p5.m - 1) > 0
                A(N + 1 - i, N + 1 - i - (2*p5.m - 1)) = A(N + 1 - i, N + 1 - i - (2*p5.m - 1)) - k__1*lambda^(p5.n)*exp(1i*q*p5.m);
            end
            if N + 1 - i - (2*p6.m - 1) <= N && N + 1 - i - (2*p6.m - 1) > 0
                A(N + 1 - i, N + 1 - i - (2*p6.m - 1)) = A(N + 1 - i, N + 1 - i - (2*p6.m - 1)) - k__2*lambda^(p6.n)*exp(1i*q*p6.m);
            end
        end
    end
    
    for i = 0:p
        k__2 = k__2min + (k__2max - k__2min)*i/p;
        k__b = 2*k__0 + k__1 +k__2;
        OBC2_lambda(:,i + 1) = double(solve(det(subs(A)) == 0, lambda));
        OBC2_k_2(i + 1) = k__2;
        i
        toc
    end
    
    pm = sprintf('k_0 = %d, k_1 = %d, e_1(n) = (%d, %d), e_2(m) = (%d, %d), (dr, ds) = (%d, %d)', k__0, k__1, a, b, c, d, dr, ds);
    figure
    for i = 1:max(size(OBC2_lambda(:,1)))
        if isempty(OBC2_lambda) == 1
            ;
        else
            plot(OBC2_k_2, OBC2_lambda(i,:))
        end
        hold on
    end
    title("|\lambda|-k_2", pm)
    xlabel("k_2")
    ylabel("|\lambda|")
    if flag == 1
        saveas(gcf, strcat('OBC2 ', pm, '.png'));
    end

end