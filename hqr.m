function [Y T R rk] = hqr(A, Bl, Br, C, level, tol_or_rank)
    if ~isempty(Bl)
        sizebl = size(Bl);
        if norm(Bl'*Bl - eye(sizebl(2))) >= 1e-13
            [Bl, R] = qr(Bl, 0);
            Br = R*Br;
        end
    end
    H = [A;Br;C];
    if level == 0
        [Y, T, R] = qr_wybased(H);
        % disp(size(H))
        YA = Y(1:size(A,1),:);
        YB = Y(size(A,1)+1:size(Br,1)+size(A,1),:);
        YC = Y(size(A,1)+size(Br,1)+1:end,:);
        Y = [YA;Bl*YB;YC];
        rk = size(R, 1);
        % rk = 0;
    else
        [m, n] = size(A);
        m2 = floor(m/2);
        n2 = floor(n/2);
        [Bl1, Br1] = myqr(H(m2+1:m, 1:n2), tol_or_rank);
        rk = min(size(Br1));
        % fprintf('%d, %d \n', rk, level)
        [Y1, T1, R1, rk1] = hqr(H(1:m2, 1:n2), Bl1, Br1, [H(m+1:m+size(Br,1), 1:n2); H(m+size(Br,1)+1:m+size(Br,1)+size(C,1), 1:n2)], level-1, tol_or_rank);
        S = Y1(1:m2,:)'*H(1:m2,n2+1:n) + Y1(m2+1:m,:)'*H(m2+1:m,n2+1:n) + Y1(m+1:m+size(Br,1),:)'*H(m+1:m+size(Br,1),n2+1:n) + Y1(m+size(Br,1)+1:end,:)'*H(m+size(Br,1)+1:m+size(Br,1)+size(C,1),n2+1:n);
        S = T1'*S;
        H(1:m2,n2+1:n) = H(1:m2,n2+1:n) - Y1(1:m2,:)*S;
        H(m2+1:m,n2+1:n) = H(m2+1:m,n2+1:n) - Y1(m2+1:m,:)*S;
        H(m+1:m+size(Br,1),n2+1:n) = H(m+1:m+size(Br,1),n2+1:n) - Y1(m+1:m+size(Br,1),:)*S;
        H(m+size(Br,1)+1:end,n2+1:n) = H(m+size(Br,1)+1:end,n2+1:n) - Y1(m+size(Br,1)+1:end,:)*S;
        
        [Y2, T2, R2, rk2] = hqr(H(m2+1:m,n2+1:n), [], [], [H(m+1:m+size(Br,1), n2+1:n); H(m+size(Br,1)+1:m+size(Br,1)+size(C,1), n2+1:n)], level-1, tol_or_rank);
        T12 = Y1(m2+1:m,:)'*Y2(1:m-m2,:) + Y1(m+1:m+size(Br,1),:)'*Y2(m-m2+1:m-m2+size(Br,1),:) + Y1(m+size(Br,1)+1:end,:)'*Y2(m-m2+size(Br,1)+1:end,:);
        T12 = -T1*T12*T2;
        T = [T1 T12;zeros(size(T2,1), size(T1,2)) T2];
        R = [R1 H(1:m2,n2+1:n); zeros(size(R2,1), size(R1,2)) R2];
        YA = [Y1(1:m2,:) zeros(m2, size(Y2,2)); Y1(m2+1:m,:) Y2(1:m-m2,:)];
        YB = [Y1(m+1:m+size(Br,1),:) Y2(m-m2+1:m-m2+size(Br,1),:)];
        YC = [Y1(m+size(Br,1)+1:end,:) Y2(m-m2+size(Br,1)+1:end,:)];
        Y = [YA;Bl*YB;YC];
        rk = max(rk, max(rk1, rk2));
    end
end