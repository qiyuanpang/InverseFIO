function [YA, YB, YC, T, R, rk] = hodlrqr(A, Bl, Br, C, lvls, norm2, tol_or_rank)
    if ~isempty(Bl)
        sizebl = size(Bl);
        if norm(Bl'*Bl - eye(sizebl(2))) >= 1e-13
            [Bl, R] = qr(Bl, 0);
            Br = R*Br;
        end
    end
    if lvls == 0
        H = [A.A11;Br;C];
        m = A.m;
        n = A.n;
        [Y, T, R] = qr_wybased(H);
        YA = Y(1:m,:);
        YA = struct('A11', YA, 'A12', [], 'A21', [], 'A22', [], 'lvl', lvls);
        YB = Y(m+1:m+size(Br,1),:);
        YB = struct('Q', Bl, 'R', YB);
        YC = Y(m+size(Br,1)+1:end,:);
        R = struct('A11', R, 'A12', [], 'A21', [], 'A22', [], 'lvl', lvls);
        T = struct('A11', T, 'A12', [], 'A21', [], 'A22', [], 'lvl', lvls);
        rk = 1;
    else
        % rk1 = A.A21.r;
        n1 = A.A21.n;
        n2 = A.A12.n;
        [YA1, YB1, YC1, T1, R1, rk] = hodlrqr(A.A11, A.A21.Q, A.A21.R, [Br(:, 1:min(n1, size(Br,2))); C(:, 1:min(n1, size(C,2)))], lvls-1, norm2, tol_or_rank);
        S = struct('Q', [], 'R', []);
        S.Q = [hodlr_adj_apply(YA1, A.A12.Q) YB1.R' YC1(1:size(Br,1), :)' YC1(size(Br,1)+1:size(Br,1)+size(C,1), :)'];
        S.R = [A.A12.R; hodlr_adj_apply(A.A22, YB1.Q)'; Br(:, n1+1:min(n1+n2,size(Br,2))); C(:, n1+1:min(n1+n2,size(C,2)))];
        [QS, RS] = recompression(S.Q, S.R, tol_or_rank*norm2);
        rk = max(rk, size(QS,2));
        QS = hodlr_adj_apply(T1, QS);
        [Q12, R12] = recompression([A.A12.Q -hodlr_apply(YA1, QS)], [A.A12.R; RS], tol_or_rank*norm2);
        rk = max(rk, size(Q12,2));
        % A22 = hodlr(hodlr_apply(A.A22, eye(A.A22.n))-YB1.Q*YB1.R*QS*RS, lvls-1, tol_or_rank*norm2);
        A22 = hodlr_sub_lr(A.A22, YB1.Q, YB1.R*QS*RS, tol_or_rank*norm2);
        rk = max(rk, A22.maxrk);
        if isempty(Br)
            Br2 = [];
        else
            Br2 = Br(:, n1+1:min(n1+n2,size(Br,2))) - (YC1(1:size(Br,1), :)*QS)*RS;
        end
        if isempty(C)
            C2 = [];
        else
            C2 = C(:, n1+1:min(n1+n2,size(C,2))) - (YC1(size(Br,1)+1:size(Br,1)+size(C,1), :)*QS)*RS;
        end
        [YA2, YB2, YC2, T2, R2, rk2] = hodlrqr(A22, [], [], [Br2;C2], lvls-1, norm2, tol_or_rank);
        T12 = struct('Q', [], 'R', []);
        T12.Q = [YB1.R' YC1(1:size(Br,1), :)' YC1(size(Br,1)+1:size(Br,1)+size(C,1), :)'];
        T12.R = [hodlr_adj_apply(YA2, YB1.Q)'; YC2(1:size(Br2,1), :); YC2(size(Br2,1)+1:size(Br2,1)+size(C2,1), :)];
        [QT12, RT12] = recompression(T12.Q, T12.R, tol_or_rank);
        rk2 = max(rk2, size(QT12, 2));
        T = struct('A11', T1, 'A12', [], 'A21', [], 'A22', T2, 'lvl', lvls);
        T.A12 = struct('Q', -hodlr_apply(T1, QT12), 'R', hodlr_adj_apply(T2, RT12')');
        R = struct('A11', R1, 'A12', [], 'A21', [], 'A22', R2, 'lvl', lvls);
        R.A12 = struct('Q', Q12, 'R', R12);
        YA = struct('A11', YA1, 'A12', [], 'A21', YB1, 'A22', YA2, 'lvl', lvls);
        YB = [YC1(1:size(Br,1), :) YC2(1:size(Br2,1), :)];
        YB = struct('Q', Bl, 'R', YB);
        YC = [YC1(size(Br,1)+1:size(Br,1)+size(C,1), :) YC2(size(Br2,1)+1:size(Br2,1)+size(C2,1), :)];
        rk = max(rk, rk2);
    end
end