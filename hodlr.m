function F = hodlr(A, lvl, tol_or_rank)
    if lvl == 0
        F = struct('A11', A, 'A12', [], 'A21', [], 'A22', [], 'lvl', lvl, 'm', size(A, 1), 'n', size(A,2), 'maxrk', 0);
    else if lvl > 0
        [m, n] = size(A);
        F = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', lvl, 'm', size(A, 1), 'n', size(A,2), 'maxrk', 0);
        m2 = floor(m/2);
        n2 = floor(n/2);
        [U1 S1 V1] = svd(A(1:m2, n2+1:n), 0);
        diagS1 = diag(S1);
        if tol_or_rank < 1
            tol = tol_or_rank*abs(S1(1));
            k = nnz(abs(diagS1) > tol);
            [Q1, R1] = recompression(U1(:, 1:k), diagS1(1:k).*V1(:, 1:k)', tol);
        else
            k = floor(tol_or_rank);
        end
        F.A12 = struct('Q', Q1, 'R', R1, 'lvl', lvl, 'm', m2, 'n', n-n2, 'r', size(Q1, 2));
        [U2 S2 V2] = svd(A(m2+1:m, 1:n2), 0);
        diagS2 = diag(S2);
        if tol_or_rank < 1
            tol = tol_or_rank*abs(S2(1));
            k = nnz(abs(diagS2) > tol);
            [Q2, R2] = recompression(U2(:, 1:k), diagS2(1:k).*V2(:, 1:k)', tol);
        else
            k = floor(tol_or_rank);
        end
        F.A21 = struct('Q', Q2, 'R', R2,  'lvl', lvl, 'm', m-m2, 'n', n2, 'r', size(Q2, 2));
        F.A11 = hodlr(A(1:m2, 1:n2), lvl-1, tol_or_rank);
        F.A22 = hodlr(A(m2+1:m, n2+1:n), lvl-1, tol_or_rank);
        F.maxrk = max([F.A12.r F.A21.r F.A11.maxrk F.A22.maxrk]);
    end
end