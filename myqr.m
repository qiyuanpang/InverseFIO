function [Q, R] = myqr(A,rank_or_tol)
    [m,n] = size(A);
    [~,R,E] = qr(A,0);
    if rank_or_tol < 1
        k = sum(abs(diag(R)) > abs(R(1))*rank_or_tol);
    else
        k = min(rank_or_tol,min(m,n));
    end
    sk = E(1:k);
    rd = E(k+1:end);
    T = R(1:k,1:k)\R(1:k,k+1:end);
    Q = A(:, sk);
    [~,I] = sort(E);
    R = [eye(k) T];
    R = R(:,I);
end