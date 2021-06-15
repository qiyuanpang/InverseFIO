% TESTOFFDIAG is a test subroutine.  First testGreen.m is run and the
% variables Kernel, DiagBlock, SampleList and RKTree has been
% constructed.  Then the diagonal blocks of Kernel is taken away and
% leaves the off-diagonal part Goffdiag.  This part is compared with
% RKTreeApplyH2 to verify the low rank matrix-vector multiplication
% subroutine.

Goffdiag = full(Kernel);
disp(' Generating off-diagonal matrix of G');
for i = 1 : length(DiagBlock)
  Goffdiag(DiagBlock{i}.indTarget, DiagBlock{i}.indSource) = 0.0;
end
disp(' Test (G_H2 - G_offdiag) * random vector');
u = zeros(N,N);
u(1,1) = 1;
v1 = RKTreeApplyH2( SampleList, RKTree, u, N, nLevel+1, 1 );
v2 = reshape(Goffdiag * reshape(u,[],1),N,N);
fprintf(' norm(v1-v2) = %25.10f\n', norm(v1 - v2));


% for i = 1 : 2: N
  % for j = 1 : 2: N
    % u = zeros(N,N);
    % u(i,j) = 1;
    % v1 = RKTreeApplyH2( SampleList, RKTree, u, N, nLevel+1, 1 );
    % v2 = reshape(Goffdiag * reshape(u,[],1),N,N);
    % fprintf(' i = %2d, j = %2d, norm(v1-v2) = %25.10f\n', i, j, ...
      % norm(v1 - v2));
  % end
% end
