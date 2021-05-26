% RSKELF_DIAG  Extract diagonal using recursive skeletonization factorization
%              via matrix unfolding.
%
%    This is a variant of the selected inversion algorithm for sparse matrices
%    via the multifrontal factorization, generalized to account for additional
%    ID-based sparsification operators. The point is that only a sparse subset
%    of all matrix entries need to be reconstructed from the top-level
%    skeletons in order to compute the diagonal.
%
%    Typical complexity: same as RSKELF.
%
%    D = RSKELF_DIAG(F) produces the diagonal D of the factored matrix F.
%
%    D = RSKELF_DIAG(F,DINV) computes D = DIAG(F) if DINV = 0 (default) and
%    D = DIAG(INV(F)) if DINV = 1.
%
%    D = RSKELF_DIAG(F,DINV,OPTS) also passes various options to the algorithm.
%    Valid options include:
%
%      - VERB: display status info if VERB = 1 (default: VERB = 0). This prints
%              to screen a table tracking extraction statistics through
%              factorization level (i.e., tree leaves are at level 1). Special
%              levels: 'A', for determining all required entries to compute.
%
%    Related references:
%
%      L. Lin, J. Lu, L. Ying, R. Car, W. E. Fast algorithm for extracting the
%        diagonal of the inverse matrix with application to the electronic
%        structure analysis of metallic systems. Commun. Math. Sci. 7 (3):
%        755-777, 2009.
%
%    See also MF_DIAG, RSKELF, RSKELF_SPDIAG.

function D = rskelf_diag(F,dinv,opts)

  % set default parameters
  if nargin < 2, dinv = 0; end
  if nargin < 3, opts = []; end
  if ~isfield(opts,'verb'), opts.verb = 0; end

  % print header
  if opts.verb
    fprintf([repmat('-',1,31) '\n'])
    fprintf('%3s | %12s | %10s\n','lvl','nnz kept','time (s)')
    fprintf([repmat('-',1,31) '\n'])
  end

  % initialize
  N = F.N;
  nlvl = F.nlvl;
  rem = true(N,1);  % which points remain?
  nz = N;           % initial capacity for sparse matrix workspace
  I = zeros(nz,1);
  J = zeros(nz,1);

  % find required entries at each level
  ts = tic;
  keep = cell(nlvl,1);  % entries to keep after unfolding at each level
  keep{1} = sparse(1:N,1:N,true(N,1),N,N);  % at leaf, just need diagonal
  for lvl = 1:nlvl-1  % current level is lvl+1
    nz = 0;

    % eliminate redundant indices
    rem([F.factors(F.lvp(lvl)+1:F.lvp(lvl+1)).rd]) = false;

    % keep entries needed directly by previous level
    [I_,J_] = find(keep{lvl});
    idx = rem(I_) & rem(J_);
    [I,J,nz] = sppush2(I,J,nz,I_(idx),J_(idx));

    % loop over nodes at previous level
    for i = F.lvp(lvl)+1:F.lvp(lvl+1)

      % keep skeleton entries
      sk = F.factors(i).sk;
      [I_,J_] = ndgrid(sk);
      [I,J,nz] = sppush2(I,J,nz,I_,J_);
    end

    % construct requirement matrix
    idx = 1:nz;
    if F.symm ~= 'n', idx = find(I(idx) >= J(idx)); end
    if isoctave(), keep{lvl+1} = sparse(I(idx),J(idx),true(size(idx)),N,N);
    else, keep{lvl+1} = logical(sparse(I(idx),J(idx),ones(size(idx)),N,N));
    end
  end
  t = toc(ts);

  % print summary
  if opts.verb
    keep_ = keep{1};
    for lvl = 1:nlvl-1, keep_ = keep_ | keep{lvl+1}; end
    fprintf('%3s | %12d | %10.2e\n','a',nnz(keep_),t)
  end

  % unfold factorization
  V = zeros(length(I),1);
  M = sparse(N,N);  % successively unfolded matrix
  for lvl = nlvl:-1:1  % loop from top-down
    ts = tic;

    % find all existing entries
    [I_,J_,V_] = find(M);
    nz = length(V_);
    I(1:nz) = I_;
    J(1:nz) = J_;
    V(1:nz) = V_;

    % loop over nodes
    for i = F.lvp(lvl)+1:F.lvp(lvl+1)
      sk = F.factors(i).sk;
      rd = F.factors(i).rd;
      nsk = length(sk);
      nrd = length(rd);

      T = F.factors(i).T;
      L = F.factors(i).L;
      p = F.factors(i).p;
      E = F.factors(i).E;
      if F.symm == 'n' || F.symm == 's'
        U = F.factors(i).U;
        G = F.factors(i).F;
      else
        U = F.factors(i).L';
        G = F.factors(i).E';
      end

      % unfold local factorization
      ird = 1:nrd;
      isk = nrd+(1:nsk);
      X = zeros(nrd+nsk);
      % redundant part
      if F.symm == 'h'
        if dinv, X(ird,ird) = inv(F.factors(i).U);
        else,    X(ird,ird) =     F.factors(i).U ;
        end
      else,      X(ird,ird) = eye(nrd);
      end
      % skeleton part
      Xsk = spsymm(spget(M,sk,sk),F.symm);
      X(isk,isk) = Xsk;
      % undo elimination and sparsification
      if dinv
        X(:,ird) = (X(:,ird) - X(:,isk)*E)/L;
        X(ird,:) = U\(X(ird,:) - G*X(isk,:));
        if ~isempty(p)
          X(:,ird(p)) = X(:,ird);
          if F.symm == 'h', X(ird(p),:) = X(ird,:); end
        end
        if F.symm == 's', X(:,isk) = X(:,isk) - X(:,ird)*T.';
        else,             X(:,isk) = X(:,isk) - X(:,ird)*T' ;
        end
        X(isk,:) = X(isk,:) - T*X(ird,:);
      else
        X(:,isk) = X(:,isk) + X(:,ird)*G;
        X(isk,:) = X(isk,:) + E*X(ird,:);
        X(:,ird) = X(:,ird)*U;
        X(ird,:) = L*X(ird,:);
        if ~isempty(p)
          if F.symm == 'h', X(:,ird(p)) = X(:,ird); end
          X(ird(p),:) = X(ird,:);
        end
        X(:,ird) = X(:,ird) + X(:,isk)*T;
        if F.symm == 's', X(ird,:) = X(ird,:) + T.'*X(isk,:);
        else,             X(ird,:) = X(ird,:) + T' *X(isk,:);
        end
      end
      X(isk,isk) = X(isk,isk) - Xsk;  % to be stored as update

      % store update to global sparse matrix
      [I_,J_] = ndgrid([rd sk]);
      [I,J,V,nz] = sppush3(I,J,V,nz,I_,J_,X);
    end

    % update unfolded sparse matrix
    M = sparse(I(1:nz),J(1:nz),V(1:nz),N,N) .* keep{lvl};
    t = toc(ts);

    % print summary
    if opts.verb, fprintf('%3d | %12d | %10.2e\n',lvl,nnz(keep{lvl}),t); end
  end

  % finish
  D = spdiags(M,0);
  if opts.verb, fprintf([repmat('-',1,31) '\n']); end
end