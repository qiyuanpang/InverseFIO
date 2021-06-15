function Tree = five_setup(V,N, h)
%N: system size is N*N, e.g. N = 512;
%V: the external potential of size (N,N). V can be real or complex. 
%   e.g. V = 1+rand(N,N);
%h: the effective mesh size. e.g. h=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the operator -(1/h^2) (Discrete Laplace) + V
% as well as the hierachical Schur complements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = 8; %smallest system
L = round(log(N/P)/log(2))+1;
%h = 1/N;
GRD = reshape(1:(N*N),N,N);
Tree = cell(L,1);
%----
if(1)
  ell = 1;
  wid = 2^(ell-1) * P;
  nck = N/wid;
  Tree{ell} = cell(nck,nck);
  %construct the matrix to be used
  sz = (P+1)^2;
  idx = reshape(1:(P+1)^2, P+1, P+1);
  low = 1:P;    upp = 2:P+1;    all = 1:P+1;
  mid = 2:P;    sss = 1;    eee = P+1;
  cur = zeros(sz,sz);
  hoh = [1/2, ones(1,P-1), 1/2];
  cur = cur + sparse(idx(low,all),idx(upp,all),-1*ones(P,1)*hoh,sz,sz);
  cur = cur + sparse(idx(upp,all),idx(low,all),-1*ones(P,1)*hoh,sz,sz);
  cur = cur + sparse(idx(all,low),idx(all,upp),-1*hoh'*ones(1,P),sz,sz);
  cur = cur + sparse(idx(all,upp),idx(all,low),-1*hoh'* ones(1,P),sz,sz);
  cur = cur + sparse(idx(all,all),idx(all,all),4*hoh'*hoh,sz,sz);
  cur = 1/(h^2) * cur;
  %construct the matrix
  for i=1:nck
    for j=1:nck
      %construct matrix
      is = (i-1)*P+[1:P+1]; is(is==N+1)=1;
      js = (j-1)*P+[1:P+1]; js(js==N+1)=1;
      coef = [0.5, ones(1,P-1), 0.5];
      mask = coef'*coef;
      vv = V(is,js) .* mask;
      tmp = cur + diag(vv(:));
      %reorganize
      idx = reshape(1:(P+1)^2, P+1,P+1);
      bd = [idx(1:end-1,1)' idx(end,1:end-1) idx(end:-1:2,end)' idx(1,end:-1:2)];
      in = idx(2:P,2:P);        in = in(:)';
      ord = [in bd];
      tmp = tmp(ord,ord);
      idx = GRD(is,js);
      BDS = [idx(1:end-1,1)' idx(end,1:end-1) idx(end:-1:2,end)' idx(1,end:-1:2)];
      INS = idx(2:P,2:P);        INS = INS(:)';
      if(numel(INS)~=numel(in) || numel(BDS)~=numel(bd)) error('wrong'); end;
	%Schur complement
	nin = 1:numel(in);
	nbd = numel(in) + [1:numel(bd)];

	invA = inv(full(tmp(nin,nin))+eps*rand(numel(in), numel(in)));
	invA(find(abs(invA)<eps)) = 0;
	B = tmp(nbd,nin);
	B(find(abs(B)<eps)) = 0;
	BinvA = B*invA;
	BinvA(find(abs(BinvA)<eps)) = 0;
	S = tmp(nbd,nbd) - BinvA*transpose(B);
	Tree{ell}{i,j} = struct('INS', INS', 'BDS', BDS', 'invA', ...
	  invA, 'BinvA', BinvA, ...
	  'S', S); clear invA B BinvA S;
      end
    end
  end
  %----
  for ell=2:L-1
    wid = 2^(ell-1) * P;
    nck = N/wid;
    Tree{ell} = cell(nck,nck);
    for i=1:nck
      for j=1:nck
	%construct matrix
	is = (i-1)*wid+[1:wid+1]; is(is==N+1)=1;
	js = (j-1)*wid+[1:wid+1]; js(js==N+1)=1;
	Q = wid/2;
	cnt = 0;
	d0 = cnt+[1];        cnt = cnt+1;
	r1 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	r2 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	r3 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	r4 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d1 = cnt+[1];        cnt = cnt+1;
	s1 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d2 = cnt+[1];        cnt = cnt+1;
	s2 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d3 = cnt+[1];        cnt = cnt+1;
	s3 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d4 = cnt+[1];        cnt = cnt+1;
	s4 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d5 = cnt+[1];        cnt = cnt+1;
	s5 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d6 = cnt+[1];        cnt = cnt+1;
	s6 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d7 = cnt+[1];        cnt = cnt+1;
	s7 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	d8 = cnt+[1];        cnt = cnt+1;
	s8 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	tmp = zeros(12*Q-3, 12*Q-3);
	idx = [d1 s1 d2 r1 d0 r4(end:-1:1) d8 s8];
	tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i-1,2*j-1}.S; Tree{ell-1}{2*i-1,2*j-1}.S=[];
	idx = [d2 s2 d3 s3 d4 r2 d0 r1(end:-1:1)];
	tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i  ,2*j-1}.S; Tree{ell-1}{2*i  ,2*j-1}.S=[];
	idx = [d8 r4 d0 r3(end:-1:1) d6 s6 d7 s7];
	tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i-1,2*j  }.S; Tree{ell-1}{2*i-1,2*j  }.S=[];
	idx = [d0 r2(end:-1:1) d4 s4 d5 s5 d6 r3];
	tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i  ,2*j  }.S; Tree{ell-1}{2*i  ,2*j  }.S=[];
	bd = [d1 s1 d2 s2 d3 s3 d4 s4 d5 s5 d6 s6 d7 s7 d8 s8];
	in = [d0 r1 r2 r3 r4];
	ord = [in bd];
	%tmp = tmp(ord,ord);
	idx = GRD(is,js);
	BDS = [idx(1:end-1,1)' idx(end,1:end-1) idx(end:-1:2,end)' idx(1,end:-1:2)];
	INS = [idx(Q+1,Q+1) idx(Q+1,2:Q) idx(2*Q:-1:Q+2,Q+1)' idx(Q+1,2*Q:-1:Q+2) idx(2:Q,Q+1)'];
	if(numel(INS)~=numel(in) || numel(BDS)~=numel(bd)) error('wrong'); end;
	  %Schur complement
	  nin = 1:numel(in);
	  nbd = numel(in) + [1:numel(bd)];
	  invA = inv(full(tmp(nin,nin))+eps*rand(size(numel(in), numel(in))));
	  invA(find(abs(invA)<eps)) = 0;
	  B = tmp(nbd,nin);
	  B(find(abs(B)<eps)) = 0;
	  BinvA = B*invA;
	  BinvA(find(abs(BinvA)<eps)) = 0;

	  S = tmp(nbd,nbd) - BinvA*transpose(B);
	  Tree{ell}{i,j} = struct('INS', INS', 'BDS', BDS', 'invA', ...
	    invA, 'BinvA', BinvA, ...
	    'S', S); clear invA B BinvA S;
	end
      end
    end
    %----
    if(1)
      ell = L;
      %wid = 2^(ell-1) * P;
      %nck = N/wid;
      %Tree{ell} = cell(nck,nck);
      wid = N;
      nck = 1;
      Tree{ell} = cell(1,1);
      for i=1:nck
	for j=1:nck
	  is = (i-1)*wid+[1:wid+1]; is(is==N+1)=1;
	  js = (j-1)*wid+[1:wid+1]; js(js==N+1)=1;
	  Q = wid/2;
	  cnt = 0;
	  d0 = cnt+[1];        cnt = cnt+1;
	  r1 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  r2 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  r3 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  r4 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  d1 = cnt+[1];        cnt = cnt+1;
	  s1 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  d2 = cnt+[1];        cnt = cnt+1;
	  s2 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  d3 = d1;
	  s3 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  d4 = cnt+[1];        cnt = cnt+1;
	  s4 = cnt+[1:Q-1];        cnt = cnt+(Q-1);
	  d5 = d1;
	  s5 = s2(end:-1:1);
	  d6 = d2;
	  s6 = s1(end:-1:1);
	  d7 = d1;
	  s7 = s4(end:-1:1);
	  d8 = d4;
	  s8 = s3(end:-1:1);

	  tmp = zeros(8*Q-4, 8*Q-4);
	  idx = [d1 s1 d2 r1 d0 r4(end:-1:1) d8 s8];
	  tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i-1,2*j-1}.S; Tree{ell-1}{2*i-1,2*j-1}.S=[];
	  idx = [d2 s2 d3 s3 d4 r2 d0 r1(end:-1:1)];
	  tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i  ,2*j-1}.S; Tree{ell-1}{2*i  ,2*j-1}.S=[];
	  idx = [d8 r4 d0 r3(end:-1:1) d6 s6 d7 s7];
	  tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i-1,2*j  }.S; Tree{ell-1}{2*i-1,2*j  }.S=[];
	  idx = [d0 r2(end:-1:1) d4 s4 d5 s5 d6 r3];
	  tmp(idx,idx) = tmp(idx,idx) + Tree{ell-1}{2*i  ,2*j  }.S; Tree{ell-1}{2*i  ,2*j  }.S=[];
	  bd = [d1 s1 d2 s2 s3 d4 s4];
	  in = [d0 r1 r2 r3 r4];
	  ord = [in bd];
	  %tmp = tmp(ord,ord);
	  idx = GRD(is,js);
	  BDS = [idx(1:end-1,1)' idx(end,2:end-1)];
	  INS = [idx(Q+1,Q+1) idx(Q+1,2:Q) idx(2*Q:-1:Q+2,Q+1)' idx(Q+1,2*Q:-1:Q+2) idx(2:Q,Q+1)'];
	  if(numel(INS)~=numel(in) || numel(BDS)~=numel(bd)) error('wrong'); end;
	    %Schur complement
	    nin = 1:numel(in);
	    nbd = numel(in) + [1:numel(bd)];
	    tmp(find(abs(tmp)<eps)) = 0;
	    invA = inv(full(tmp(nin,nin))+eps*rand(numel(in), numel(in)));
	    invA(find(abs(invA)<eps)) = 0;
	    B = tmp(nbd,nin);
	    B(find(abs(B)<eps)) = 0;
	    BinvA = B*invA;
	    BinvA(find(abs(BinvA)<eps)) = 0;

	    S = tmp(nbd,nbd) - BinvA*transpose(B);
	    S = inv(full(S)+eps*rand(size(S))); 
	    Tree{ell}{i,j} = struct('INS', INS', 'BDS', BDS', 'invA', ...
	      invA, 'BinvA', BinvA, ...
	      'S', S); clear invA B BinvA S;
	  end
	end
      end

