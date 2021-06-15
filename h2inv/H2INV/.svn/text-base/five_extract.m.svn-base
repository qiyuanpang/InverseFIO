function M = five_extractall(N, Tree)
%N: system size is N*N, e.g. N = 512;
%V: the external potential of size (N,N). V can be real or complex. 
%   e.g. V = 1+rand(N,N);
%h: the effective mesh size. e.g. h=1
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the diagonal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  P = 8; %smallest system
  L = round(log(N/P)/log(2))+1;
  
  Lmat = cell(L,1);
  Lmat{L} = cell(1,1);
  Lmat{L}{1,1} = Tree{L}{1,1}.S; %already inverted
  %  clear {Tree}{L}{1,1}.S;
  if(1)
    ell=L;
    wid = 2^(ell-1) * P;
    nck = N/wid;
    Lmat{ell-1} = cell(nck*2,nck*2);
    for i=1:nck
      for j=1:nck
        invA = Tree{ell}{i,j}.invA;        
        BinvA = Tree{ell}{i,j}.BinvA;
        
        invS = Lmat{ell}{i,j};
        
        invA(find(abs(invA)<eps)) = 0;
        BinvA(find(abs(BinvA)<eps)) = 0;
        invS(find(abs(invS)<eps)) = 0;

        G12 = -transpose(BinvA) * invS;
        G12(find(abs(G12)<eps)) = 0;
        G21 = transpose(G12);
        G11 = invA + (-G12)*BinvA;
        G22 = invS;

        G = [G11 G12; G21 G22];        

	%order
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
        
        idx = [d1 s1 d2 r1 d0 r4(end:-1:1) d8 s8];
        Lmat{ell-1}{2*i-1,2*j-1} = G(idx,idx);
        idx = [d2 s2 d3 s3 d4 r2 d0 r1(end:-1:1)];
        Lmat{ell-1}{2*i  ,2*j-1} = G(idx,idx);
        idx = [d8 r4 d0 r3(end:-1:1) d6 s6 d7 s7];
        Lmat{ell-1}{2*i-1,2*j  } = G(idx,idx);
        idx = [d0 r2(end:-1:1) d4 s4 d5 s5 d6 r3];
        Lmat{ell-1}{2*i  ,2*j  } = G(idx,idx);
      end
    end
  end
  %---
  for ell=L-1:-1:2
    wid = 2^(ell-1) * P;
    nck = N/wid;
    Lmat{ell-1} = cell(nck*2,nck*2);
    for i=1:nck
      for j=1:nck
        invA = Tree{ell}{i,j}.invA;        
        BinvA = Tree{ell}{i,j}.BinvA;
        
        invS = Lmat{ell}{i,j};

        invA(find(abs(invA)<eps)) = 0;
        BinvA(find(abs(BinvA)<eps)) = 0;
        invS(find(abs(invS)<eps)) = 0;
       
        G12 = - transpose(BinvA) * invS;
	G12(find(abs(G12)<eps)) = 0;
        G21 = transpose(G12);
        G11 = invA + (-G12)*BinvA;
        G22 = invS;

        G = [G11 G12; G21 G22];        

        %order
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
        
        idx = [d1 s1 d2 r1 d0 r4(end:-1:1) d8 s8];
        Lmat{ell-1}{2*i-1,2*j-1} = G(idx,idx);
        idx = [d2 s2 d3 s3 d4 r2 d0 r1(end:-1:1)];
        Lmat{ell-1}{2*i  ,2*j-1} = G(idx,idx);
        idx = [d8 r4 d0 r3(end:-1:1) d6 s6 d7 s7];
        Lmat{ell-1}{2*i-1,2*j  } = G(idx,idx);
        idx = [d0 r2(end:-1:1) d4 s4 d5 s5 d6 r3];
        Lmat{ell-1}{2*i  ,2*j  } = G(idx,idx);
      end
    end
  end
    
  if(1)
    M = zeros(N,N);
    ell = 1;
    wid = 2^(ell-1) * P;
    nck = N/wid;
    for i=1:nck
      for j=1:nck
        INS = Tree{ell}{i,j}.INS;        
        BDS = Tree{ell}{i,j}.BDS; 
        invA = Tree{ell}{i,j}.invA;
        BinvA = Tree{ell}{i,j}.BinvA;
        invS = Lmat{ell}{i,j};
        
        invA(find(abs(invA)<eps)) = 0;
        BinvA(find(abs(BinvA)<eps)) = 0;
        invS(find(abs(invS)<eps)) = 0;
        
        nbds = numel(BDS);
        nins = numel(INS);
        nttl = nbds + nins;
        G12 = -transpose(BinvA) * invS;
	G12(find(abs(G12)<eps)) = 0;
        G21 = transpose(G12);
        G11 = invA + (-G12)*BinvA;
        G22 = invS;
	% clear invA B TT Tree{ell}{i,j};
        G = [G11 G12; G21 G22];        
        % clear G12 G21 G11 G22;
        M([INS;BDS]) = diag(G);
      end
    end
  end
  
  
