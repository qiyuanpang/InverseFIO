function u=five_solve(N,Tree, f)
  
  P = 8; %smallest system
  L = round(log(N/P)/log(2))+1;
  
  %-----------------
  u = f;
  for ell=1:L
    wid = 2^(ell-1)*P;
    nck = N/wid;
    for i=1:nck
      for j=1:nck
	INS = Tree{ell}{i,j}.INS;        BDS = Tree{ell}{i,j}.BDS;
	invA = Tree{ell}{i,j}.invA;        BinvA = Tree{ell}{i,j}.BinvA;
	S = Tree{ell}{i,j}.S;
	%load
        in = u(INS,:);
        bd = u(BDS,:);
        %L
        bd = bd - BinvA*in;
        %solve
        in = invA*in;
        if(ell==L)
          bd = S*bd; %LEXING
        end
        %save
        u(INS,:) = in;
        u(BDS,:) = bd;
      end
    end
  end
  %-----------------
  for ell=L:-1:1
    wid = 2^(ell-1)*P;
    nck = N/wid;
    for i=1:nck
      for j=1:nck
	INS = Tree{ell}{i,j}.INS;        BDS = Tree{ell}{i,j}.BDS;
	invA = Tree{ell}{i,j}.invA;        BinvA = Tree{ell}{i,j}.BinvA;
	S = Tree{ell}{i,j}.S;
	%load
        in = u(INS,:);
        bd = u(BDS,:);
        %Lt
        in = in - transpose(BinvA)*bd;
        %save
        u(INS,:) = in;
        u(BDS,:) = bd;
      end
    end
  end
  
  
