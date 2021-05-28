function [F, rank] = hifie2my(A,x,occ,rank_or_tol)
  [F, rank] = hifie2_base(A,x,occ,rank_or_tol,@hifie_id);
end
