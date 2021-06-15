#nohup nice -19 matlab -nodisplay -nodesktop -nosplash < testnorm2single.m > testnorm2_output_single_256_svd3L4_prof.txt  & 
# For profiling
nohup nice -19 matlab -nodisplay -nodesktop -nojvm -nosplash < maintest.m > res.txt 2>&1 &
# For non-profiling
