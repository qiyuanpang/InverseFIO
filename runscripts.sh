rm -rf ./results/hqr/1D/kernel1/results.log;
nohup matlab -nodisplay -nodesktop -r FIOs_1D_inverse_qr > ./results/hqr/1D/kernel1/results.log < /dev/null &;
