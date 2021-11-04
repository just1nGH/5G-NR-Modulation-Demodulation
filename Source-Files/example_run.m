addpath('Functions')  
msg = randi([0 1],100,1);
symb = nrModuMapper(msg,'16QAM');
N0 = 0.5;
rxsymb = symb + sqrt(N0/2)*randn(size(symb)) ;
msg_hat = nrSoftModuDemapper(rxsymb,'16QAM',N0,'max-log-map');
numErr = sum(msg ~=(msg_hat < 0))