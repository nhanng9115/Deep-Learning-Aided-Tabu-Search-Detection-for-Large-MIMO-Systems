
N = 64;
M = 128;
L = 15;

C_detnet_QPSK = (18*L+2*M-1)*N^2 + (2*M-1+L)*N
C_detnet_QAM = (42*L+2*M-1)*N^2 + (2*M-1+L)*N;

C_scnet = (2*M-1+2*L)*N^2 + (2*M-1+5*L)*N
C_fsnet = (2*M-1+2*L)*N^2 + (2*M-1+4*L)*N