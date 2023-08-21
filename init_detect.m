function [x_init, n_oper] = init_detect(y, H, sigma2, modScheme, detector)
N = size(H,2);

switch detector
    case 'ZF'
        %% ZF
        W = (H'*H)^(-1) * H';
        x_hat = W*y;
        n_mul = 1/3*N^3 + 3/2*N^2 + 1/6*N;
        n_add = 1/3*N^3 + 3/2*N^2 + 1/6*N - 1;
        
    case 'MMSE'
        %% MMSE
        W = (H'*H + sigma2*eye(N))^(-1) * H';
        x_hat = W*y;
        
        n_add = 3*N^3 - N^2 - 2*N + 2;
        n_mul = 3*N^3 + N^2 - 2*N;

        n_add = n_add + N^2 - N;
        n_mul = n_mul + N^2;
        
    case 'OSIC'
        OSIC_type = 1;
        [x_hat, n_add, n_mul] = OSIC_detector(y,H,sigma2,N,OSIC_type,modScheme);
end
% slice
x_init = zeros(N,1);
for kk = 1:N
    x_init(kk) = mod_slicer(x_hat(kk), modScheme);
end
n_oper = n_add + n_mul;
end % eof
