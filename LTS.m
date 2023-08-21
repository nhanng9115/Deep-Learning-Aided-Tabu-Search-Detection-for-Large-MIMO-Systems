function result = LTS(y, H, M, P, modScheme, Nt, xInit_c,...
    mOrder, nBit, xBin)
global sim_mod

nMul = 0; nAdd = 0; nComp = 0;
nOper_TS = [0,0,0,0];

H_pseudo = (H'*H)^(-1)*H';
x_bar = H_pseudo*y;
x_quan = xInit_c;

symset = loadAlphabet(modScheme);
delta = 0.25;
d_min = 2;
Th = delta*d_min;
% count
nMul = nMul + 1;

[~, U] = qr(H);

x_hat = zeros(Nt,1);
first_time = 1;
for k = Nt:-1:1
    sum = 0;
    for l = k+1 : Nt
        sum = sum + U(k,l)/U(k,k) * (x_hat(l) - x_bar(l)); % 1 add, 1 div, 1 mul, 1 sub
        nAdd = nAdd + 1;
        nMul = nMul + 1;
    end
    r_k = x_bar(k) - sum; % 1 sub
    nAdd = nAdd + 1;
    
    % find symbol in A closest to rk
    dis_tmp = abs(symset - r_k);
    [dis_min, ind_min] = min(dis_tmp);
    a_q = symset(ind_min);
    
    nComp = nComp + 1;
    if dis_min < Th % 1 comp
        x_hat(k) = a_q;
    else
        x_hat(k) = x_quan(k);
        
        % tabu search
        x_tild_init = x_hat(k:end); % init solution
        H_tild = U(k:end, k:end);
        y_tild = H_tild * x_bar(k:end); % k(k+1)/2 mul, k(k-1)/2 add
        % count: H_tild and x_bar dont change => only need to compute the
        % first time
        if first_time == 1
            nMul = nMul + k*(k+1)/2;
            nAdd = nAdd + k*(k-1)/2;
            first_time = 0;
        end
        
        % real system
        x_tild_init_r = [real(x_tild_init); imag(x_tild_init)];
        H_tild_r = [real(H_tild), -imag(H_tild); imag(H_tild), real(H_tild)];
        y_tild_r = [real(y_tild); imag(y_tild)];
        
        Nt_tild = Nt - k + 1; % 1 sub
        [xHatTS, nOperTS_tmp] = LTS_TS(y_tild_r, H_tild_r, M, P, modScheme, Nt_tild, x_tild_init_r);
        x_hat(k:end) = xHatTS;
        % count: convert to real mul. and add comp. of TS
        nOper_TS = nOper_TS + nOperTS_tmp;
    end
end
nMul_r = nMul*4;
nAdd_r = nAdd*2 + nMul*2;

if strcmp(sim_mod,'ber')
    x_hat_demod = qamdemod(x_hat,mOrder,0,'gray');
    x_hat_bin = de2bi(x_hat_demod, nBit);
    result = biterr(xBin,x_hat_bin);
else
    result = nMul_r + nAdd_r + nComp;
    % nOper = [nMul, nAdd, nComp, nTot];
end
end % end function