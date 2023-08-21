function [ber, com] = TS(y, H, M, P, m_scheme, s_init, m_order, n_bit, x_bin)
% global test
N = size(H,2);
K = size(H,1);
Nt = N/2;
% counting number of complex additions, multiplications, comparisons
n_mul = 0; n_add = 0; n_comp = 0;

cand = s_init;
cand_c = cand(1:Nt) + 1i*cand(Nt+1:end);
tb_list = cand_c; % Put the initial solution to tabu list
dis_cand = y - H*cand;
tb_metric = norm(dis_cand)^2; % = disCand'*disCand => need Nr muls and Nr-1 adds
    
%% test
x_hat = cand;
m_hat = tb_metric;
p = 1;
count_iter = 0;
for mm = 1:M
    
    if count_iter > ceil(M*1)
        break
    end
    
    %% Find neighbor set
    nb_set = find_neighbors(cand_c, m_scheme,tb_list, Nt);
    L = size(nb_set, 2);
    if L == 0
        break
    end
    list_m = [];
    list_dis = [];
    for ii = 1:L
        x = nb_set(:,ii);
        d = find(x - cand);
        delta_x = cand(d) - x(d); % 1 sub
        dis = dis_cand + H(:,d).*delta_x; % K muls and K adds
        m_x = norm(dis)^2; % K muls and K-1 adds
        % counting
        n_mul = n_mul + 2*K;
        n_add = n_add + 2*K-1;
        
        list_m = cat(2,list_m,m_x);
        list_dis =  cat(2,list_dis,dis);
    end
    
    %% Find the best neighbor
    [m_xbest, i_min] = min(list_m); % L comps
    n_comp = n_comp + L;
    x_best = nb_set(:,i_min);
    cand = x_best;
    dis_cand = list_dis(:,i_min); % update distance of this candidate
    
    %% test
    p = p + 1; % increase the tabu lenth
    
    %% save the best solution
    if m_xbest < m_hat % 1 comp
        x_hat = x_best;
        m_hat = m_xbest;
        count_iter = 0;
    else
        count_iter = count_iter + 1;
    end
    n_comp = n_comp + 1;
    
    %% release the first solution if tabulist is full and push a new solution
    if p > P
        tb_list = tb_list(:,2:end);
        tb_metric = tb_metric(:,2:end);
        p = p - 1;
    end
    cand_c = cand(1:Nt) + 1i*cand(Nt+1:end);
    tb_list = cat(2, tb_list, cand_c);
    tb_metric = cat(2, tb_metric, m_xbest);
end

x_hat = x_hat(1:Nt) + 1i*x_hat(Nt+1:end);
x_hat_demod = qamdemod(x_hat,m_order);
x_hat_bin = de2bi(x_hat_demod, n_bit);
ber = biterr(x_bin,x_hat_bin); % ber

com = n_mul + n_add; % complexity
end % end function