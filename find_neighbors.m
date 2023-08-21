function [nb_set, avoid_nb] = find_neighbors(cand, modScheme, tbList, Nt)

[constel, ~] = loadAlphabet(modScheme);
tbList_row = tbList.';
avoid_nb = zeros(Nt,1);

if strcmp(modScheme,'QPSK')
    K = 2;
else
    K = 4;
end

Neigh = repmat(cand,1,K*Nt);
for i=1:Nt
    dis = abs(constel - cand(i));
    idx_2 = (dis==2); % 6 for 16-QAM
    dis(idx_2 == 0) = Inf;
    
    [dis_incr,idx_incr] = sort(dis);
    
    for nn = 1:K
        if dis_incr(nn) == 2
            Neigh(i,i+(nn-1)*Nt) = constel(idx_incr(nn));
            if ~isempty(tbList_row)
                [~,idx] = ismember(Neigh(:,i+(nn-1)*Nt).',tbList_row, 'rows');
                if idx ~= 0
                    Neigh(:,i+(nn-1)*Nt) = zeros(Nt,1);
                end
            end
        else
            Neigh(:,i+(nn-1)*Nt) = zeros(Nt,1);
        end
    end
end

nb_set_c = Neigh(:,any(Neigh));
nb_set = [real(nb_set_c); imag(nb_set_c)];
end % eof