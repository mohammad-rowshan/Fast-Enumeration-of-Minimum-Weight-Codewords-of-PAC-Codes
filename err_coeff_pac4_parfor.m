N = 2^8;            % Code length
R = 0.75;            % Code rate
K = N * R;          % code dimention: the number of information bits
design_snr_db = 4;  % for the optimized code
g = [1 0 1 1 0 1 1];%[1 0 1 1];% 0 1 1 0 1 1];%[1 0 0 0 0 1 1];%[1 0 1 0 1 0 1];%
I = construct_dega(design_snr_db, N, K); % the indices of K information bits
% Rem = [56,52];%[60,58,57]; For (64,32)
% Add = [22,25];%[27,29,30];
% I0 = I;
% j = 1;
% for i = Rem
%     idx = find(I==i);
%     I(idx) = Add(j);
%     j = j + 1;
% end
% I = sort(I);
% find_set_M1(3,[5,10],5)

tic
[d_min,A_dmin,A_dminPAC] = err_coeff(I,N,g) 
toc
%set_Ki = find_Ki(7,5)


% returns the minimum distance of the code and the number of codewords with minumum distance (error coefficient)
function [dmin, A_dmin, A_dminPAC] = err_coeff(I,N,g)
    d = min(sum(dec2bin(I)-'0',2));
    dmin = 2^d; n = log2(N); A_dmin = 0; A_dminPAC = 0; A_dminPAC2 = 0;
    B = find(sum(dec2bin(I)-'0',2)==d); 
    set_F = find_F(I); % F indicator
    cosetsAdmin = zeros(length(B),4);
    cosetsAdmin(:,1) = I(B); 
    for i = B' % a set of indices of elements in I, not elements of I
        %i=find(I==205); % for jumoping to the test coset
        Ki_size = n - d;
        for x = find(dec2bin(I(i),n)-'0'==1) % support of I(i): positions of 1's in the binary representation of the sub-channel index
                ii = dec2bin(bitxor(N-1,I(i)),n)-'0';
                Ki_size = Ki_size + sum(ii(1:x-1));
        end
        
        %set_Kix = find_Kix(I(i),set_F,n); Ki_size = Ki_size - length(set_Kix); %For modified codes. To find rozen rows that could be member of Ki.
        A_dmin = A_dmin + 2^Ki_size;
        cosetsAdmin(find(B==i),2) = 2^Ki_size;
        A_dminPAC = A_dminPAC + 2^Ki_size;
        set_Fxi = find_Fxi(I(i),set_F,n); % elements in F that differ more than one element with I(i)
        if isempty(set_Fxi)
            A_dminPAC2 = A_dminPAC2 + 2^Ki_size;
            continue;
        end
        i_lead = I(i);
        max_f = max(set_Fxi);
        set_KiRed = find_KiRed(i_lead,max_f,n); %reduced Ki set up to max(Fxi)
        set_KiRed_gen = find_KiRed_gen(i_lead,max_f,n); %generalized

        % Subtracting the removed min_weight codewords by precoding
        %M_subKiRed = zeros(2^length(set_KiRed), max_f-i_lead);
        %M_subKiRed = find_ones_full(i_lead,max_f,set_KiRed_gen,n);

        cnt_comb = 0; cnt_comb_1 = 0; cnt_comb_2 = 0;
        A_dmin_minus = 0;
        A_Kif = 0;
        A_dmin_minus_f = 0;
        for x = 0:length(set_KiRed)
            A = nchoosek(set_KiRed,x); %each row, one combination
            if x==0
                A=[0];
                %A=[A,0];
            end
            %[n_combs,~] = size(A); %size(A,1)
            for y = 1:size(A,1) %parfor
                discounted = false;
                cnt_comb = cnt_comb + 1;
                u_supp = [i_lead];
                J = [];
                F1 = []; % the idex set of positions where u_i=1 in Ki
                F1xi = []; % the idex set of positions where u_i=1 not in Ki
                v = zeros(N,1); u = zeros(N,1);
                cur_state = zeros(1,length(g)-1);
                v(i_lead+1) = 1; 
                [u(i_lead+1),cur_state] = conv_1bit_nxtState(v(i_lead+1),cur_state,g); %the i-th bit first, +1 is due to starting index in MATLAB.
                if A(y,:)==0
                    J0 = []; % different with moving J.
                else 
                    J0 = A(y,:);
                end
                %[M_subKiRed, ones] = find_ones(i_lead,max_f,J0,set_KiRed,M_subKiRed,n);
                %ones = M_subKiRed(get_idx_MsubKiRed(J0,set_KiRed_gen),:); % retruns partial set M only for Ki, not Ki_gen.
                ones = zeros(1,max_f-i_lead);
                M = [];
                cnt_i = 1; cnt_f = 1;
                for j = i_lead+1:max_f
                    if set_F(j+1) == 0
                        if A(y,cnt_i)==j
                            u_supp = [u_supp j];
                            %ones = M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:);
                            %ones = bitxor(ones, find_ones_Fx(i_lead,max_f,J,[j],n));
                            [M, ones1] = find_ones_1R(i_lead,max_f,J,M,j,n);
                            ones = bitxor(ones, ones1);
                            J = [J j];
                            if conv_1bit(v(j+1),cur_state,g) == 0
                                v(j+1) = 1;
                            end
                            %[u(j+1),cur_state] = conv_1bit_nxtState(v(j+1),cur_state,g);
                            if cnt_i<x
                                cnt_i = cnt_i + 1;
                            end
                        else
                            if ones(j-i_lead)==1
                                u_supp = [u_supp j];
                                if conv_1bit(v(j+1),cur_state,g) == 0
                                    v(j+1) = 1;
                                end
                            else
                                if conv_1bit(v(j+1),cur_state,g) == 1
                                    v(j+1) = 1;
                                end

                            end
%                             if conv_1bit(v(j+1),cur_state,g) == 1 && x>1  
%                                 if sum(dec2bin(bitxor(i_lead,bitor(j,i_lead)),n)-'0'==1)>1
%                                     if ((length(J)>1  && ismember(j,find_set_M(i_lead,J,n)) ~= 1) || length(J)<=1)
%                                         v(j+1) = 1;
%                                     else
%                                         u_supp = [u_supp j];
%                                     end
%                                 else
%                                     v(j+1) = 1;
%                                 end
%                             elseif sum(dec2bin(bitxor(i_lead,bitor(j,i_lead)),n)-'0'==1)>1 && x>1
%                                 if (length(J)>1  && ismember(j,find_set_M(i_lead,J,n)) == 1)
%                                     v(j+1) = 1;
%                                     u_supp = [u_supp j];
%                                 end
%                             end
                            %[u(j+1),cur_state] = conv_1bit_nxtState(v(j+1),cur_state,g);
                        end
                        [u(j+1),cur_state] = conv_1bit_nxtState(v(j+1),cur_state,g);
                    else
                        [u(j+1),cur_state] = conv_1bit_nxtState(v(j+1),cur_state,g);
                        if set_Fxi(cnt_f)==j
                            if u(j+1)==1
                                F1xi = [F1xi j];
                                u_supp = [u_supp j];
                                %set_M0 = find_set_M(i_lead,J,n);
                                %if isempty(F1) || ((length(J)>0 && ~isempty(F1)) && ismember(j,find_set_M(i_lead,J,n)) ~= 1)
                                if ones(j-i_lead)==0
                                    %discounted = true;
                                    %cnt_comb_1 = cnt_comb_1 + 1;
                                    A_dmin_minus = A_dmin_minus + 2^(Ki_size-length(set_KiRed));
                                    A_Kif = A_Kif + 1;
                                    break;
                                end
                            else
                                if ones(j-i_lead)==1
                                    %discounted = true;
                                    %cnt_comb_1 = cnt_comb_1 + 1;
                                    A_dmin_minus = A_dmin_minus + 2^(Ki_size-length(set_KiRed));
                                    A_Kif = A_Kif + 1;
                                    break;
                                end
                            end
                            if cnt_f<length(set_Fxi)
                                cnt_f = cnt_f + 1;
                            end
                        else
                            if u(j+1)==1
                                F1 = [F1 j];
                                %ones = bitxor(ones, find_ones_Fx(i_lead,max_f,J,[j],n));
                                [M, ones1] = find_ones_1R(i_lead,max_f,J,M,j,n);
                                ones = bitxor(ones, ones1);
                                J = [J j];
                                u_supp = [u_supp j];
                                %ones = M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:);
                            end
                        end
                    end
                end
                %%{
                % intersect(set_M,J) always empty?It is possible when |J|=1
%                 if discounted == false %&& length(J)>1 % After max_f, J might get larger by including more from F and then it might require Fxi members as well
%                     set_M = find_set_M(i_lead,J,n);
%                     F11 = union(F1,F1xi);
%                     %len_M = length(set_M)
%                     %len_M_I = length(intersect(set_M,I))
%                     if length(intersect(set_M,set_Fxi)) > length(F1xi) || (isempty(intersect(set_M,F1xi)) && ~isempty(F1xi)) % intersection with F
%                     %if length(intersect(set_M,set_Fxi)) > length(F1xi) || (isempty(intersect(set_M,F1xi)) && ~isempty(F1xi)) % intersection with F
%                     %if length(intersect(set_M,I))<length(set_M) &&
%                     %length(intersect(set_M,J)) ~= length(set_M) %
%                     %intersection with F, intersect(set_M,J) always empty?
%                     %It is possible when |J|=1
%                         %if nnz(set_M>max_f)==0
%                         A_dmin_minus_f = A_dmin_minus_f + 2^(Ki_size-length(set_KiRed));
%                         discounted = true;
%                         cnt_comb_2 = cnt_comb_2 + 1;
%                         %end
%                     end
%                 end
%                 if discounted == false
%                     discounted = false;
%                 end
            end
        end
        %A_Kif
        %set_KiRed
        A_dminPAC2 = A_dminPAC2 + (2^length(set_KiRed) - A_Kif)  * 2^(Ki_size-length(set_KiRed));
        A_dminPAC = A_dminPAC - A_dmin_minus - A_dmin_minus_f;
        cosetsAdmin(find(B==i),3) = A_dmin_minus + A_dmin_minus_f;
        cosetsAdmin(find(B==i),4) = cosetsAdmin(find(B==i),2) - cosetsAdmin(find(B==i),3);
        %i_lead
        %Ki_size
        %A_dmin_minus
        %A_dmin_minus_f
        %A_dminPAC

    end
    cosetsAdmin(:,4) = cosetsAdmin(:,2) - cosetsAdmin(:,3);
    cosetsAdmin
    A_dminPAC2
end



function set_Fxi = find_Fxi(index,F,n) 
    i = index + 2;
    set_Fxi = [];
    while i < 2^n
        if F(i) == 1 && sum(dec2bin(bitxor(index,bitor(i-1,index)),n)-'0'==1)>1
            set_Fxi = [set_Fxi, i-1];
        end
        i = i + 1;
    end

end

function set_Kix = find_Kix(index,F,n) 
    i = index + 2;
    set_Kix = [];
    while i < 2^n
        if F(i) == 1 && sum(dec2bin(bitxor(index,bitor(i-1,index)),n)-'0'==1)==1
            set_Kix = [set_Kix, i-1];
        end
        i = i + 1;
    end

end


function set_F = find_F(I) % indicator vector, not index vector
    N = max(I)+1;
    set_F = ones(N,1);
    j=1;
    for i = 1:N
        if I(j)==i-1
            set_F(i) = 0;
            j = j + 1;
        end
    end
end


function set_KiRed = find_KiRed(index,max_Fxi,n)
    set_KiRed = [];
    Ki = find_Ki(index,n);
    for j = 1:length(Ki)
        if Ki(j) < max_Fxi
            set_KiRed = [set_KiRed, Ki(j)];
        end
    end
end

function set_KiRed_gen = find_KiRed_gen(i_lead,max_f,n)
    set_KiRed_gen = [];
    for j = i_lead+1:max_f
        if sum(dec2bin(bitxor(i_lead,bitor(j,i_lead)),n)-'0'==1) == 1
            set_KiRed_gen = [set_KiRed_gen, j];
        end
    end
end


function set_Ki = find_Ki(index,n)
    set_Ki = [];
    supp_index = find_supp(index,n);
    supp_zeros = find(dec2bin(bitxor(2^n-1,index),n)-'0'==1);  %complement of supp
    %Ki_size = n - length(supp_index); %due to signle-bit addition
    index_bin = dec2bin(index, n);
    for j = supp_zeros
        index_add = index_bin;
        index_add(j) = '1';
        set_Ki = [set_Ki, bin2dec(index_add)];
    end
    for x = supp_index
        for j = supp_zeros
            if j<x
                index_leftSW = index_bin;
                index_leftSW(x) = '0';
                index_leftSW(j) = '1';
                set_Ki = [set_Ki, bin2dec(index_leftSW)];
            end
        end
    end
    set_Ki = sort(set_Ki);
end



function M = find_set_M(index,J,n)
    M = [];
    supp_index = find_supp(index,n);
    supp_zeros = find(dec2bin(bitxor(2^n-1,index),n)-'0'==1);  %complement of supp
    
    %Sm_Si = set(supp_index)
    Sm_Ti = supp_zeros;
    
    for s = 2:length(J)
        row_sets = nchoosek(J, s);
        for row_idx = 1:size(row_sets,1)
            valid_set_for_M = true;
            Sm_S = supp_index;
            Sm_T = [];
            column_elmnt_cnt = zeros(length(supp_zeros),1);
            for r = row_sets(row_idx,:)
                r_bin = dec2bin(r,n);
                supp_r = find_supp(r,n);
                Sm_S = intersect(Sm_S, supp_r);
                if ~isempty(intersect(Sm_T,intersect(Sm_Ti, supp_r)))
                    valid_set_for_M = false;
                    break;
                end
                Sm_T = union(Sm_T,intersect(Sm_Ti, supp_r));
                cnt=1;
                for t = supp_zeros
                    if r_bin(t)=='1'  %%%%%%%% 1 WAS INCORRECT. %%%%%%%%%%
                        column_elmnt_cnt(cnt) = column_elmnt_cnt(cnt) + 1;
                        if column_elmnt_cnt(cnt) > 1
                            valid_set_for_M = false;
                            break;
                        end
                    end
                    cnt = cnt + 1;
                end
                if valid_set_for_M == false
                    break;
                end
            end
            if valid_set_for_M == true
                Sm = union(Sm_T,Sm_S);
                M = [M supp2dec(Sm,n)];
            end
        end
    end
    % Removeing the repeated elements by even times.
    M1 = unique(M);
    count_M = histc(M,M1);
    indices = [];
    for cnt = 1:length(count_M)
        if mod(count_M(cnt),2) == 0
            indices = [indices cnt];
        end
    end
    M1(indices) = [];
    M = M1;
end

%%
function [M_subKiRed,ones] = find_ones(i_lead,max_f,J,set_KiRed,M_subKiRed,n)
    ones = zeros(1,max_f-i_lead);
    x = size(J,1);
    if x == 2
        M  = find_set_M1(i_lead,J,n);
        if ~isempty(M)
            ones(1,M(1)-i_lead) = 1;
        end
    elseif x > 2
        M  = find_set_M1(i_lead,J,n);
        if ~isempty(M)
            ones(1,M(1)-i_lead) = 1;
        end
        B = nchoosek(J,x-1);
        for z = 1:size(B,1)
            ones = bitxor(ones, M_subKiRed(get_idx_MsubKiRed(B(z,:),set_KiRed),:));
        end
        M_subKiRed(get_idx_MsubKiRed(J,set_KiRed),:) = ones;
    end
    for z=1:x
        ones(1,J(z)-i_lead) = 1;
    end
end

function M_subKiRed = find_ones_full0(i_lead,max_f,set_KiRed_gen,n)
    M_subKiRed = zeros(2^length(set_KiRed_gen), max_f-i_lead);
    for x = 2:length(set_KiRed_gen)
        A = nchoosek(set_KiRed_gen,x); %each row, one combination
        for y = 1:size(A,1)
            J = A(y,:);
            %[M_subKiRed, ones] = find_ones(i_lead,max_f,J0,set_KiRed_gen,M_subKiRed,n);

            ones = zeros(1,max_f-i_lead);
            %ones = zeros(1,2^n);
            
            if x == 2
                M  = find_set_M1(i_lead,J,n);
                if ~isempty(M) && M(1) <= max_f
                    ones(1,M(1)-i_lead) = 1;
                    %M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:) = ones;
                end
            elseif x > 2
                M  = find_set_M1(i_lead,J,n);
                if ~isempty(M) && M(1) <= max_f
                    ones(1,M(1)-i_lead) = 1;
                end
                for xx=x-1:2
                    %B = nchoosek(J,x-1);
                    B = nchoosek(J,xx);
                    for z = 1:size(B,1)
                        ones = bitxor(ones, M_subKiRed(get_idx_MsubKiRed(B(z,:),set_KiRed_gen),:));
                    end
                end
            end
            M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:) = ones;
            %for z=1:x
                %ones(1,J(z)-i_lead) = 1;
            %end

        end
    end
    %%%%%%%%%%%%
end


function M_subKiRed = find_ones_full(i_lead,max_f,set_KiRed_gen,n)
    M_subKiRed_base = zeros(2^length(set_KiRed_gen), max_f-i_lead);
    M_subKiRed = zeros(2^length(set_KiRed_gen), max_f-i_lead);
    for x = 2:length(set_KiRed_gen)
        A = nchoosek(set_KiRed_gen,x); %each row, one combination
        for y = 1:size(A,1)
            J = A(y,:);
            %[M_subKiRed, ones] = find_ones(i_lead,max_f,J0,set_KiRed_gen,M_subKiRed,n);

            ones = zeros(1,max_f-i_lead);
            %ones = zeros(1,2^n);
            M  = find_set_M1(i_lead,J,n);
            if ~isempty(M) && M(1) <= max_f
                ones(1,M(1)-i_lead) = 1;
                %M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:) = ones;
            end
            
            M_subKiRed_base(get_idx_MsubKiRed(J,set_KiRed_gen),:) = ones;
            if x == 2
                M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:) = ones;
            end
        end
    end
    for x = 3:length(set_KiRed_gen)
        A = nchoosek(set_KiRed_gen,x); %each row, one combination
        for y = 1:size(A,1)
            J = A(y,:);
            %[M_subKiRed, ones] = find_ones(i_lead,max_f,J0,set_KiRed_gen,M_subKiRed,n);

            ones = zeros(1,max_f-i_lead);
            M  = find_set_M1(i_lead,J,n);
            if ~isempty(M) && M(1) <= max_f
                ones(1,M(1)-i_lead) = 1;
            end
            for xx=x-1:-1:2
                %B = nchoosek(J,x-1);
                B = nchoosek(J,xx);
                for z = 1:size(B,1)
                    ones = bitxor(ones, M_subKiRed_base(get_idx_MsubKiRed(B(z,:),set_KiRed_gen),:));
                end
            end
            M_subKiRed(get_idx_MsubKiRed(J,set_KiRed_gen),:) = ones;
            %for z=1:x
                %ones(1,J(z)-i_lead) = 1;
            %end

        end
    end
    %%%%%%%%%%%%
end

function ones = find_ones_Fx(i_lead,max_f,J,F1,n)
    ones = zeros(1,max_f-i_lead);
    for x = 1:length(F1)
        A = nchoosek(F1,x); %each row, one combination
        for y = 1:size(A,1)
            F = A(y,:);
            for xx=1:length(J)
                B = nchoosek(J,xx);
                for z = 1:size(B,1)
                    JF = union(B(z,:),F);
                    M  = find_set_M1(i_lead,JF,n);
                    if ~isempty(M) && M(1) <= max_f
                        ones(1,M(1)-i_lead) = bitxor(ones(1,M(1)-i_lead), 1);
                    end
                end
            end
        end
    end
end


function [M, ones] = find_ones_1R(i_lead,max_f,J,M,j,n) %Progressively
    ones = zeros(1,max_f-i_lead);
    for x=1:length(M)
        %Mj = union(M(x),j);
        m  = find_set_M1R(i_lead,M(x),j,n);
        if m>0 && m <= max_f
            M = [M m]; %Make sure M will not have a repetitive elements. It will cancelled in the next line anyways.
            ones(1,m-i_lead) = bitxor(ones(1,m-i_lead), 1);
        end
    end
    for x=1:length(J)
        Jj = union(J(x),j);
        m  = find_set_M1(i_lead,Jj,n);
        if m>0 && m <= max_f
            M = [M m];
            ones(1,m-i_lead) = bitxor(ones(1,m-i_lead), 1);
        end
    end
end



function m = find_set_M1(index,J,n)
    %M = [];
    m = 0;
    supp_index = find_supp(index,n);
    supp_zeros = find(dec2bin(bitxor(2^n-1,index),n)-'0'==1);  %complement of supp
    
    %Sm_Si = set(supp_index)
    Sm_Ti = supp_zeros;
    
    %for s = 2:length(J)
    row_sets = J; %nchoosek(J, s);
    %for row_idx = 1:size(row_sets,1)
    valid_set_for_M = true;
    Sm_S = supp_index;
    Sm_T = [];
    %column_elmnt_cnt = zeros(length(supp_zeros),1);
    for r = row_sets %(row_idx,:)
        %r_bin = dec2bin(r,n);
        supp_r = find_supp(r,n);
        Sm_S = intersect(Sm_S, supp_r);
        if ~isempty(intersect(Sm_T,intersect(Sm_Ti, supp_r))) %distinctness condition
            valid_set_for_M = false; % No need to set M.
            break;
        end
        Sm_T = union(Sm_T,intersect(Sm_Ti, supp_r));
        %{
        cnt=0;
        for t = supp_zeros % The difference with the upper if block? 
            if r_bin(t)=='1'
                column_elmnt_cnt(cnt) = column_elmnt_cnt(cnt) + 1;
                if column_elmnt_cnt(cnt) > 1    % No memebr of J would satisfy this condition.
                    valid_set_for_M = false; % No need to set M.
                    break;
                end
            end
            cnt = cnt + 1;
        end
        %}
        %if valid_set_for_M == false
            %break;
        %end
    end
    if valid_set_for_M == true
        Sm = union(Sm_T,Sm_S);
        m = supp2dec(Sm,n);
    end
    %end
    %end
    % Removeing the repeated elements by even times.
    %{
    M1 = unique(M); %It should be done with the previous summation on level elements and this Sm.
    count_M = histc(M,M1);
    indices = [];
    for cnt = 1:length(count_M)
        if mod(count_M(cnt),2) == 0
            indices = [indices cnt];
        end
    end
    M1(indices) = [];
    M = M1;
    %}
end


function m1 = find_set_M1R(index,m,j,n)
    %M = [];
    m1 = 0;
    supp_index = find_supp(index,n);
    supp_zeros = find(dec2bin(bitxor(2^n-1,index),n)-'0'==1);  %complement of supp
    
    Ti = supp_zeros;
    Si = supp_index;

    supp_j = find_supp(j,n);

    supp_m = find_supp(m,n);
    Sm_S = intersect(intersect(Si, supp_j), supp_m);
    Tj = intersect(Ti, supp_j);
    if isempty(intersect(intersect(Ti, supp_m), Tj)) %distinctness condition
        Sm_T = union(Tj,intersect(Ti, supp_m));
        Sm = union(Sm_T,Sm_S);
        m1 = supp2dec(Sm,n);
    end
end



function supp = find_supp(i,n)
    supp = find(dec2bin(i,n)-'0'==1);
end

function decimal = supp2dec(indices,n)
    decimal = 0;
    for i = indices 
        decimal = decimal + 2^(n-i);
    end
end

function decimal = get_idx_MsubKiRed(row_indices,KiRed)
    %n = size(KiRed,1);
    decimal = 0;
    for i = row_indices 
        decimal = decimal + 2^(find(KiRed==i)-1);
    end
    decimal = decimal + 1;
end


%%

% retruns Mean-LLRs (a measure for sub-channels' reliability) obtained from Density Evolution by Gaussian Approximation (DEGA) 
function I =  construct_dega(design_snr_db, N, K)
    mllr = zeros(N,1);
    sigma_sq = 1/(2*K/N*power(10,design_snr_db/10));
    mllr(1) = 2/sigma_sq;
    for level = 1:log2(N)
        B = 2^level;
        for j = 1:B / 2
            T = mllr(j);
            mllr(j) = calc_phi_inv(T);
            mllr(B / 2 + j) = 2 * T;
        end
    end
    
    mask = zeros(N,3);
    nat = bitrevorder(0:N-1);
    for i = 1:N
        mask(i,:) = [nat(i), mllr(i), 1];
    end
    % sort sub-channels by mllr
    mask = sortrows(mask,2); %direction: ascend (default)
    % set info bits to 1 for sub-channels with K largest mllr values
    for i = 1:N-K
        mask(i,3) = 0;
    end
    % sort channels with respect to index (in bitreversal order; line 42
    mask = sortrows(mask,1); %direction: ascend (default)
    I = find(mask(:,3)==1)-1;
end

% returns Phi inverse based on piece-wise linear approximation
function phi_inv = calc_phi_inv(x)
    if (x>12)
        phi_inv = 0.9861 * x - 2.3152;
    elseif (x<=12 && x>3.5)
        phi_inv = x*(0.009005 * x + 0.7694) - 0.9507;
    elseif (x<=3.5 && x>1)
        phi_inv = x*(0.062883*x + 0.3678)- 0.1627;
    else
        phi_inv = x*(0.2202*x + 0.06448);
    end
end




function u = conv_1bit(in_bit,cur_state,g)
    u = in_bit * g(1); % by convention, we have always g(1)=1
    for i = 2:length(g)
        if g(i)==1
            u = bitxor(u,cur_state(i-1));
        end
    end
end

function [u,state] = conv_1bit_nxtState(in_bit,cur_state,g)
    u = in_bit * g(1); % by convention, we have always g(1)=1
    for i = 2:length(g)
        if g(i)==1
            u = bitxor(u,cur_state(i-1));
        end
    end
    m = length(cur_state);
    if in_bit == 0
        state = [0,cur_state(1:m-1)];
    else
        state = [1,cur_state(1:m-1)];
    end        
end
