N = 2^7;            % Code length
R = 0.8594;            % Code rate
K = N * R;          % code dimention: the number of information bits
design_snr_db = 4;  % for the optimized code
g = [1 0 1 1 0 1 1];%[1 0 1 1];% 0 1 1 0 1 1];%[1 0 0 0 0 1 1];%[1 0 1 0 1 0 1];%
I = construct_dega(design_snr_db, N, K); % the indices of K information bits

tic
[d_min,A_dminPolar,A_dminPAC] = err_coeff(I,N,g) 
toc

%% The main function

% returns the minimum distance of the code and the number of codewords with minumum distance (error coefficient)
function [dmin, A_dmin, A_dminPAC] = err_coeff(I,N,g)
    d = min(sum(dec2bin(I)-'0',2));
    dmin = 2^d; n = log2(N); A_dmin = 0; A_dminPAC = 0; A_dminPAC2 = 0;
    B = find(sum(dec2bin(I)-'0',2)==d); 
    set_F = find_F(I); % F indicator
    cosetsAdmin = zeros(length(B),4);
    cosetsAdmin(:,1) = I(B); 
    for i = B' % a set of indices of elements in I, not elements of I
        Ki_size = n - d;
        for x = find(dec2bin(I(i),n)-'0'==1) % support of I(i): positions of 1's in the binary representation of the sub-channel index
                ii = dec2bin(bitxor(N-1,I(i)),n)-'0';
                Ki_size = Ki_size + sum(ii(1:x-1));
        end
        
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

        A_dmin_minus = 0;
        A_Kif = 0;
        A_dmin_minus_f = 0;
        for x = 0:length(set_KiRed)
            A = nchoosek(set_KiRed,x); %each row, one combination
            if x==0
                A=[0];
            end
            for y = 1:size(A,1) %parfor
                J = [];
                v = zeros(N,1); u = zeros(N,1);
                cur_state = zeros(1,length(g)-1);
                v(i_lead+1) = 1; 
                [u(i_lead+1),cur_state] = conv_1bit_nxtState(v(i_lead+1),cur_state,g); %the i-th bit first, +1 is due to starting index in MATLAB.
                ones = zeros(1,max_f-i_lead);
                M = [];
                cnt_i = 1; cnt_f = 1;
                for j = i_lead+1:max_f
                    if set_F(j+1) == 0
                        if A(y,cnt_i)==j
                            [M, ones1] = find_ones_1R(i_lead,max_f,J,M,j,n);
                            ones = bitxor(ones, ones1);
                            J = [J j];
                            if conv_1bit(v(j+1),cur_state,g) == 0
                                v(j+1) = 1;
                            end
                            if cnt_i<x
                                cnt_i = cnt_i + 1;
                            end
                        else
                            if ones(j-i_lead)==1
                                if conv_1bit(v(j+1),cur_state,g) == 0
                                    v(j+1) = 1;
                                end
                            else
                                if conv_1bit(v(j+1),cur_state,g) == 1
                                    v(j+1) = 1;
                                end

                            end
                        end
                        [u(j+1),cur_state] = conv_1bit_nxtState(v(j+1),cur_state,g);
                    else
                        [u(j+1),cur_state] = conv_1bit_nxtState(v(j+1),cur_state,g);
                        if set_Fxi(cnt_f)==j
                            if u(j+1)==1
                                if ones(j-i_lead)==0
                                    A_dmin_minus = A_dmin_minus + 2^(Ki_size-length(set_KiRed));
                                    A_Kif = A_Kif + 1;
                                    break;
                                end
                            else
                                if ones(j-i_lead)==1
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
                                [M, ones1] = find_ones_1R(i_lead,max_f,J,M,j,n);
                                ones = bitxor(ones, ones1);
                                J = [J j];
                            end
                        end
                    end
                end
            end
        end
        A_dminPAC = A_dminPAC - (A_dmin_minus + A_dmin_minus_f);
        cosetsAdmin(find(B==i),3) = A_dmin_minus + A_dmin_minus_f;

    end
    cosetsAdmin(:,4) = cosetsAdmin(:,2) - cosetsAdmin(:,3);
    disp('          i  A_dmin(Polar)  Reduction  A_dmin(PAC)')
    cosetsAdmin
end



%% Fuctions used to form auxiliary sets

function set_Fxi = find_Fxi(index,F,n) 
    i = index + 2;  % Knowing MATLAB does not support index 0 for arrays and starting from i+1
    set_Fxi = [];
    while i < 2^n
        if F(i) == 1 && sum(dec2bin(bitxor(index,bitor(i-1,index)),n)-'0'==1)>1
            set_Fxi = [set_Fxi, i-1];
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


%% Functions related to formation of set M

function [M, ones] = find_ones_1R(i_lead,max_f,J,M,j,n) %Progressively
    ones = zeros(1,max_f-i_lead);
    for x=1:length(M)
        %Mj = union(M(x),j);
        m  = find_set_M1R(i_lead,M(x),j,n);
        if m>0 && m <= max_f
            M = [M m]; %Make sure M will not have a repetitive elements. It will be cancelled in the next line anyways.
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
    m = 0;
    supp_index = find_supp(index,n);
    supp_zeros = find(dec2bin(bitxor(2^n-1,index),n)-'0'==1);  %complement of supp
    
    Sm_Ti = supp_zeros;
    
    row_sets = J; %nchoosek(J, s);
    valid_set_for_M = true;
    Sm_S = supp_index;
    Sm_T = [];
    for r = row_sets %(row_idx,:)
        supp_r = find_supp(r,n);
        Sm_S = intersect(Sm_S, supp_r);
        if ~isempty(intersect(Sm_T,intersect(Sm_Ti, supp_r))) %distinctness condition
            valid_set_for_M = false; % No need to set M.
            break;
        end
        Sm_T = union(Sm_T,intersect(Sm_Ti, supp_r));
    end
    if valid_set_for_M == true
        Sm = union(Sm_T,Sm_S);
        m = supp2dec(Sm,n);
    end
end


function m1 = find_set_M1R(index,m,j,n)
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


%% Polar COde Construction: Density Evolution based on Guassian Approximation

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
    for i = 0:N-1
        nat(i+1) = bitreversed(i,uint8(log2(N)));
    end
    %nat = bitrevorder(0:N-1);
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

function dec = bitreversed(i,n) % bitrevorder() is in singal processing toolbox.
    dec = bin2dec(fliplr(dec2bin(i,n)));
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
