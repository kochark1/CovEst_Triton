load(strcat(interm_folder,'/RQWfile.mat'));

ch_samples = 2000;

dim = 2;

SII = eye(number_of_antennas)*(1/mu_val);
ZZ_2= zeros(number_of_antennas, ch_samples);
ZZ_3= zeros(ch_samples,1);
ZZ_4= zeros(number_of_antennas, number_of_antennas, ch_samples);

h_mat = generateh(R_sqrt_root,number_of_antennas,number_of_cells,number_of_users,ch_samples); %  (M, L, K, ch_samples)


num_hat_t = 0;
num_hat_diag_t = 0;
num_hat_diag_reg_t = 0;

num_hat_DL_t = 0;
num_hat_diag_DL_t = 0;
num_hat_diag_reg_DL_t = 0;

den1_hat_t = 0;
den1_hat_diag_t = 0;
den1_hat_diag_reg_t = 0;

den1_hat_DL_t = 0;
den1_hat_diag_DL_t = 0;
den1_hat_diag_reg_DL_t = 0;

den2_hat_t = 0;
den2_hat_diag_t = 0;
den2_hat_diag_reg_t = 0;

for ii = 1:ch_samples
    sum_mat=0;
    sum_mat_diag=0;
    sum_mat_diag_reg=0;
    for k= 1:number_of_users
        [s1, h1] = mmse_cal(Wmatrix, number_of_antennas, dim, pilotSequenceLength, mu_val, h_mat(:,:,k,ii));
        [s2, h2] = mmse_cal(W_diag, number_of_antennas, dim, pilotSequenceLength, mu_val, h_mat(:,:,k,ii));
        
        if k== targetUser
            h_target = h1;
            h_target_diag = h2;
        end
        sum_mat = sum_mat + s1;
        sum_mat_diag = sum_mat_diag + s2;
    end

    v = compute_v_or_b(sum_mat, h_target, SII, 0); % M, 1
    v_diag = compute_v_or_b(sum_mat_diag, h_target_diag, SII, 0); % M, 1
    
    b = compute_v_or_b(sum_mat, h_target, SII, 1);
    b_diag = compute_v_or_b(sum_mat_diag, h_target_diag, SII, 1);
    
    num_hat_t = num_hat_t + (1/ch_samples) * (v'*h_mat(:,targetCell,targetUser, ii));
    num_hat_diag_t = num_hat_diag_t + (1/ch_samples) * (v_diag'*h_mat(:,targetCell,targetUser, ii));

    num_hat_DL_t = num_hat_DL_t + (1/ch_samples) * (b'*h_mat(:,targetCell,targetUser, ii));
    num_hat_diag_DL_t = num_hat_diag_DL_t + (1/ch_samples) * (b_diag'*h_mat(:,targetCell,targetUser, ii));

    hh = zeros(number_of_antennas, number_of_antennas);

    for l=1:number_of_cells
        for k=1:number_of_users
            hh = hh + (h_mat(:,l,k, ii)*h_mat(:,l,k, ii)');
        end
    end

    den1_hat_t = den1_hat_t + (1/ch_samples) * trace(hh * (v*v'));
    den1_hat_diag_t = den1_hat_diag_t + (1/ch_samples) * trace(hh * (v_diag*v_diag'));

    den1_hat_DL_t = den1_hat_DL_t + (1/ch_samples) * trace(hh * (b*b'));
    den1_hat_diag_DL_t = den1_hat_diag_DL_t + (1/ch_samples) * trace(hh * (b_diag*b_diag'));

    den2_hat_t = den2_hat_t + (1/ch_samples) * (1/mu_val) * (v'*v);
    den2_hat_diag_t = den2_hat_diag_t + (1/ch_samples) * (1/mu_val) * (v_diag'*v_diag);
end

num_hat = abs(num_hat_t)^2;
num_hat_diag = abs(num_hat_diag_t)^2;

den_hat = (den1_hat_t + den2_hat_t) - num_hat;
den_hat_diag = (den1_hat_diag_t + den2_hat_diag_t) - num_hat_diag;

gamma_hat = real(num_hat/den_hat);
gamma_hat_diag = real(num_hat_diag/den_hat_diag);

SE_const_ZF = log2(1+gamma_hat);
SE_const_diag_ZF = log2(1+gamma_hat_diag);
SE_const_diag_reg_ZF = SE_const_diag_ZF;

num_hat_DL = abs(num_hat_DL_t)^2;
num_hat_diag_DL = abs(num_hat_diag_DL_t)^2;

den_hat_DL = (den1_hat_DL_t) - num_hat_DL  +(1/lambda_DL);
den_hat_diag_DL = (den1_hat_diag_DL_t) - num_hat_diag_DL +(1/lambda_DL);

gamma_hat_DL = real(num_hat_DL/den_hat_DL);
gamma_hat_diag_DL = real(num_hat_diag_DL/den_hat_diag_DL);

SE_const_DL_ZF  = log2(1+gamma_hat_DL);
SE_const_diag_DL_ZF  = log2(1+gamma_hat_diag_DL);
SE_const_diag_reg_DL_ZF  = SE_const_diag_DL_ZF;


function v_or_b = compute_v_or_b(pre_mat, h_target, SII, ULDL_flag)
    v = (pre_mat+SII)\h_target;
%     v = h_target;
    const = 1; % TODO
    if ULDL_flag == 0
        v_or_b = v;
    else
        v_or_b = v/const;
    end
end

function [sum_mat, h_out]= mmse_cal(W, number_of_antennas, dim, P, mu, h_mat)

    h_LS = sum(h_mat, dim) + (1/sqrt(P*mu)) * crandn(number_of_antennas,1);

    h_mmse = W*h_LS;
    sum_mat = h_mmse*h_mmse'; % sum over k
    h_out = h_mmse;

end