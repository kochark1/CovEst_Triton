tic;
load(strcat(interm_folder,'/RQWfile.mat'));

RR = Rmatrices(:,:,targetCell,targetUser);
SS = diag(diag(RR));
PP = diag(diag(Qmatrix));
Ssum = diag(diag(Rsum));
Ssum_DL = diag(diag(Rsum_DL));

% Terms in Theorem 1 without the multiplicative factors involving NR & NQ
% Details given in the function
[Enum_const, Eden1_const, Eden2_const] = preComputExpectation(Rmatrices, Qmatrix, Wmatrix, W_bar, Rsum, number_of_cells, number_of_antennas, targetCell, targetUser, alpha_R);
[Enum_const_DL, Eden1_const_DL, Eden2_const_DL] = preComputExpectation(Rmatrices, Qmatrix, Wmatrix, W_bar, Rsum_DL, number_of_cells, number_of_antennas, targetCell, targetUser, alpha_R);

% Terms in Theorem 2 without the multiplicative factors involving NR & NQ
% Details given in the function
[Enum_const_diag, Eden1_const_diag, Eden2_const_diag] = preComputExpectation_daig(Rmatrices, Qmatrix, W_diag, W_bar_diag, Rsum, number_of_cells, targetCell, targetUser, alpha_R);
[Enum_const_diag_DL, Eden1_const_diag_DL, Eden2_const_diag_DL] = preComputExpectation_daig(Rmatrices, Qmatrix, W_diag, W_bar_diag, Rsum_DL, number_of_cells, targetCell, targetUser, alpha_R);

prelogFactor = (1 - pilotSequenceLength/ulCoherenceSymbols - NR_vec*pilotSequenceLength/(ulCoherenceSymbols*fullFramecoherenceBlocks));

% Terms in Theorem 3 are generated in the following for loops since they
% are not in closed-form

for nr = 1:length(NR_vec)
    for nq = 1:length(NQ_vec)
        NQ = NQ_vec(nq);
        NR = NR_vec(nr);

        %--- Theoretical SE for LMMSE
        kappa1 = NQ^3/((NQ-number_of_antennas)^3-(NQ-number_of_antennas)); % defined in Theorem 1
        kappa2 = NQ^2/((NQ-number_of_antennas)^2-1); % defined in Theorem 1

        num_sc = (NQ/(NQ-number_of_antennas))^2; % multiplicative factor(s) to match values returned by preComputExpectation()
        den_v1 = [kappa1; kappa1/NR]; % multiplicative factor(s) to match values returned by preComputExpectation()
        den_v2 = [kappa2; kappa2/NR; kappa1/NQ; kappa1/(NQ*NR)]; % multiplicative factor(s) to match values returned by preComputExpectation()

        gamma_theo = real((num_sc*Enum_const)/(den_v1'*Eden1_const + den_v2'*Eden2_const - num_sc*Enum_const)); % Eq 14 with LMMSE
        gamma_theo_DL = real((num_sc*Enum_const_DL)/(den_v1'*Eden1_const_DL + den_v2'*Eden2_const_DL - num_sc*Enum_const_DL +(1/lambda_DL))); % Eq 15 with LMMSE
        
        SE_theo(nr,nq) = prelogFactor(nr)*log2(1 + gamma_theo);
        SE_theo_DL(nr,nq) = log2(1 + gamma_theo_DL);
        
        if nr==1 && nq == nq_set(2)
            a = num_sc*Enum_const;
            c = [0 kappa1]*Eden1_const + [0 kappa2 0 kappa1/NQ]*Eden2_const;
            d = 0;
            b = [kappa1 0]*Eden1_const + [kappa2 0 kappa1/NQ 0]*Eden2_const - a + d;

            a_DL = num_sc*Enum_const_DL;
            c_DL = [0 kappa1]*Eden1_const_DL + [0 kappa2 0 kappa1/NQ]*Eden2_const_DL;
            d_DL = 1/lambda_DL;
            b_DL = [kappa1 0]*Eden1_const_DL + [kappa2 0 kappa1/NQ 0]*Eden2_const_DL - a_DL + d_DL;
        end
        %--- Theoretical SE for elementwise LMMSE
        kappa3 = (NQ/(NQ-1))^2; % defined in Theorem 2
        kappa4 = kappa3/(NQ-2); % defined in Theorem 2

        num_sc = kappa3; % multiplicative factor(s) to match values returned by preComputExpectation_diag
        den_v1 = [kappa3; kappa3/NR; kappa4; kappa4/NR]; % multiplicative factor(s) to match values returned by preComputExpectation_diag
        den_v2 = [kappa3; kappa3/NR; kappa4; kappa4/NR]; % multiplicative factor(s) to match values returned by preComputExpectation_diag

        gamma_theo = real((num_sc*Enum_const_diag)/(den_v1'*Eden1_const_diag + den_v2'*Eden2_const_diag - num_sc*Enum_const_diag)); % Eq 14 with elementwise LMMSE
        gamma_theo_DL = real((num_sc*Enum_const_diag_DL)/(den_v1'*Eden1_const_diag_DL + den_v2'*Eden2_const_diag_DL - num_sc*Enum_const_diag_DL+(1/lambda_DL))); % Eq 14 with elementwise LMMSE
        
        SE_theo_diag(nr, nq) = prelogFactor(nr)*log2(1 + gamma_theo);
        SE_theo_diag_DL(nr, nq) = log2(1 + gamma_theo_DL);

        if nr==1 && nq == nq_set(2)
            f = num_sc*Enum_const_diag;
            h = [0 kappa3 0 kappa4]*Eden1_const_diag + [0 kappa3 0 kappa4]*Eden2_const_diag;
            d = 0;
            g = [kappa3 0 kappa4 0]*Eden1_const_diag + [kappa3 0 kappa4 0]*Eden2_const_diag - f + d;

            f_DL = num_sc*Enum_const_diag_DL;
            h_DL = [0 kappa3 0 kappa4]*Eden1_const_diag_DL + [0 kappa3 0 kappa4]*Eden2_const_diag_DL;
            d_DL = 1/lambda_DL;
            g_DL = [kappa3 0 kappa4 0]*Eden1_const_diag_DL + [kappa3 0 kappa4 0]*Eden2_const_diag_DL - f_DL + d_DL;
            
            N_thr_UL = real((f*c-a*h)/(a*g-f*b));
            N_thr_DL = real((f_DL*c_DL-a_DL*h_DL)/(a_DL*g_DL-f_DL*b_DL));  
%             NQ_vec(nq)
%             [N_thr_UL N_thr_DL]
        end
        %--- Theoretical SE for elementwise LMMSE with regularized P
        E = E_vec(:,:,nq);
        G = G_vec(:,:,nq);
        EQE = E*Qmatrix*E;
        GmE2 = G-E^2;
        Enum_const_diag_reg = (abs(trace(E*SS^2)))^2;
        Eden1_const_diag_reg = trace(EQE*SS*Rsum*SS) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum.*Qmatrix.*Qmatrix)) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum.*RR.*RR)) ...
            + (1 + alpha_R^2/(2*NR))*trace(GmE2*PP*SS^2*Ssum) + (alpha_R^2/(2*NR))*trace(GmE2*PP^3*Ssum);
        Eden1_const_diag_reg_DL = trace(EQE*SS*Rsum_DL*SS) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum_DL.*Qmatrix.*Qmatrix)) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum_DL.*RR.*RR)) ...
            + (1 + alpha_R^2/(2*NR))*trace(GmE2*PP*SS^2*Ssum_DL) + (alpha_R^2/(2*NR))*trace(GmE2*PP^3*Ssum_DL);

        Eden2_const_diag_reg = 0;
        for l =1:number_of_cells
            Sjlu = diag(diag(Rmatrices(:,:,l,targetUser)));
            Eden2_const_diag_reg = Eden2_const_diag_reg + abs(trace(E*SS*Sjlu))^2+ (alpha_R^2/(2*NR))*sum(sum(E*Sjlu*(Qmatrix.*Qmatrix)*E*Sjlu))...
                + (alpha_R^2/(2*NR))*sum(sum(E*Sjlu*(RR.*RR)*E*Sjlu)) + (1 + alpha_R^2/(2*NR))*trace(GmE2*SS^2*Sjlu^2)...
                + (alpha_R^2/(2*NR))*trace(GmE2*PP^2*Sjlu^2);
        end

        gamma_theo = real((Enum_const_diag_reg)/(Eden1_const_diag_reg + Eden2_const_diag_reg - Enum_const_diag_reg)); % Eq 14 with elementwise LMMSE + regularized P
        gamma_theo_DL = real((Enum_const_diag_reg)/(Eden1_const_diag_reg_DL + Eden2_const_diag_reg - Enum_const_diag_reg +(1/lambda_DL))); % Eq 15 with elementwise LMMSE + regularized P
        
        SE_theo_diag_reg(nr, nq) = prelogFactor(nr)*log2(1 + gamma_theo);
        SE_theo_diag_reg_DL(nr, nq) = log2(1 + gamma_theo_DL);
        
        [nr nq toc]
    end
end
if ZF_FLAG
    true_ZF
end

save(strcat(out_folder,'/true_SE.mat'), 'SE_theo', 'SE_theo_DL', 'SE_theo_diag', 'SE_theo_diag_DL',...
    'N_thr_UL', 'N_thr_DL', 'SE_theo_diag_reg', 'SE_theo_diag_reg_DL', 'SE_const_ZF', 'SE_const_diag_ZF',...
    'SE_const_diag_reg_ZF','SE_const_DL_ZF','SE_const_diag_DL_ZF','SE_const_diag_reg_DL_ZF', 'prelogFactor');
disp('True and theoretical Done!')
% SE_theo(nr,nq)
% SE_theo_DL(nr,nq)
% 
% SE_theo_diag(nr, nq)
% SE_theo_diag_DL(nr, nq)
% 
% SE_theo_diag_reg(nr, nq)
% SE_theo_diag_reg_DL(nr, nq)