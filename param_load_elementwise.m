load 'matfiles/RQWfile.mat';
[Enum_const_diag, Eden1_const_diag, Eden2_const_diag] = preComputExpectation_daig(Rmatrices, Qmatrix, W_diag, W_bar_diag, Rsum, number_of_cells, targetCell, targetUser, alpha_R);
[Enum_const_diag_DL, Eden1_const_diag_DL, Eden2_const_diag_DL] = preComputExpectation_daig(Rmatrices, Qmatrix, W_diag, W_bar_diag, Rsum_DL, number_of_cells, targetCell, targetUser, alpha_R);
prelogFactor = (1 - pilotSequenceLength/ulCoherenceSymbols - NR_vec*number_of_users/(ulCoherenceSymbols*fullFramecoherenceBlocks));

for nr = 1:length(NR_vec)
    for nq = nq_sim:nq_sim
        clear W_Est_diag;
        NQ = NQ_vec(nq);
        NR = NR_vec(nr);

        load(strcat('matfiles/RQest_',string(nr), '_', string(nq),'_', string(targetUser), '.mat'), 'W_Est_diag');

        num_hat_t = zeros(TT,1);
        den1_hat_t = zeros(TT,1);
        den1_hat_DL_t = zeros(TT,1);
        den2_hat_t = zeros(TT,1);

        for t = 1:TT
            W_est = W_Est_diag(:,:,t);
            num_hat_t(t) = trace(W_est'*Rmatrices(:,:,targetCell, targetUser));
            den1_hat_t(t) = trace(W_est*Qmatrix*W_est'*Rsum);
            den1_hat_DL_t(t) = trace(W_est*Qmatrix*W_est'*Rsum_DL);
            for l=1:number_of_cells
                den2_hat_t(t) = den2_hat_t(t) + abs(trace(W_est'*Rmatrices(:,:,l, targetUser)))^2;
            end
        end

        num_hat = abs(mean(num_hat_t))^2;
        den_hat = mean(den1_hat_t + den2_hat_t) - num_hat;
        gamma_hat = real(num_hat/den_hat);
        SE_hat(nr, nq) = prelogFactor(nr)*log2(1+gamma_hat);

        den_hat_DL = mean(den1_hat_DL_t + den2_hat_t) - num_hat;
        gamma_hat_DL = real(num_hat/den_hat_DL);
        SE_hat_DL(nr, nq) = log2(1+gamma_hat_DL);
        
        %--- Expected SE Theoretical
        kappa3 = (NQ/(NQ-1))^2;
        kappa4 = kappa3/(NQ-2);

        num_sc = kappa3;
        den_v1 = [kappa3; kappa3/NR; kappa4; kappa4/NR];
        den_v2 = [kappa3; kappa3/NR; kappa4; kappa4/NR];

        gamma_theo = real((num_sc*Enum_const_diag)/(den_v1'*Eden1_const_diag + den_v2'*Eden2_const_diag - num_sc*Enum_const_diag));
        SE_theo(nr, nq) = prelogFactor(nr)*log2(1 + gamma_theo);

        gamma_theo_DL = real((num_sc*Enum_const_diag_DL)/(den_v1'*Eden1_const_diag_DL + den_v2'*Eden2_const_diag_DL - num_sc*Enum_const_diag_DL+(1/lambda_DL)));
        SE_theo_DL(nr, nq) = log2(1 + gamma_theo_DL);

        num_verify1(nr, nq) = real(num_hat);
        num_verify2(nr, nq) = real(num_sc*Enum_const_diag);

        den1_verify1(nr, nq) =  real(mean(den1_hat_t));
        den1_verify2(nr, nq) = real(den_v1'*Eden1_const_diag);

        den2_verify1(nr, nq) =  real(mean(den2_hat_t));
        den2_verify2(nr, nq) = real(den_v2'*Eden2_const_diag);

        [nr nq toc]
    end
end

figure, plot(NR_vec, (1 - number_of_users/ulCoherenceSymbols)*SE_const_diag*ones(size(SE_theo(:,nq_sim))),'^-');
% ylim([0, SE_const_diag*1.1]);
hold on;
plot(NR_vec, SE_theo(:,nq_sim),'o-');
plot(NR_vec, SE_hat(:,nq_sim),'g*-');
hold off;

xlabel('N_{R}');
ylabel('UL Spectral efficiency per user [bits/s/Hz]');
legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Spectral efficiency (Numerical) - estimated covariance','Location','SouthEast');

figure, plot(NR_vec, SE_const_diag_DL*ones(size(SE_theo_DL(:,nq_sim))),'^-');
% ylim([0, SE_const_diag_DL*1.1]);
hold on;
plot(NR_vec, SE_theo_DL(:,nq_sim),'o-');
plot(NR_vec, SE_hat_DL(:,nq_sim),'g*-');
hold off;

xlabel('N_{R}');
ylabel('DL Spectral efficiency per user [bits/s/Hz]');
legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Spectral efficiency (Numerical) - estimated covariance','Location','SouthEast');

% -----------------------------------------------
figure, plot(NR_vec, num_verify1(:,nq_sim),'o-');
ylim([0, max(num_verify1(:,nq_sim))*1.1]);
hold on;
plot(NR_vec, num_verify2(:,nq_sim),'*-');
hold off;
norm(num_verify1(:,nq_sim)-num_verify2(:,nq_sim))/norm(num_verify2(:,nq_sim))*100

figure, plot(NR_vec, den1_verify1(:,nq_sim),'o-');
ylim([0, max(den1_verify1(:,nq_sim))*1.1]);
hold on;
plot(NR_vec, den1_verify2(:,nq_sim),'*-');
hold off;
norm(den1_verify1(:,nq_sim)-den1_verify2(:,nq_sim))/norm(den1_verify2(:,nq_sim))*100

figure, plot(NR_vec, den2_verify1(:,nq_sim),'o-');
ylim([0, max(den2_verify1(:,nq_sim))*1.1]);
hold on;
plot(NR_vec, den2_verify2(:,nq_sim),'*-');
hold off;
norm(den2_verify1(:,nq_sim)-den2_verify2(:,nq_sim))/norm(den2_verify2(:,nq_sim))*100