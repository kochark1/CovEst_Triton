load 'matfiles/RQWfile.mat';


[Enum_const, Eden1_const, Eden2_const] = preComputExpectation(Rmatrices, Qmatrix, Wmatrix, W_bar, Rsum, number_of_cells, number_of_antennas, targetCell, targetUser, alpha_R);
[Enum_const_DL, Eden1_const_DL, Eden2_const_DL] = preComputExpectation(Rmatrices, Qmatrix, Wmatrix, W_bar, Rsum_DL, number_of_cells, number_of_antennas, targetCell, targetUser, alpha_R);
prelogFactor = (1 - pilotSequenceLength/ulCoherenceSymbols - NR_vec*pilotSequenceLength/(ulCoherenceSymbols*fullFramecoherenceBlocks));

for nr = 1:length(NR_vec)
    for nq = nq_sim:nq_sim %1:length(NQ_vec)
        clear W_Est;
        NQ = NQ_vec(nq);
        NR = NR_vec(nr);

        load(strcat('matfiles/RQest_',string(nr), '_', string(nq),'_', string(targetUser), '.mat'), 'W_Est');

        num_hat_t = zeros(TT,1);
        den1_hat_t = zeros(TT,1);
        den1_hat_DL_t = zeros(TT,1);
        den2_hat_t = zeros(TT,1);

        for t = 1:TT
            W_est = W_Est(:,:,t);

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
        SE_hat(nr,nq) = prelogFactor(nr)*log2(1+gamma_hat);

        den_hat_DL = mean(den1_hat_DL_t + den2_hat_t) - num_hat;
        gamma_hat_DL = real(num_hat/den_hat_DL);
        SE_hat_DL(nr,nq) = log2(1+gamma_hat_DL);

        %--- Expected SE Theoretical
        kappa1 = NQ^3/((NQ-number_of_antennas)^3-(NQ-number_of_antennas));
        kappa2 = NQ^2/((NQ-number_of_antennas)^2-1);

        num_sc = (NQ/(NQ-number_of_antennas))^2;
        den_v1 = [kappa1; kappa1/NR];
        den_v2 = [kappa2; kappa2/NR; kappa1/NQ; kappa1/(NQ*NR)];

        gamma_theo = real((num_sc*Enum_const)/(den_v1'*Eden1_const + den_v2'*Eden2_const - num_sc*Enum_const));
        gamma_theo_DL = real((num_sc*Enum_const_DL)/(den_v1'*Eden1_const_DL + den_v2'*Eden2_const_DL - num_sc*Enum_const_DL +(1/lambda_DL)));
        SE_theo(nr,nq) = prelogFactor(nr)*log2(1 + gamma_theo);
        SE_theo_DL(nr,nq) = log2(1 + gamma_theo_DL);

        num_verify1(nr,nq) = real(num_hat);
        num_verify2(nr,nq) = real(num_sc*Enum_const);

        den1_verify1(nr,nq) =  real(mean(den1_hat_t));
        den1_verify2(nr,nq) = real(den_v1'*Eden1_const);

        den2_verify1(nr,nq) =  real(mean(den2_hat_t));
        den2_verify2(nr,nq) = real(den_v2'*Eden2_const);

        [nr nq toc]
    end
end
% figure, surf(NQ_vec, NR_vec, (1 - number_of_users/ulCoherenceSymbols)*SE_const*ones(size(SE_theo)));
% hold on;
% surf(NQ_vec, NR_vec, SE_theo);
% hold off;
% xlabel('N_{Q}');
% ylabel('N_{R}');
% zlabel('UL Spectral efficiency per user [bits/s/Hz]');
% legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Location','NorthEast');
% 
% figure, surf(NQ_vec, NR_vec, SE_const_DL*ones(size(SE_theo)));
% hold on;
% surf(NQ_vec, NR_vec, SE_theo_DL);
% hold off;
% xlabel('N_{Q}');
% ylabel('N_{R}');
% zlabel('DL Spectral efficiency per user [bits/s/Hz]');
% legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Location','NorthEast');

% % % figure, plot(NR_vec, (1 - number_of_users/ulCoherenceSymbols)*SE_const*ones(size(SE_theo(:,nq))),'^-');
figure, plot(NR_vec, prelogFactor*SE_const,'^-');
hold on;
plot(NR_vec, SE_theo(:,nq_sim),'o-');
plot(NR_vec, SE_hat(:,nq_sim),'g*-');
hold off;

xlabel('N_{R}');
ylabel('UL Spectral efficiency per user [bits/s/Hz]');
legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Spectral efficiency (Numerical) - estimated covariance','Location','SouthEast');
% 
% 
% figure, plot(NR_vec, SE_const_DL*ones(size(SE_theo_DL(:,nq))),'^-');
% hold on;
% plot(NR_vec, SE_theo_DL(:,nq),'o-');
% plot(NR_vec, SE_hat_DL(:,nq),'g*-');
% hold off;

% xlabel('N_{R}');
% ylabel('DL Spectral efficiency per user [bits/s/Hz]');
% legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Spectral efficiency (Numerical) - estimated covariance','Location','SouthEast');

figure, plot(NR_vec, num_verify1(:,nq_sim),'o-');
ylim([0, max(num_verify1(:,nq_sim))*1.1]);
hold on;
plot(NR_vec, num_verify2(:,nq_sim),'*-');
hold off;

figure, plot(NR_vec, den1_verify1(:,nq_sim),'o-');
ylim([0, max(den1_verify1(:,nq_sim))*1.1]);
hold on;
plot(NR_vec, den1_verify2(:,nq_sim),'*-');
hold off;

figure, plot(NR_vec, den2_verify1(:,nq_sim),'o-');
ylim([0, max(den2_verify1(:,nq_sim))*1.1]);
hold on;
plot(NR_vec, den2_verify2(:,nq_sim),'*-');
hold off;
