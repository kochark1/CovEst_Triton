load 'matfiles/RQWfile.mat';
RR = Rmatrices(:,:,targetCell,targetUser);
SS = diag(diag(RR));
PP = diag(diag(Qmatrix));
Ssum = diag(diag(Rsum));
Ssum_DL = diag(diag(Rsum_DL));
prelogFactor = (1 - pilotSequenceLength/ulCoherenceSymbols - NR_vec*number_of_users/(ulCoherenceSymbols*fullFramecoherenceBlocks));

for nr = 1:length(NR_vec)
    for nq = nq_sim:nq_sim
        clear W_Est_diag_reg;
        NR = NR_vec(nr);

        load(strcat('matfiles/RQest_',string(nr), '_', string(nq),'_', string(targetUser), '.mat'), 'W_Est_diag_reg');

        num_hat_t = zeros(TT,1);
        den1_hat_t = zeros(TT,1);
        den1_hat_DL_t = zeros(TT,1);
        den2_hat_t = zeros(TT,1);

        for t = 1:TT
            W_est = W_Est_diag_reg(:,:,t);
            num_hat_t(t) = trace(W_est'*RR);
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

        den_hat_DL = mean(den1_hat_DL_t + den2_hat_t) - num_hat +(1/lambda_DL);
        gamma_hat_DL = real(num_hat/den_hat_DL);
        SE_hat_DL(nr,nq) = log2(1+gamma_hat_DL);
        
        %--- Expected SE Theoretical
        E = E_vec(:,:,nq);
        G = G_vec(:,:,nq);
        EQE = E*Qmatrix*E;
        GmE2 = G-E^2;
        Enum_const = (abs(trace(E*SS^2)))^2;
        Eden1_const = trace(EQE*SS*Rsum*SS) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum.*Qmatrix.*Qmatrix)) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum.*RR.*RR)) ...
            + (1 + alpha_R^2/(2*NR))*trace(GmE2*PP*SS^2*Ssum) + (alpha_R^2/(2*NR))*trace(GmE2*PP^3*Ssum);
        Eden1_const_DL = trace(EQE*SS*Rsum_DL*SS) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum_DL.*Qmatrix.*Qmatrix)) + (alpha_R^2/(2*NR))*trace(EQE*(Rsum_DL.*RR.*RR)) ...
            + (1 + alpha_R^2/(2*NR))*trace(GmE2*PP*SS^2*Ssum_DL) + (alpha_R^2/(2*NR))*trace(GmE2*PP^3*Ssum_DL);

        Eden2_const = 0;
        for l =1:number_of_cells
            Sjlu = diag(diag(Rmatrices(:,:,l,targetUser)));
            Eden2_const = Eden2_const + abs(trace(E*SS*Sjlu))^2+ (alpha_R^2/(2*NR))*sum(sum(E*Sjlu*(Qmatrix.*Qmatrix)*E*Sjlu))...
                + (alpha_R^2/(2*NR))*sum(sum(E*Sjlu*(RR.*RR)*E*Sjlu)) + (1 + alpha_R^2/(2*NR))*trace(GmE2*SS^2*Sjlu^2)...
                + (alpha_R^2/(2*NR))*trace(GmE2*PP^2*Sjlu^2);
        end

        gamma_theo = real((Enum_const)/(Eden1_const + Eden2_const - Enum_const));
        SE_theo(nr, nq) = prelogFactor(nr)*log2(1 + gamma_theo);

        gamma_theo_DL = real((Enum_const)/(Eden1_const_DL + Eden2_const - Enum_const +(1/lambda_DL)));
        SE_theo_DL(nr, nq) = log2(1 + gamma_theo_DL);

        num_verify1(nr, nq) = real(num_hat);
        num_verify2(nr, nq) = real(Enum_const);

        den1_verify1(nr, nq) =  real(mean(den1_hat_t));
        den1_verify2(nr, nq) = real(Eden1_const);

        den2_verify1(nr, nq) =  real(mean(den2_hat_t));
        den2_verify2(nr, nq) = real(Eden2_const);

        [nr nq toc]
    end
end

figure, plot(NR_vec, (1 - pilotSequenceLength/ulCoherenceSymbols)*SE_const_diag*ones(size(SE_theo(:,nq_sim))),'^-');
% ylim([0, SE_const_diag*1.1]);
hold on;
plot(NR_vec, SE_theo(:,nq_sim),'o-');
plot(NR_vec, SE_hat(:,nq_sim),'g*-');
hold off;

xlabel('N_{R}');
ylabel('UL Spectral efficiency per user [bits/s/Hz]');
legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Spectral efficiency (Numerical) - estimated covariance','Location','SouthEast');

figure, plot(NR_vec, SE_const_diag_DL*ones(size(SE_theo_DL(:,nq_sim))),'^-');
ylim([0, SE_const_diag*1.1]);
hold on;
plot(NR_vec, SE_theo_DL(:,nq_sim),'o-');
plot(NR_vec, SE_hat_DL(:,nq_sim),'g*-');
hold off;

xlabel('N_{R}');
ylabel('DL Spectral efficiency per user [bits/s/Hz]');
legend('Spectral efficiency - known covariance','Spectral efficiency (Theoretical) - estimated covariance','Spectral efficiency (Numerical) - estimated covariance','Location','SouthEast');

% figure, plot(NR_vec, num_verify1(3,:),'o-');
% ylim([0, max(num_verify1(:,3))*1.1]);
% hold on;
% plot(NR_vec, num_verify2(3,:),'*-');
% hold off;
% 
% figure, plot(NR_vec, den1_verify1(:,3),'o-');
% ylim([0, max(den1_verify1(:,3))*1.1]);
% hold on;
% plot(NR_vec, den1_verify2(:,3),'*-');
% hold off;
% 
% figure, plot(NR_vec, den2_verify1(:,3),'o-');
% ylim([0, max(den2_verify1(:,3))*1.1]);
% hold on;
% plot(NR_vec, den2_verify2(:,3),'*-');
% hold off;