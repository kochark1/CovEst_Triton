% close all;
figure, plot(NR_vec, prelogFactor_LMMSE*SE_const*ones(size(SE_theo(:,nq_sim))),'b^-');
ylim([0, prelogFactor_LMMSE*SE_const*1.1]);
hold on;
plot(NR_vec, SE_theo(:,nq_sim),'go-');
plot(NR_vec, SE_hat,'r*-');

plot(NR_vec, prelogFactor_LMMSE*SE_const_diag*ones(size(SE_theo_diag(:,nq_sim))),'b^--');
plot(NR_vec, SE_theo_diag(:,nq_sim),'g+-');
plot(NR_vec, SE_hat_diag,'ch-.');

plot(NR_vec, prelogFactor_LMMSE*SE_const_diag*ones(size(SE_theo_diag_reg(:,nq_sim))),'m*--');
plot(NR_vec, SE_theo_diag_reg(:,nq_sim),'ro:');
plot(NR_vec, SE_hat_diag_reg,'k.-.');

hold off;

xlabel('N_{R}');
ylabel('UL Spectral efficiency per user [bits/s/Hz]');
legend('SE - known R and Q','SE (Theoretical) - estimated R and Q','SE (Numerical) - estimated R and Q', ...
       'SE - known S and P','SE (Theoretical) - estimated S and P','SE (Numerical) - estimated S and P', ...
       'SE - known S and P','SE (Theoretical) - estimated S and regularized P','SE (Numerical) - estimated S and regularized P','Location','southEast');

%-----------------------

figure, plot(NR_vec, SE_const_DL*ones(size(SE_theo_DL(:,nq_sim))),'b^-');
ylim([0, SE_const_DL*1.1]);
hold on;
plot(NR_vec, SE_theo_DL(:,nq_sim),'go-');
plot(NR_vec, SE_hat_DL,'r*-');

plot(NR_vec, SE_const_diag_DL*ones(size(SE_theo_diag_DL(:,nq_sim))),'b^--');
plot(NR_vec, SE_theo_diag_DL(:,nq_sim),'g+-');
plot(NR_vec, SE_hat_diag_DL,'ch-.');

plot(NR_vec, SE_const_diag_DL*ones(size(SE_theo_diag_reg_DL(:,nq_sim))),'m*--');
plot(NR_vec, SE_theo_diag_reg_DL(:,nq_sim),'ro:');
plot(NR_vec, SE_hat_diag_reg_DL,'k.-.');

hold off;

xlabel('N_{R}');
ylabel('DL Spectral efficiency per user [bits/s/Hz]');
legend('SE - known R and Q','SE (Theoretical) - estimated R and Q','SE (Numerical) - estimated R and Q', ...
       'SE - known S and P','SE (Theoretical) - estimated S and P','SE (Numerical) - estimated S and P', ...
       'SE - known S and P','SE (Theoretical) - estimated S and regularized P','SE (Numerical) - estimated S and regularized P','Location','southEast');