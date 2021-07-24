
figure, surf(NQ_vec, NR_vec, prelogFactor_LMMSE*SE_const*ones(size(SE_theo)));
hold on;
surf(NQ_vec, NR_vec, SE_theo);
hold off;
xlabel('N_{Q}');
ylabel('N_{R}');
zlabel('UL Spectral efficiency per user [bits/s/Hz]');
legend('LMMSE','LMMSE-type (Theoretical)','Location','NorthEast');
print('ulSE_3D', '-dpdf', '-bestfit');

figure, surf(NQ_vec, NR_vec, SE_const_DL*ones(size(SE_theo)));
hold on;
surf(NQ_vec, NR_vec, SE_theo_DL);
hold off;
xlabel('N_{Q}');
ylabel('N_{R}');
zlabel('DL Spectral efficiency per user [bits/s/Hz]');
legend('LMMSE','LMMSE-type (Theoretical)','Location','NorthEast');
print('dlSE_3D', '-dpdf', '-bestfit');