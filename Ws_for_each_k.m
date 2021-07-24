function Ws_for_each_k(NR_vec, NQ_vec, nq_set, R_sqrt_root, number_of_cells, number_of_users, number_of_antennas, targetCell, kk, pilotSequenceLength, mu, phi, alpha_R, alpha_Q, Rb, TT, Pb, interm_folder)
    % compute a set of TT =500 number of covariance estimates for different NQ and NR values 
    for nr = 1:length(NR_vec)
        for nq = nq_set
            clear W_Est;
            clear W_Est_diag;
            clear W_Est_diag_reg;
            NR = NR_vec(nr);
            NQ = NQ_vec(nq);
%             tempp = 0;
%             R_est_mean = zeros(number_of_antennas,number_of_antennas);
%             Q_est_mean = zeros(number_of_antennas,number_of_antennas);

            for trail = 1:TT
                [R_est, Q_est] = getRQest(R_sqrt_root, number_of_cells, number_of_users, number_of_antennas, targetCell, kk, pilotSequenceLength, NR, NQ, mu, phi, alpha_R, Rb);
    %             R_est_mean = R_est_mean + (1/TT)*R_est;
    %             Q_est_mean = Q_est_mean + (1/TT)*Q_est;

                S_est = diag(diag(R_est));
                P_est = diag(diag(Q_est));

                W_Est(:,:,trail) = R_est/Q_est;
                W_Est_diag(:,:,trail) = real(S_est/P_est);
                W_Est_diag_reg(:,:,trail) = real(S_est/(alpha_Q*P_est+(1-alpha_Q)*Pb));
            end
    %         norm(R_est_mean-R_bar)/norm(R_bar)*100
    %         norm(Q_est_mean-Qmatrix)/norm(Qmatrix)*100

            save(strcat(interm_folder,'\RQest_',string(nr), '_', string(nq),'_', string(kk), '.mat'), 'W_Est', 'W_Est_diag', 'W_Est_diag_reg');

            [nr  nq toc]
        end
    end
end