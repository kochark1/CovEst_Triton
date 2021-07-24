function [SE_const] = computeCap(Wmatrix, Rmatrices, Qmatrix, L, targetCell, targetUser, Rsum_UL_DL, dl, lambda_DL)
    den2_const = 0;
    for l =1:L
        den2_const = den2_const + abs(trace(Wmatrix'*Rmatrices(:,:,l,targetUser)))^2;
    end
    
    num_const = (abs(trace(Wmatrix'*Rmatrices(:,:,targetCell,targetUser))))^2;
    den1_const = trace(Wmatrix*Qmatrix*Wmatrix'*Rsum_UL_DL);
    
    den_const = real(den1_const + den2_const-num_const);
    if dl == 1
        den_const = den_const+(1/lambda_DL);
    end
    
    
    gamma_snr_const = num_const/den_const;
    SE_const = log2(1+gamma_snr_const); % SE without prelog factor
end