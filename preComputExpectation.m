function [Enum, Eden1, Eden2] = preComputExpectation(R, Q, W, W_bar, Rsum, L, M, targetCell, targetUser, alpha_R)
%     Enum is the term in eq 24 without multiplicative factor
%     Eden1 is group of terms in eq 25 without multiplicative factor
%     Eden2 is group of terms in eq 26 without multiplicative factor
    RR = R(:,:,targetCell,targetUser);
    Q_inv = inv(Q);
    Eden2 = zeros(4,1);
    for l =1:L
        Rjlu = R(:,:,l,targetUser);
        Wjlu = Rjlu*Q_inv;
        Wl_sqr = Wjlu^2;
        Eden2(1,:) = Eden2(1,:) + abs(trace(W_bar'*Rjlu))^2; % term in eq 26 that has (only) kappa_2 as multiplicative factor
        Eden2(2,:) = Eden2(2,:) + (alpha_R^2/2)*(trace(Wjlu*Q*Wjlu'*Q) + trace(Wjlu*RR*Wjlu'*RR)); % terms in eq 26 that have kappa_2/NR as multiplicative factor
        Eden2(3,:) = Eden2(3,:) + trace(W_bar'*W_bar*(Q*Wjlu')*Wjlu*Q); % term in eq 26 that has kappa_1/NQ as multiplicative factor
        Eden2(4,:) = Eden2(4,:) + (alpha_R^2/2)*(M*trace(Wl_sqr*Q^2)+ trace(W)*trace(Wl_sqr*Q*RR)); % terms in eq 26 that have kappa_1/(NQ*NR) as multiplicative factor
    end
    Enum = (abs(trace(W_bar'*RR)))^2; % square of term in eq 24 without NQ/(NQ-M)
    Eden1(1,:) = trace(W_bar*Q*W_bar'*Rsum); % term in eq 25 that has (only) kappa_1 as multiplicative factor
    Eden1(2,:) = (alpha_R^2/2)*(M*trace(Rsum*Q) + trace(W)*trace(Rsum*RR)); % terms in eq 25 that have kappa_1/NR as multiplicative factor
end