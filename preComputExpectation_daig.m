function [Enum, Eden1, Eden2] = preComputExpectation_daig(R, Q, W_diag, W_bar, Rsum, L, targetCell, targetUser, alpha_R)
%     Enum is the term in eq 36 without multiplicative factor
%     Eden1 is group of terms in eq 37 without multiplicative factor
%     Eden2 is group of terms in eq 38 without multiplicative factor

    RR = R(:,:,targetCell,targetUser);
    SS = diag(diag(RR));
    PP = diag(diag(Q));
    P_inv = inv(PP);
    Ssum = diag(diag(Rsum));
    Eden2 = zeros(4,1);
    for l =1:L
        Sjlu = diag(diag(R(:,:,l,targetUser)));
        Wjlu = Sjlu*P_inv;
        Wl_sqr = Wjlu^2;
        Eden2(1,:) = Eden2(1,:) + abs(trace(W_bar'*Sjlu))^2; % term in eq 38 that has (only) kappa_3 as multiplicative factor
        Eden2(2,:) = Eden2(2,:) + (alpha_R^2/2)*(sum(sum(Wjlu*(Q.*Q)*Wjlu)) + sum(sum(Wjlu*(RR.*RR)*Wjlu))); % terms in eq 38 that have kappa_3/NR as multiplicative factor
        Eden2(3,:) = Eden2(3,:) + trace(W_bar^2*Sjlu^2); % term in eq 38 that has (only) kappa_4 as multiplicative factor
        Eden2(4,:) = Eden2(4,:) + (alpha_R^2/2)*(sum(sum(Wjlu^2*PP^2)) + sum(sum(Wjlu^2*Sjlu^2))); % terms in eq 38 that have kappa_4/NR as multiplicative factor
    end
    Enum = (abs(trace(W_bar'*RR)))^2;  % square of term in eq 36 without NQ/(NQ-1)
    Eden1(1,:) = trace(W_bar*Q*W_bar'*Rsum); % term in eq 37 that has (only) kappa_3 as multiplicative factor
    Eden1(2,:) = (alpha_R^2/2)*(trace(P_inv*Q*P_inv*(Rsum.*RR.*RR))+trace(P_inv*Q*P_inv*(Rsum.*Q.*Q))); % terms in eq 37 that have kappa_3/NR as multiplicative factor
    Eden1(3,:) = trace(W_bar*PP*W_bar'*Ssum); % term in eq 37 that has (only) kappa_4 as multiplicative factor
    Eden1(4,:) = (alpha_R^2/2)*(trace(Ssum*PP)+trace(W_diag*Ssum*SS)); % terms in eq 37 that have kappa_4/NR as multiplicative factor
end