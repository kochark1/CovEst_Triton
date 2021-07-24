function [Rmatrices, Rsum, Rsum_DL, R_sqrt_root, Qmatrix, Wmatrix, W_bar, SE_const, SE_const_DL, W_diag, W_bar_diag, SE_const_diag, SE_const_diag_DL] = ConfigRQWC(L, K, M, targetCell, targetUser, mu, lambda_DL, alpha_R, Rb, pilotSequenceLength)
    % - Generate user locations and parameters for the covariance matrices - %
    center  = zeros(1,L);
    userRadius                  = 120; % User radius in m (Users located on a ring around the BS)
    cellRadius                  = 150; % Half of the distance between BSs in m
    for ll = 1:L
%         center(ll) = ceil((ll-1)/6) * 2 *  cellRadius * exp(1i * (pi/3 * mod(ll-1,6) ) ); % Location of the remaining BSs
        center(ll) = ceil((ll-1)/6) * 2 * sqrt(3)/2 *  cellRadius * exp(1i * (pi/3 * mod(ll-1,6) ) ); % Location of the remaining BSs
        userLocations  = userRadius * exp(1i * 2 * pi / K * (0:K-1)) + center(ll); % User locations in cell ii
        rxPowerInDb    = snrModel(userLocations); % Received power in dB corresponding to each user in cell ii
        rxPower(ll,:)  = 10.^(rxPowerInDb/10); % Received power corresponding to each user in cell ii
        meanAngle(ll,:)= angle(userLocations) * 180 / pi; % Mean angle of channel cluster
    end

    angleSpread     = 20 * ones(L,K); % Angular spread of the channel cluster in Degrees ( sys.numCell x sys.numUser )
    normAntennaSpacing          = 1/2; % Antenna spacing normalized by wavelength
    sumR = 0;
    Q_temp = 0;
    I = eye(M);
    integralSpacing = 0.01;
    
    for ll = 1:L
        for kk = 1:K
            % - Compute covariance matrix of user (ll,kk) at the reference BS - % 
            userMeanAngle   = meanAngle(ll,kk);
            userAngleSpread = angleSpread(ll,kk);
            
            thetaRange      = -180:integralSpacing:360; % Account for wrap around
            pTheta          = 1/userAngleSpread * ((thetaRange >= userMeanAngle - userAngleSpread/2) & (thetaRange <= userMeanAngle + userAngleSpread/2));

            aVector         = exp(1i * 2 * pi * normAntennaSpacing * (0:M-1).' * cos(pi/180 * thetaRange(pTheta>0)));
            Rmatrices(:,:,ll,kk)         = rxPower(ll,kk) * (aVector * diag(pTheta(pTheta>0)) * aVector') * integralSpacing;
            R_sqrt_root(:,:,ll,kk)    = (Rmatrices(:,:,ll,kk))^(1/2);

            sumR = sumR + Rmatrices(:,:,ll,kk);
        end
        Q_temp = Q_temp + Rmatrices(:,:,ll,targetUser);
    end
    Qmatrix = Q_temp + (1/(pilotSequenceLength*mu))*I;
    Q_inv = inv(Qmatrix);
    P_inv = inv(diag(diag(Qmatrix)));
    R_bar = alpha_R*Rmatrices(:,:,targetCell,targetUser)+ (1-alpha_R)*Rb;
    W_bar = R_bar*Q_inv;
    W_bar_diag = diag(diag(R_bar))*P_inv;
    Rsum = sumR+(1/mu)*I;
    Rsum_DL = sumR;
    Rtarget = Rmatrices(:,:,targetCell, targetUser);
    
    Wmatrix = Rtarget*Q_inv;
    SE_const = computeCap(Wmatrix, Rmatrices, Qmatrix, L, targetCell, targetUser, Rsum, 0, 0);
    SE_const_DL = computeCap(Wmatrix, Rmatrices, Qmatrix, L, targetCell, targetUser, Rsum_DL, 1, lambda_DL);
    
    W_diag = diag(diag(Rtarget))*P_inv;
    SE_const_diag = computeCap(W_diag, Rmatrices, Qmatrix, L, targetCell, targetUser, Rsum, 0, 0);
    SE_const_diag_DL = computeCap(W_diag, Rmatrices, Qmatrix, L, targetCell, targetUser, Rsum_DL, 1, lambda_DL);
end