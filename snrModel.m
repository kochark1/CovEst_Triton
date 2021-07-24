function snr_dB = snrModel(userLocations)
    model_id = 1;
    if model_id == 0
        snr_dB =  (71.89 - 37.6*log10(abs(userLocations))); % Bad modedl
    elseif model_id == 1
        snr_dB =  (93.95 - pathLoss(abs(userLocations))); % 6dBm as transmit power
    end
end

function PL = pathLoss(userLocations)
    h_BS = 25;
    h_UT = 1.5;
    h_diff = h_BS-h_UT;
    d_3D = sqrt(userLocations.^2 + h_diff^2);
    f=3.4;
    
    PL = 32.4 + 20*log10(f) + 30*log10(d_3D);
end