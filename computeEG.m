function [E, G] = computeEG(p_vec, NQ_vec, alpha_Q, pb_vec)
    p_vec = real(p_vec);
    NQ_index = 1;
    E = zeros(length(p_vec), length(p_vec), length(NQ_vec));
    G = E;
    for NQ = NQ_vec
        a = (1/(2*NQ)) * alpha_Q * p_vec;
        b = (1-alpha_Q)*pb_vec;
        for ii = 1:length(a)
            [E(ii,ii,NQ_index), G(ii,ii,NQ_index)] = numExpectation(a(ii), b(ii), 2*NQ);
        end
        NQ_index = NQ_index + 1;
    end
end

function [E, G] = numExpectation(a, b, dg)
% Numerical expectation of 1/(ax+b) and 1/(ax+b)^2 over x, where x is chi
% squared distribution with dg degrees of freedom. But we approximate the
% chi squared pdf with normal pdf since dg is large.

    mn = dg; %mean
    sd = sqrt(2*dg); %standard deviation
    x_resolution = 8*sd/200;
    
    x = mn-(5*sd):x_resolution:mn+(5*sd);
    
    x_pdf = normpdf(x, mn, sd);
    x_pdf = pdfChisqr(x, mn);
    if isnan(sum(x_pdf)) || isinf(sum(x_pdf))
        x_pdf = normpdf(x, mn, sd);
    end
    
    
    E = sum((1./(a*x + b)).*x_pdf*x_resolution);
    G = sum((1./(a*x + b).^2).*x_pdf*x_resolution);
end
