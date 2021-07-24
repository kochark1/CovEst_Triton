function f = pdfChisqr(x, k)
f = x.^(k/2-1).* exp(-x/2)./(2^(k/2)*gamma(k/2));
end