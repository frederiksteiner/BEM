function ws = circular_quadrature_weights_log(n)
    ws = zeros(1, n);
    m = floor(n/2);
    af = 1+mod(n, 2);
    vals = 2*pi/n*(1-(1:n));
    for kk=1:m-1
        ws = ws + 2/kk * cos(kk*vals);
    end
    if m > 0
        ws = ws + af/m * cos(m*vals);
    end
    ws = (log(4) + ws)/n;
end