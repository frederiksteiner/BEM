function ws = circular_quadrature_weights(n)
    ws = zeros(1, n);
    ws(:) = 1/n;
end