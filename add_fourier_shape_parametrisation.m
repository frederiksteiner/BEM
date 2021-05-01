function mbp = add_fourier_shape_parametrisation(mbp, vs, no)
%adds to mbp trigonometric poly of the N points in vs
    if isempty(mbp)
        l = 1;
    else
        l = size(mbp,2) + 1;
    end

    fs = fourier_compute_coefficients(vs);
    dfs = fourier_differentiate_coefficients(fs);
    d2fs = fourier_differentiate_coefficients(dfs);

    mbp(l).n0 = no;
    mbp(l).gamma = @(n) fourier_evaluate_on_uniform(fs, n);
    mbp(l).dgamma = @(n) fourier_evaluate_on_uniform(dfs, n);
    mbp(l).d2gamma = @(n) fourier_evaluate_on_uniform(d2fs, n);
    mbp(l).ngamma = @(n) -mbp(l).n0*[0,-1;1,0]*normalize(mbp(l).dgamma(n), 'norm');
end

