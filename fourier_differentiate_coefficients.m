function dfs = fourier_differentiate_coefficients(fs)
    % fs needs form  [b_m, ..., b_1, a_0, a_1, ..., a_m]  for n = 2m+1
    % or  [b_m, ..., b_1, a_0, a_1, ..., a_m, a_m+1]  when n = 2m+2
    % representing  q(x) = a_0 + sum a_i cos(...) + sum b_i sin(...)
    % where the a_k & b_k can be complex column vectors
    nf = size(fs, 2);
    if mod(nf, 2) == 0
        if all(fs(:, nf) == 0)
            fs = fs(:, 1:nf-1);
            nf = nf-1;
        else
            fs = [zeros(size(fs, 1), 1), fs];
            nf = nf+1;
        end
    end
    mf = floor((nf-1)/2);
    dfs = (2*pi) * (-mf:mf) .* fs(:, nf:-1:1);
end