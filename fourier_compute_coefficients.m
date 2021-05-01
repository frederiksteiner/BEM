function fs = fourier_compute_coefficients(vs, nf)
    % fs takes form  [b_m, ..., b_1, a_0, a_1, ..., a_m]  for n = 2m+1
    % or  [b_m, ..., b_1, a_0, a_1, ..., a_m, a_m+1]  when n = 2m+2
    % representing  q(x) = a_0 + sum a_i cos(...) + sum b_i sin(...)
    % where the a_k & b_k can be complex column vectors
    nv = size(vs, 2);
    if nargin < 2
        nf = nv;
    end
    mv = floor((nv-1)/2);
    mf = floor((nf-1)/2);
    m = min(mv, mf);
    l = min(nv-mv-1, nf-mf-1);

    vs = fft(vs, [], 2);
    vs(:, [mv+1:nv, 1:mv]) = vs;
    vs(:, 1:mv) = -1i * vs(:, 1:mv);
    vs(:, [1:mv, mv+2:2*mv+1]) = vs(:, [1:mv, mv+2:2*mv+1]) + 1i * vs(:, [2*mv+1:-1:mv+2, mv:-1:1]);

    fs = zeros(size(vs, 1), nf);
    fs(:, mf+1-m:mf+1+l) = vs(:, mv+1-m:mv+1+l) / nv;
end