function vs = fourier_evaluate_on_uniform(fs, nv)
    % fs needs form  [b_m, ..., b_1, a_0, a_1, ..., a_m]  for n = 2m+1
    % or  [b_m, ..., b_1, a_0, a_1, ..., a_m, a_m+1]  when n = 2m+2
    % representing  q(x) = a_0 + sum a_i cos(...) + sum b_i sin(...)
    % where the a_k & b_k can be complex column vectors
    nf = size(fs, 2);
    if nargin < 2
        nv = nf;
    end
    mf = floor((nf-1)/2);
    lf = nf-mf-1;
    fct = ceil(nf/nv);
    n = fct * nv;

    vse = zeros(size(fs, 1), n);
    vse(:, 1) = fs(:, mf+1);
    vse(:, 2:mf+1) = (fs(:, mf+2:2*mf+1) - 1i * fs(:, mf:-1:1)) / 2;
    if lf ~= mf
        vse(:, lf+1) = fs(:, nf) / 2;
        vse(:, n-mf) = vse(:, n-mf) + fs(:, nf) / 2; 
    end
    vse(:, n-mf+1:n) = (fs(:, 2*mf+1:-1:mf+2) + 1i * fs(:, 1:mf)) / 2;
    vse = ifft(vse, [], 2)*n;

    vs = vse(:, 1:fct:n);
end