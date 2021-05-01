function [zs, err] = cauchy_data_onto_quad_mesh_bw(qm, dcd, kappa, u, eta, gap)
    r = numel(dcd);
    zs = zeros(1, size(qm.vs, 2));
    zs(1, qm.bs < 1 | qm.bs > r) = potential_cauchy_bw(dcd, kappa, qm.vs(:, qm.bs < 1 | qm.bs > r), eta);
    for ll=1:r
        fs = fourier_compute_coefficients(dcd(ll).rhoDs);
        vs = fourier_evaluate_on_uniform(fs, qm.ns(ll));
        zs(1, qm.bs == ll) = vs(1 + qm.ds(1, qm.bs == ll));
    end

    if nargout > 1
        if nargin < 6
            gap = 2;
        end
        fltr = qm.bs < 1 & qm.ds > gap;
        err = max(abs(zs(fltr) - u(qm.vs(:, fltr))));
    end
end