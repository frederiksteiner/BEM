%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to interpolate onto quad mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zs, err] = cauchy_data_onto_quad_mesh(qm, dcd, kappa, u, gap)
    r = numel(dcd);
    zs = zeros(1, size(qm.vs, 2));
    zs(1, qm.bs < 1 | qm.bs > r) = potential_cauchy(dcd, kappa, qm.vs(:, qm.bs < 1 | qm.bs > r));
    for ll=1:r
        fs = fourier_compute_coefficients(dcd(ll).rhoDs);
        vs = fourier_evaluate_on_uniform(fs, qm.ns(ll));
        zs(1, qm.bs == ll) = vs(1 + qm.ds(1, qm.bs == ll));
    end

    if nargout > 1
        if nargin < 5
            gap = 2;
        end
        fltr = qm.bs < 1 & qm.ds > gap;
        err = max(abs(zs(fltr) - u(qm.vs(:, fltr))));
    end
end
