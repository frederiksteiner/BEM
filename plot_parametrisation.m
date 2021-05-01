function plot_parametrisation(mbp, np, h, dn)
    if nargin < 4
        dn = 1e-1;
    end
    if nargin < 3
        h = 1e-2;
    end
    if nargin < 2
        np = 2^8;
    end

    r = numel(mbp);
    hold on;
    for q=1:r
        xs = mbp(q).gamma(np);
        dxs = mbp(q).dgamma(np);
        nxs = mbp(q).ngamma(np);

        pdxs = [xs(1, :)-h*dxs(1, :); xs(1, :)+h*dxs(1, :)];
        pdys = [xs(2, :)-h*dxs(2, :); xs(2, :)+h*dxs(2, :)];

        pnxs = [xs(1, :); xs(1, :)+dn*nxs(1, :)];
        pnys = [xs(2, :); xs(2, :)+dn*nxs(2, :)];

        plot(pdxs, pdys, 'r', ...
             pnxs, pnys, 'b', ...
             xs(1, :), xs(2, :), 'k');
    end
    hold off;
    axis 'equal';
end