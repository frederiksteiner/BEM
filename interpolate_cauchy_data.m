function dcd = interpolate_cauchy_data(mbp, ns, u, dnu)
    r = numel(mbp);

    dcd(r).n = 0;
    dcd(r).ts = [];
    dcd(r).ws = [];
    dcd(r).wls = [];

    dcd(r).gammas = [];
    dcd(r).dgammas = [];
    dcd(r).normdgammas = [];
    dcd(r).d2gammas = [];
    dcd(r).ngammas = [];

    dcd(r).rhoDs = [];
    dcd(r).rhoNs = [];

    for ll=1:r
        n = ns(ll);
        dcd(ll).n = n;
        dcd(ll).ts = circular_quadrature_nodes(n);
        dcd(ll).ws = circular_quadrature_weights(n);
        dcd(ll).wls = circular_quadrature_weights_log(n);

        dcd(ll).gammas = mbp(ll).gamma(n);
        dcd(ll).dgammas = mbp(ll).dgamma(n);
        dcd(ll).normdgammas = vecnorm(dcd(ll).dgammas);
        dcd(ll).d2gammas = mbp(ll).d2gamma(n);
        dcd(ll).ngammas = mbp(ll).ngamma(n);

        dcd(ll).rhoDs = u(dcd(ll).gammas);
        dcd(ll).rhoNs = dnu(dcd(ll).gammas, dcd(ll).ngammas);
    end
end