function dcd = nystroem_cauchy_bw(mbp, ns, g,kappa, eta)
    % calculate dcd
        r = numel(mbp);

    % pre-allocate dcd fields
    dcd(r).n = 0;
    dcd(r).ts = [];
    dcd(r).ws = [];
    dcd(r).wls = [];

    dcd(r).gammas = [];
    dcd(r).normdgammas = [];
    dcd(r).dgammas = [];
    dcd(r).d2gammas = [];
    dcd(r).ngammas = [];

    dcd(r).rhoDs = [];
    dcd(r).rhoNs = [];

    % compute dcd fields
    for ll=1:r
        n = ns(ll);
        llts = circular_quadrature_nodes(n);
        dcd(ll).n = n;
        dcd(ll).ts = llts;
        dcd(ll).ws = circular_quadrature_weights(n);
        dcd(ll).wls = circular_quadrature_weights_log(n);

        dcd(ll).gammas = mbp(ll).gamma(n);
        dcd(ll).dgammas = mbp(ll).dgamma(n);
        dcd(ll).normdgammas = vecnorm(dcd(ll).dgammas);
        dcd(ll).d2gammas = mbp(ll).d2gamma(n);
        dcd(ll).ngammas = mbp(ll).ngamma(n);
        dcd(ll).rhoDs = g(dcd(ll).gammas);
               
    end   
    
    matlength = sum(ns);
    A = zeros(matlength);
    plD = zeros(matlength,1);
%     plN = zeros(matlength,1);
%     b = zeros(matlength,1);

    xrange = 1;
    
    for m = 1:r
        yrange = 1;
        nm = dcd(m).n;
        mgammas = dcd(m).gammas;
            for l = 1:r
                nl = dcd(l).n;
                if l ~= m
                    A(yrange:(yrange + nl)-1,xrange:(xrange +nm-1)) = ...
                        - calc_Klm(dcd,l,m,kappa) - 1i*eta*calc_Vlm(dcd,l,m, kappa);
                else
                    A(yrange:(yrange + nl)-1,xrange:(xrange +nm-1)) = eye(nm)...
                        - calc_Klm(dcd,l,m,kappa) - 1i*eta*calc_Vlm(dcd,l,m, kappa);
                end
                yrange = yrange + nl;
            end
        plD(xrange:(xrange + nm-1)) = g(mgammas).';
        xrange = xrange + nm;
    end
    plN = A\plD;
    xrange = 1;
    for i = 1:r
        n = dcd(i).n;
        dcd(i).rhoDs = plD(xrange:(xrange + n-1)).';
        dcd(i).rhoNs = plN(xrange:(xrange + n-1)).';
        xrange = xrange + n;
    end
end

