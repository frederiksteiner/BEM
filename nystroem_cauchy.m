function dcd = nystroem_cauchy(mbp, ns, g, h,kappa, btf)
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
               
    end   
    
    matlength = sum(ns);
    A = zeros(matlength);
    plD = zeros(matlength,1);
    plN = zeros(matlength,1);
    b = zeros(matlength,1);

    xrange = 1;
    for m = 1:r
        yrange = 1;
        nm = dcd(m).n;
        % wenn true, dann Dirichlet, sonst Neumann bekannt
        mgammas = dcd(m).gammas;
        mngammas = dcd(m).ngammas;
        if btf(m)
            for l = 1:r
                nl = dcd(l).n; 
                A(yrange:(yrange + nl)-1,xrange:(xrange +nm-1)) = calc_Vlm(dcd,l,m, kappa);
                b(yrange:(yrange + nl-1)) =b(yrange:(yrange + nl-1)) + calc_Klm(dcd,l,m, kappa) * g(mgammas).';
                yrange = yrange + nl;
            end
        else
            for l = 1:r
                nl = dcd(l).n; 
                A(yrange:(yrange + nl)-1,xrange:(xrange +nm-1)) = - calc_Klm(dcd,l,m, kappa);
                b(yrange:(yrange + nl-1)) =b(yrange:(yrange + nl-1)) - calc_Vlm(dcd,l,m, kappa) * h(mgammas, mngammas).';
                yrange = yrange + nl;
            end
        end
        plD(xrange:(xrange + nm-1)) = g(mgammas).';
        plN(xrange:(xrange + nm-1)) = h(mgammas,mngammas).';
        xrange = xrange + nm;
    end

    x = A\b;
    xrange = 1;
    for i = 1:r
        n = dcd(i).n;
        if btf(i) 
            dcd(i).rhoDs = plD(xrange:(xrange + n-1)).';
            dcd(i).rhoNs = x(xrange:(xrange + n-1)).';
        else
            dcd(i).rhoDs = x(xrange:(xrange + n-1)).';
            dcd(i).rhoNs = plN(xrange:(xrange + n-1)).';
        end
        xrange = xrange + n;
    end
end
