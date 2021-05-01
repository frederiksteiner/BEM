function qm = generate_quad_mesh(mbp, h, xl, r)
    p = numel(mbp);

    d = struct([]);
    xne = [NaN; NaN];
    xsw = [NaN; NaN];
    for q=p:-1:1
        % compute a first discrete boundary with 2^x points
        n = max(2^ceil(log2(ceil(norm(mbp(q).dgamma(1)) / h))), 8);
        %ks = 0:n-1;
        gammas = mbp(q).gamma(n);
        deltas = gammas(:, [2:n, 1]) - gammas;
        hr = max(vecnorm(deltas));

        % refine it until it is a sufficiently resolved discrete boundary
        % vis-a-vis the meshsize h of the lattice
        while hr > h
            n = 2*n;
            %ks = 0:n-1;
            gammas = mbp(q).gamma(n);
            deltas = gammas(:, [2:n, 1]) - gammas;
            hr = max(vecnorm(deltas));
        end
        % pad the 'just' sufficient refinement by a factor 2^r
        if nargin < 4
            r = 4;
        end
        d(q).n = 2^r * n;
        d(q).ks = 0:d(q).n-1;
        d(q).gammas = mbp(q).gamma(d(q).n);

        xne = max(xne, max(d(q).gammas, [], 2));
        xsw = min(xsw, min(d(q).gammas, [], 2));
    end

    % compute the base stuff for the lattice
    if nargin < 3 || isempty(xl)
        xl = [0; 0];
        nxl = 0;
        for q=1:p
            xl = xl + sum(d(q).gammas, 2);
            nxl = nxl + d(q).n;
        end
        xl = xl / nxl;
    end
    ms = ceil((xl - xsw) / h)+4;
    xo = xl - h * ms;
    ms = (ms + ceil((xne-xl) / h)+4).';

    for q=1:p
        % discretise the discrete boundary onto the coarse lattice (2*h)
        gammaijs = round((d(q).gammas - xo) / (2*h));
        isnewijs = any(gammaijs ~= gammaijs(:, [d(q).n, 1:d(q).n-1]));

        % compute the lattice boundary and a "best" mapping onto the
        % discrete boundary
        rgammaijs = gammaijs(:, isnewijs);
        nr = size(rgammaijs, 2);
        rsks = d(q).ks(:, isnewijs);
        reks = d(q).ks(:, isnewijs);
        reks(reks == 0) = d(q).n;
        reks = reks(:, [2:nr, 1]) - 1;
        rlks = reks - rsks;
        rlks(rlks < 0) = rlks(rlks < 0) + d(q).n;
        rn = 2*d(q).n;
        rks = 2*rsks + rlks;
        rks(rks >= rn) = rks(rks >= rn) - rn;

        % then we double up to get onto lattice (h)
        rgammaijs = reshape([2*rgammaijs; rgammaijs+rgammaijs(:, [2:nr, 1])], [], 2*nr);
        rlks = rks(:, [2:nr, 1]) - rks;
        rlks(rlks < 0) = rlks(rlks < 0) + rn;
        rn = 2*rn;
        reks = 2*rks;
        rsks = reks+rlks;
        rsks(rsks >= rn) = rsks(rsks >= rn) - rn;
        rks = reshape([reks; rsks], [], 2*nr);
        crgammas = mbp(q).gamma(rn); 
        rgammas = crgammas(:, 1+rks);

        % remove ears and hairs (ieratively)
        haschanged = true;
        while haschanged
            % first ears
            nr = size(rgammaijs, 2);
            idxl = [nr, 1:nr-1];
            idxr = [2:nr, 1];
            isear = all(abs(rgammaijs(:, idxl)-rgammaijs(:, idxr)) <= 1) ...
                    & any(abs(rgammaijs(:, idxl)-rgammaijs(:, idxr)) == 1);
            keep = ~isear;
            rgammaijs = rgammaijs(:, keep);
            rks = rks(:, keep);
            rgammas = rgammas(:, keep);

            % then hairs
            nr = size(rgammaijs, 2);
            idxl = [nr, 1:nr-1];
            idxr = [2:nr, 1];
            ishair = all(rgammaijs(:, idxl) == rgammaijs(:, idxr));
            rks(:, idxl(1, ishair)) = rks(:, ishair);
            rgammas(:, idxl(1, ishair)) = rgammas(:, ishair);
            keep = ~(ishair | ishair(1, idxl));
            rgammaijs = rgammaijs(:, keep);
            rks = rks(1, keep);
            rgammas = rgammas(:, keep);

            haschanged = any(ishair) || any(isear);
        end

        d(q).rgammaijs = rgammaijs;
        d(q).rks = rks;
        d(q).rn = rn;
        d(q).rgammas = rgammas;
    end

    % compute a lattice and further stuffs
    [latx, laty] = ndgrid(xo(1)+h*(0:ms(1)), xo(2)+h*(0:ms(2)));
    latb = zeros(size(latx));
    latb(2:ms(1), 2:ms(2)) = -1;
    nd = 0;
    for q=1:p
        bijs = 1 + d(q).rgammaijs(1, :) + (ms(1)+1) * d(q).rgammaijs(2, :);
        latb(bijs) = 0;
        latx(bijs) = d(q).rgammas(1, :);
        laty(bijs) = d(q).rgammas(2, :);

        nd = nd + size(d(q).rgammaijs, 2);
    end

    % now we can smooth the lattice towards fitting the boundary
    tosmooth = (latb ~= 0);
    fltr = [0, 1, 0; 1, 4, 1; 0, 1, 0];
    fltr = fltr / sum(fltr, 'all');
    for jjj=1:10
        lats = conv2(latx, fltr, 'same');
        latx(tosmooth) = lats(tosmooth);
        lats = conv2(laty, fltr, 'same');
        laty(tosmooth) = lats(tosmooth);
    end

    % double up to the fine lattice (h/2)  ->  cell diagonal should be <= h
    % this is the base for cutting triangles into 3 quadrangles !
    ns = zeros(size(mbp));
    flatx = zeros(2*ms+1);
    flaty = zeros(2*ms+1);
    flatd = -ones(2*ms+1);
    flatb = -ones(2*ms+1);
    flatb([1, 2*ms(1)+1], :) = p+1;
    flatb(:, [1, 2*ms(2)+1]) = p+1;
    flatx(1:2:2*ms(1)+1, 1:2:2*ms(2)+1) = latx;
    flaty(1:2:2*ms(1)+1, 1:2:2*ms(2)+1) = laty;
    flatx(2:2:2*ms(1), 1:2:2*ms(2)+1) = 0.5*(flatx(1:2:2*ms(1)-1, 1:2:2*ms(2)+1) + flatx(3:2:2*ms(1)+1, 1:2:2*ms(2)+1));
    flaty(2:2:2*ms(1), 1:2:2*ms(2)+1) = 0.5*(flaty(1:2:2*ms(1)-1, 1:2:2*ms(2)+1) + flaty(3:2:2*ms(1)+1, 1:2:2*ms(2)+1));
    flatx(:, 2:2:2*ms(2)) = 0.5*(flatx(:, 1:2:2*ms(2)-1) + flatx(:, 3:2:2*ms(2)+1));
    flaty(:, 2:2:2*ms(2)) = 0.5*(flaty(:, 1:2:2*ms(2)-1) + flaty(:, 3:2:2*ms(2)+1));

    % compute stuff on the lattice (h/2)
    for q=1:p
        nr = size(d(q).rgammaijs, 2);
        fgammaijs = reshape([2*d(q).rgammaijs; d(q).rgammaijs+d(q).rgammaijs(:, [2:nr, 1])], [], 2*nr);
        flks = d(q).rks(:, [2:nr, 1]) - d(q).rks;
        flks(flks < 0) = flks(flks < 0) + d(q).rn;
        fn = 2*d(q).rn;
        feks = 2*d(q).rks;
        fsks = feks+flks;
        fsks(fsks >= fn) = fsks(fsks >= fn) - fn;
        fks = reshape([feks; fsks], [], 2*nr);
        cfgammas = mbp(q).gamma(fn);
        fgammas = cfgammas(:, 1+fks);
        cfngammas = mbp(q).ngamma(fn);
        fngammas = cfngammas(:, 1+fks);

        ns(q) = fn;

        % set boundary points
        bijs = 1 + fgammaijs(1, :) + (2*ms(1)+1) * fgammaijs(2, :);
        flatb(bijs) = q;
        flatx(bijs) = fgammas(1, :);
        flaty(bijs) = fgammas(2, :);
        flatd(bijs) = fks;

        % set some inner points
        ndelta = round(fngammas ./ max(abs(fngammas), [], 1));
        bijs = 1 + fgammaijs(1, :) - ndelta(1, :) + (2*ms(1)+1) * (fgammaijs(2, :) - ndelta(2, :));
        flatb(bijs) = max(flatb(bijs), 0);
    end

    % now we can smooth the lattice towards fitting the boundary
    tosmooth = (flatb < 1);
    fltr = [0, 1, 0; 1, 4, 1; 0, 1, 0];
    fltr = fltr / sum(fltr, 'all');
    for jjj=1:10
        flats = conv2(flatx, fltr, 'same');
        flatx(tosmooth) = flats(tosmooth);
        flats = conv2(flaty, fltr, 'same');
        flaty(tosmooth) = flats(tosmooth);
    end

    % find inner
    % for flatb: -1 is outside, 0 is inside, and {1, ..., p} are boundary
    % for flatd: on boundary is numerator for boundary parameterisation,
    %            else is grid-distance to boundary
    flatb(flatb == p+1) = -1;
    flatm = bwdist(flatb > 0);
    flatd(flatb < 1) = flatm(flatb < 1);
    flats = bwlabel(flatb < 1, 4);
    for kk=unique(flats(flatb == 0)).'
       flatb(flats == kk) = 0;
    end
    flatbb = flatb > 0;
    flatb(flatb == -1) = NaN;

    % find boundary triangles
    du = false(size(flatb));
    du(2:2:2*ms(1), 2:2:2*ms(2)) = flatbb(2:2:2*ms(1), 2:2:2*ms(2)) ...
                                   & flatbb(1:2:2*ms(1)-1, 1:2:2*ms(2)-1) ...
                                   & flatbb(3:2:2*ms(1)+1, 3:2:2*ms(2)+1);
    dd = false(size(flatb));
    dd(2:2:2*ms(1), 2:2:2*ms(2)) = flatbb(2:2:2*ms(1), 2:2:2*ms(2)) ...
                                   & flatbb(1:2:2*ms(1)-1, 3:2:2*ms(2)+1) ...
                                   & flatbb(3:2:2*ms(1)+1, 1:2:2*ms(2)-1);
    dua = du;
    dua(2:2:2*ms(1), 2:2:2*ms(2)) = dua(2:2:2*ms(1), 2:2:2*ms(2)) ...
                                    & ~isnan(flatb(1:2:2*ms(1)-1, 3:2:2*ms(2)+1));
    dub = du;
    dub(2:2:2*ms(1), 2:2:2*ms(2)) = dub(2:2:2*ms(1), 2:2:2*ms(2)) ...
                                    & ~isnan(flatb(3:2:2*ms(1)+1, 1:2:2*ms(2)-1));
    dda = dd;
    dda(2:2:2*ms(1), 2:2:2*ms(2)) = dda(2:2:2*ms(1), 2:2:2*ms(2)) ...
                                    & ~isnan(flatb(3:2:2*ms(1)+1, 3:2:2*ms(2)+1));
    ddb = dd;
    ddb(2:2:2*ms(1), 2:2:2*ms(2)) = ddb(2:2:2*ms(1), 2:2:2*ms(2)) ...
                                    & ~isnan(flatb(1:2:2*ms(1)-1, 1:2:2*ms(2)-1));

    % finally build mesh
    vs = [flatx(~isnan(flatb)).', 0.8*flatx(dua).'+0.2*flatx(dua([2:end, 1],[end, 1:end-1])).', 0.8*flatx(dub).'+0.2*flatx(dub([end, 1:end-1],[2:end, 1])).', 0.8*flatx(dda).'+0.2*flatx(dda([end, 1:end-1],[end, 1:end-1])).', 0.8*flatx(ddb).'+0.2*flatx(ddb([2:end, 1],[2:end, 1])).'; ...
          flaty(~isnan(flatb)).', 0.8*flaty(dua).'+0.2*flaty(dua([2:end, 1],[end, 1:end-1])).', 0.8*flaty(dub).'+0.2*flaty(dub([end, 1:end-1],[2:end, 1])).', 0.8*flaty(dda).'+0.2*flaty(dda([end, 1:end-1],[end, 1:end-1])).', 0.8*flaty(ddb).'+0.2*flaty(ddb([2:end, 1],[2:end, 1])).'];
    bs = [flatb(~isnan(flatb)).', flatb(dua([2:end, 1],[end, 1:end-1])).', flatb(dub([end, 1:end-1],[2:end, 1])).', flatb(dda([end, 1:end-1],[end, 1:end-1])).', flatb(ddb([2:end, 1],[2:end, 1])).'];
    ds = [flatd(~isnan(flatb)).', .5+0*flatx(dua).', .5+0*flatx(dub).', .5+0*flatx(dda).', .5+0*flatx(ddb).'];
    flatv = zeros(2*ms+1);
    flatv(~isnan(flatb)) = 1:sum(~isnan(flatb), 'all');
    flatvua = zeros(2*ms+1);
    flatvua(dua) = (1:sum(dua, 'all')) + sum(~isnan(flatb), 'all');
    flatvub = zeros(2*ms+1);
    flatvub(dub) = (1:sum(dub, 'all')) + sum(dua, 'all') + sum(~isnan(flatb), 'all');
    flatvda = zeros(2*ms+1);
    flatvda(dda) = (1:sum(dda, 'all')) + sum(dub, 'all') + sum(dua, 'all') +sum(~isnan(flatb), 'all');
    flatvdb = zeros(2*ms+1);
    flatvdb(ddb) = (1:sum(ddb, 'all')) + sum(dda, 'all') + sum(dub, 'all') + sum(dua, 'all') +sum(~isnan(flatb), 'all');
    flatvm = flatv;
    flatvm(du | dd) = 0;
    cs = [reshape(flatvm(1:2*ms(1), 1:2*ms(2)), 1, []); ...
          reshape(flatvm(2:2*ms(1)+1, 1:2*ms(2)), 1, []); ...
          reshape(flatvm(2:2*ms(1)+1, 2:2*ms(2)+1), 1, []); ...
          reshape(flatvm(1:2*ms(1), 2:2*ms(2)+1), 1, [])];
    ele = min(cs, [], 1) > 0;
    csua = [flatv(dua([2:end, 1], [2:end, 1])).', flatv(dua).', flatv(dua([2:end, 1], :)).'; ...
            flatv(dua).', flatv(dua([end, 1:end-1], [end, 1:end-1])).', flatvua(dua).'; ...
            flatvua(dua).', flatv(dua(:, [end, 1:end-1])).', flatv(dua(:, [end, 1:end-1])).'; ...
            flatv(dua([2:end, 1], :)).', flatvua(dua).', flatv(dua([2:end, 1], [end, 1:end-1])).'];
    csub = [flatv(dub([2:end, 1], [2:end, 1])).', flatvub(dub).', flatv(dub(:, [2:end, 1])).'; ...
            flatv(dub(:, [2:end, 1])).', flatv(dub([end, 1:end-1], :)).', flatv(dub([end, 1:end-1], [2:end, 1])).'; ...
            flatvub(dub).', flatv(dub([end, 1:end-1], [end, 1:end-1])).', flatv(dub([end, 1:end-1], :)).'; ...
            flatv(dub).', flatv(dub).', flatvub(dub).'];
    csda = [flatv(dda([end, 1:end-1], [2:end, 1])).', flatvda(dda).', flatvda(dda).'; ...
            flatv(dda([end, 1:end-1], :)).', flatv(dda(:, [end, 1:end-1])).', flatv(dda([end, 1:end-1], :)).'; ...
            flatvda(dda).', flatv(dda([2:end, 1], [end, 1:end-1])).', flatv(dda([end, 1:end-1], [end, 1:end-1])).',; ...
            flatv(dda).', flatv(dda).', flatv(dda(:, [end, 1:end-1])).'];
    csdb = [flatv(ddb(:, [2:end, 1])).', flatv(ddb([2:end, 1], :)).', flatv(ddb([2:end, 1], [2:end, 1])).'; ...
            flatv(ddb([end, 1:end-1], [2:end, 1])).', flatvdb(ddb).', flatv(ddb(:, [2:end, 1])).'; ...
            flatv(ddb).', flatv(ddb).', flatvdb(ddb).'; ...
            flatvdb(ddb).', flatv(ddb([2:end, 1], [end, 1:end-1])).', flatv(ddb([2:end, 1], :)).'];
    qm.ns = ns;
    qm.vs = vs;
    qm.bs = bs;
    qm.ds = ds;
    qm.cs = [cs(:, ele), csua, csub, csda, csdb];
end