function us = potential_cauchy(dcd, kappa, xs)
% calculates approx. values of us in xs


r = numel(dcd);
    us = 0;
    for i = 1:r
        nm = dcd (i).n;
        omega = 1/nm;
        gammas = dcd(i).gammas;
        normdgammas = dcd(i).normdgammas;
        ngammas = dcd(i).ngammas;
        rhoDs = dcd(i).rhoDs;
        rhoNs = dcd(i).rhoNs;

        fctrsl = (omega*1i/4)* normdgammas.*rhoNs;
        fctrdl = (omega*1i/4)* normdgammas.*rhoDs*kappa;
        for j = 1:nm
            vecnormsq = sqrt((xs(1,:) - gammas(1,j)).^2 + (xs(2,:) - gammas(2,j)).^2);
            us = us + fctrsl(j).* besselh(0,kappa*vecnormsq)...
                    - fctrdl(j) ...
                      .* ( (xs(1,:)-gammas(1,j)).*ngammas(1,j) + (xs(2,:)-gammas(2,j)).*ngammas(2,j) ) .* besselh(1, kappa*vecnormsq)./vecnormsq;
        end
    end



end