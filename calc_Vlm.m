function Vlm = calc_Vlm(dcd, l, m, kappa)
nl = dcd(l).n;
nm = dcd(m).n;
Vlm =zeros(nl,nm);
ts = dcd(l).ts;
ss = dcd(m).ts;
gammasl = dcd(l).gammas;
gammasm = dcd(m).gammas;
normdgammasm = dcd(m).normdgammas;


if l ~= m
    for i = 0:nm-1
        
        Vlm(:,i+1) = 1i/(4*nm)*besselh(0, kappa*sqrt((gammasl(1,:) - gammasm(1,i+1)).^2 + (gammasl(2,:) - gammasm(2,i+1)).^2))...
            .*normdgammasm(i+1);
    end
else
    W2 = toeplitz(dcd(l).wls);
    C = 0.57721566490153286060651209008240243104215933593992;
    V1 = zeros(nm,nm);
    V2 = zeros(nm,nm);
    for i = 0:nm-1
       vecnorm = sqrt((gammasl(1,:) - gammasm(1,i+1)).^2 + (gammasl(2,:) - gammasm(2,i+1)).^2);
       besselj0 = besselj(0, kappa*vecnorm);
       omegast = -log(sin(pi*(ts - ss(i+1))).^2); 
       V1(:,i+1) = ((1i/4*besselh(0, kappa*vecnorm) + 1/(4*pi)*besselj0.*-omegast)*normdgammasm(i+1)).';
       V2(:,i+1) = (1/(4*pi)*besselj0*normdgammasm(i+1)).';
       V1(i+1,i+1) = (1i/4 - 1/(2*pi)*(C + log(kappa*normdgammasm(i+1)/(2*pi))))*normdgammasm(i+1);
       
    end
    Vlm = 1/nm*V1 + W2.*V2;
end
end

