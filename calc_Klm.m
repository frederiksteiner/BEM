%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates part of stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Klm = calc_Klm(dcd, l, m, kappa)
nl = dcd(l).n;
nm = dcd(m).n; 
Klm =zeros(nl,nm);
ts = dcd(l).ts;
ss = dcd(m).ts;
gammasl = dcd(l).gammas;
gammasm = dcd(m).gammas;
ngammasm = dcd(m).ngammas;
normdgammasm = dcd(m).normdgammas;


if l ~= m
    
    for i = 0:nm-1
%         Klm(:,i+1) = 1i/(4*nm)*(gammasl - gammasm(:,i+1))'*ngammasm(:,i+1)*normdgammasm(i+1)./(2*pi*((gammasl(1,:) - gammasm(1,i+1)).^2 + (gammasl(2,:) - gammasm(2,i+1)).^2)');
        vecnorms = sqrt((gammasl(1,:) - gammasm(1,i+1)).^2 + (gammasl(2,:) - gammasm(2,i+1)).^2);
        Klm(:,i+1) = 1i*kappa/(4*nm)*besselh(1, kappa*vecnorms).*((gammasl - gammasm(:,i+1))'*ngammasm(:,i+1))'*normdgammasm(i+1)./(vecnorms);
        
    end
else
    d2gammas = dcd(l).d2gammas;
    W2 = toeplitz(dcd(l).wls);
    K1 = zeros(nm,nm);
    K2 = zeros(nm,nm);
    for i = 0:nm-1
        omegast = -log(sin(pi*(ts - ss(i+1))).^2);
        vecnorms = sqrt((gammasl(1,:) - gammasm(1,i+1)).^2 + (gammasl(2,:) - gammasm(2,i+1)).^2);
        besselj1 = besselj(1, kappa*vecnorms);

        K1(:,i+1) = (kappa/4*(1i * besselh(1, kappa*vecnorms) + 1/pi*besselj1.*-omegast)...
            .*((gammasl - gammasm(:,i+1)).'*ngammasm(:,i+1)).'*normdgammasm(i+1)./(vecnorms)).';
        K2(:,i+1) = (kappa/(4*pi)*besselj1.*((gammasl - gammasm(:,i+1)).'*ngammasm(:,i+1)).'*normdgammasm(i+1)./(vecnorms)).';
        K1(i+1,i+1) = 1/(4*pi)*d2gammas(:,i+1).'*ngammasm(:,i+1)/normdgammasm(i+1); 
        K2(i+1,i+1) = 0;
    end
    Klm = 1/nm*K1 + K2.*W2 + 1/2*eye(nm);
end
end
