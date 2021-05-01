%%% to test potential_cauchy

[mbp, hmax] = generate_inner_test_parametrisation('neutraly');
d = [1;1];
kappa = 1 + 2i;
h = hmax;
eta = real(kappa)/2;
u = @(x) 1i/4*besselh(0, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2));
dnu = @(x, nx) -1i/4*kappa*besselh(1, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2))...
    .*((x(1,:)-d(1)).*nx(1,:) + (x(2,:)-d(2)).*nx(2,:))...
    ./sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2);


nop = 100:200:1100;
n = length(nop);
err= zeros(n,1);
qm = generate_quad_mesh(mbp, h);
for i = 1:n
    fprintf('Start of mesh with %d points \n', nop(i));
    dcd = interpolate_cauchy_data(mbp, nop(i)*ones(size(mbp)), u, dnu);
    [zs, err(i)] = cauchy_data_onto_quad_mesh(qm, dcd, kappa, u);
end
fs = u(qm.vs);
semilogy(nop, err, '-x')
title('Konvergenzplot Potential cauchy')
xlabel('Anzahl Punkte im Mesh')
ylabel('Ausgerechneter Fehler')

