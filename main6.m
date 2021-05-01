%%% to test nystroem_cauchy_bw
[mbp, hmax] = generate_outer_test_parametrisation('neutraly');
nop = 50:100:650;

n = size(nop, 2);
d = [0; -0.35]; 
kappa = 2 + 1i;
err = zeros(1, n);
eta = real(kappa)/2;

u = @(x) 1i/4*besselh(0, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2));
qm = generate_quad_mesh(mbp, hmax);
for i = 1:n
    fprintf('Start of mesh with %d points \n', nop(i));
    dcd = nystroem_cauchy_bw(mbp, nop(i)*ones(size(mbp)), u, kappa, eta);
    [~, err(i)] = cauchy_data_onto_quad_mesh_bw(qm, dcd, kappa, u, eta);
end

semilogy(nop, err, '-x')
title('Konvergenzplot nystroem cauchy Brakhage-Werner')
xlabel('Anzahl Punkte im Mesh')
ylabel('Ausgerechneter Fehler')
