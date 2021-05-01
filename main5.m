%%% to test out nystroem cauchy

nop = 300:200:1100;

[mbp, hmax] = generate_inner_test_parametrisation('neutraly');
d = [1;1];
kappa = 1 + 2i;
h = hmax;
n = length(nop);
err= zeros(n,1);
qm = generate_quad_mesh(mbp, h);
bool = [true, false,true, false];
u = @(x) 1i/4*besselh(0, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2));
dnu = @(x, nx) -1i/4*kappa*besselh(1, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2))...
    .*((x(1,:)-d(1)).*nx(1,:) + (x(2,:)-d(2)).*nx(2,:))...
    ./sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2);

for i = 1:n
    fprintf('Start of mesh with %d points \n', nop(i));
    dcd2 = nystroem_cauchy(mbp, nop(i)*ones(size(mbp)), u, dnu,kappa, bool);
    [zs2, err(i)] = cauchy_data_onto_quad_mesh(qm, dcd2, kappa, u);
    
end

figure()
semilogy(nop, err)
semilogy(nop, err, '-x')
title('Konvergenzplot nystroem cauchy')
xlabel('Anzahl Punkte im Mesh')
ylabel('Ausgerechneter Fehler')

disp(['Error is: ' num2str(err(n))]);