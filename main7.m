%%% main 7
clear all
[mbp, hmax] = generate_inner_test_parametrisation('neutraly');
nop = 1000;
n = size(nop, 2);
omega = 2*pi; gamma = 0.2; c = 1;
kappa = omega*(omega + 1i*gamma)/c^2;
h = hmax;
eta = real(kappa)/2;
bool = [true, true,true, true]; 

d=[-0.6;0.5];
f = @(x) zeros(1,size(x,2));
ui = @(x) (1i/4*besselh(0, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2)));
dnui = @(x, nx) (-1i/4*kappa*besselh(1, kappa*sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2))...
    .*((x(1,:)-d(1)).*nx(1,:) + (x(2,:)-d(2)).*nx(2,:))...
    ./sqrt((x(1,:)-d(1)).^2 + (x(2,:) - d(2)).^2));
us = @(x) f(x) - ui(x);
dnus = @(x, nx) - dnui(x, nx);
% ui = @(x) exp(1i*kappa*(x(1,:)*d(1) + x(2,:)*d(2)));
    pause('on'); clim = [-0.1, 0.1];
    qm = generate_quad_mesh(mbp, hmax);
    dcd = nystroem_cauchy(mbp, nop(1)*ones(size(mbp)), us, dnus,kappa, bool);
    [usv, err] = cauchy_data_onto_quad_mesh(qm, dcd, kappa, us);
    fprintf('Error is:  %d \n', err);
    uiv = ui(qm.vs);
    utot = uiv + usv;
    figure()
    subplot(1,3,1)
    ptot0 = setup_animate_quad_mesh(qm, real(utot), clim);
    title('U(x,0)')
    subplot(1,3,2) 
    pui0 = setup_animate_quad_mesh(qm, real(uiv), clim);
    title('U_i(x,0)')
    subplot(1,3,3)
    pus0 = setup_animate_quad_mesh(qm, real(usv), clim);
    title('U_s(x,0)')
    figure();
    subplot(1,3,1)
    ptot = setup_animate_quad_mesh(qm, real(utot), clim);
    subplot(1,3,2) 
    pui = setup_animate_quad_mesh(qm, real(uiv), clim);
    subplot(1,3,3)
    pus = setup_animate_quad_mesh(qm, real(usv), clim);
    
    
    
for t = 0.1:0.01:4
    fprintf('Start of timestep t =  %d \n', t);
    usvt = usv*exp(-1i*omega*t);
    uivt = uiv*exp(-1i*omega*t);
    utot = usvt + uivt;
    
    update_animate_quad_mesh(ptot, real(utot))
    update_animate_quad_mesh(pui, real(uivt))
    update_animate_quad_mesh(pus, real(usvt))
    
    pause(0.05)
    
    
    
    
end




