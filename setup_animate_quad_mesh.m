function p = setup_animate_quad_mesh(qm, cs, clim)
    xl = [min(qm.vs(1, :)), max(qm.vs(1, :))] * [1.05 -0.05; -0.05, 1.05] + [-eps, eps];
    yl = [min(qm.vs(2, :)), max(qm.vs(2, :))] * [1.05 -0.05; -0.05, 1.05] + [-eps, eps];
    p = patch('Faces', qm.cs.', ...
              'Vertices', qm.vs.', ...
              'FaceVertexCData', cs.', ...
              'FaceColor', 'interp', ...
              'EdgeColor', 'none');
    xlim(xl);
    ylim(yl);
    d = daspect;
    d(1:2) = geomean(d(1:2));
    daspect(d);
    ax = gca;
    ax.CLim = clim;
end