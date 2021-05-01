function plot_quad_mesh(qm, zs, cs)
    if nargin < 2
        patch('Faces', qm.cs.', ...
              'Vertices', qm.vs.', ...
              'FaceColor', 'none');
        axis 'equal';
    else
        if nargin < 3
            cs = zs;
        end
        xl = [min(qm.vs(1, :)), max(qm.vs(1, :))] * [1.05 -0.05; -0.05, 1.05] + [-eps, eps];
        yl = [min(qm.vs(2, :)), max(qm.vs(2, :))] * [1.05 -0.05; -0.05, 1.05] + [-eps, eps];
        hold off;
        patch('Faces', qm.cs.', ...
              'Vertices', [qm.vs.', zs.'], ...
              'FaceVertexCData', cs.', ...
              'FaceColor', 'interp', ...
              'EdgeColor', 'none');
        xlim(xl);
        ylim(yl);
        d = daspect;
        d(1:2) = geomean(d(1:2));
        d(3) = 2.5 * d(3);
        daspect(d);
    end
end