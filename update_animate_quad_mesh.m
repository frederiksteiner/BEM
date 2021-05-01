function update_animate_quad_mesh(p, cs)
    p.FaceVertexCData = cs.';
    drawnow;
end