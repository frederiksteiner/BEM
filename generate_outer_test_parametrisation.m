function [mbp, hmax] = generate_outer_test_parametrisation(name)
    if nargin < 1
        name = 'smirky';
    end
    switch name
        case 'kite'
            mbp = struct([]);
            no = -1;
            vs = [.175, .325, -.825, .325;
                     0,  .75,     0, -.75];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.11;
        case 'smirky'
            mbp = struct([]);
            no = -1;
            vs = [.47,  .3, .13,   .3;
                  .15, .32, .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [-.13, -.3, -.47,  -.3;
                   .15, .32,  .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [  .45,  .30,    0, -.30, -.45, -.30,    0,  .30;
                   -.22, -.27, -.30, -.32, -.48, -.44, -.50, -.40];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.0148;
        case 'smiley'
            mbp = struct([]);
            no = -1;
            vs = [.47,  .3, .13,   .3;
                  .15, .32, .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [-.13, -.3, -.47,  -.3;
                   .15, .32,  .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [ .45,   0, -.45,   0;
                  -.3, -.35, -.3, -.55];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.03;
        case 'happy'
            mbp = struct([]);
            no = -1;
            vs = [.47,  .3, .13,   .3;
                  .15, .32, .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [-.13, -.3, -.47,  -.3;
                   .15, .32,  .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [ .45,   0, -.45,   0;
                  -.35, -.3, -.35, -.5];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.03;
        case 'neutraly'
            mbp = struct([]);
            no = -1;
            vs = [.47,  .3, .13,   .3;
                  .15, .32, .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [-.13, -.3, -.47,  -.3;
                   .15, .32,  .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [.45,   0, -.45,   0;
                  -.4, -.3,  -.4, -.5];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.03;
        case 'saddy'
            mbp = struct([]);
            no = -1;
            vs = [.47,  .3, .13,   .3;
                  .15, .32, .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [-.13, -.3, -.47,  -.3;
                   .15, .32,  .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [ .45,   0, -.45,   0;
                  -.45, -.3, -.45, -.5];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.03;
        case 'criey'
            mbp = struct([]);
            no = -1;
            vs = [.47,  .3, .13,   .3;
                  .15, .32, .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [-.13, -.3, -.47,  -.3;
                   .15, .32,  .15, -.02];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            vs = [ .45,    0, -.45,    0;
                   -.5, -.25,  -.5, -.45];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.03;
        case 'times'
            mbp = struct([]);
            no = -1;
            vs = [ .2,  .3,  .4,  .3,  .2,  .1,  .0, -.1, -.2, -.3, -.4, -.3, -.2, -.3, -.4, -.3, -.2, -.1,  .0,  .1,  .2,  .3,  .4,  .3;
                   .0,  .1,  .2,  .3,  .4,  .3,  .2,  .3,  .4,  .3,  .2,  .1,  .0, -.1, -.2, -.3, -.4, -.3, -.2, -.3, -.4, -.3, -.2, -.1];
            mbp = add_fourier_shape_parametrisation(mbp, vs, no);
            hmax = 0.03;
        otherwise
            mbp = struct([]);
            no = -1;
            cvs = [.1, 0, -.1, 0; 0, .1, 0, -.1];
            step = [0; 0.4];
            for jj=-3:3
                mbp = add_fourier_shape_parametrisation(mbp, jj*step+cvs, no);
            end
            hmax = 0.03;
    end
end