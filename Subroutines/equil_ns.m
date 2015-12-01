function [p21, r21, h21, T21, dstate] = equil_ns(ustate, varargin)
% Normal Shock w/ Reacting Flow
err = 0;
if(ustate.u_flag)
    % imp2si is not implemented, don't be stupid and use imperial units!
    ustate = imp2si(ustate);
end

opt = max(size(varargin));
scplot = 0;
if opt
    cell2mat(varargin(opt));
    for i = 1:1:opt
        switch(cell2mat(varargin(i)))
            case 'ShowConverge'
                scplot = 1;
            otherwise
                error('equil_ns: Invalid option\n')
        end
    end
end

r1 = ustate.r;
v1 = ustate.V;
h1 = ustate.h;
p1 = ustate.p;
T1 = ustate.T;

eps0 = 0.001;
eps1 = .1;

err = 1e-6;

% Old values
r2 = r1 / eps0;
u2 = eps0 * v1;
p2 = p1 + r1 * v1 * (1 - eps0);
h2 = h1 + v1 ^ 2 / 2 * (1 - eps0 ^ 2);
[p, s, rho, e, a, h2_, T] = tgasM(2,p2,r2);

fepso = h2 - h2_;


delta = 1;
if(scplot)
    it = 0;
    figure
    hold on
end
it = 0;

while err < delta
	r2 = r1 / eps1;
	u2 = eps1 * v1;
	p2 = p1 + r1 * v1 ^ 2 * (1 - eps1);
	h2 = h1 + v1 ^ 2 / 2 * (1 - eps1 ^ 2);

	[p, s, rho, e, a, h2_, T] = tgasM(2,p2,r2);

	feps = h2 - h2_;

	tmp = eps1;

	eps1 = eps1 - feps/((feps - fepso)/(eps1-eps0));

	eps0 = tmp;
	fepso = feps;
	delta = abs(feps);
    if(scplot)
        plot(it,delta)
        pause(0.001)
    end
    it = it + 1;
    if(it>100 || imag(fepso))
        error('equil_ns: Could not converge to a valid solution\n')
        err = 1;
        break;
    end
end
T2 = T;
p21 = p2/p1;
T21 = T2/T1;
r21 = r2/r1;
h21 = h2/h1;

dstate = state_init();
dstate.p = p;
dstate.s = s;
dstate.r = rho;
dstate.e = e;
dstate.a = a;
dstate.h = h2_;
dstate.T = T;
dstate.V = u2;
dstate.M = u2/a;