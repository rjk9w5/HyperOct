function [p21, r21, h21, T21, downstream_state, bd] = equil_obs(upstream_state, delta, varargin)
err = 0;
if(upstream_state.u_flag)
    % imp2si is not implemented, don't be stupid and use imperial units!
    upstream_state = imp2si(ustate);
end

opt = max(size(varargin))
scplot = 0;
if opt
    cell2mat(varargin(opt))
    for i = 1:1:opt
        switch(cell2mat(varargin(i)))
            case "ShowConverge"
                scplot = 1;
            otherwise
                fprintf("Error: Invalid option\n")
        end
    end
end

% ------------------------------------------------------------------------------
% Inputs
% upstream_state:
v1 = upstream_state.V;
r1 = upstream_state.r;
h1 = upstream_state.h;
p1 = upstream_state.p;
T1 = upstream_state.T;


eps0 = .001;
eps1 = .1;

err = 1e-7;

% First iteration, no eps update
r2 = r1 / eps0;
ee = 1-eps0;
B = atan((ee - sqrt(ee ^ 2 - 4 * eps0 * tan(delta) ^ 2))/(2 * eps0 * ...
	tan(delta)));
u1n = v1*sin(B);
u2n = eps0 * u1n;
p2 = p1 + r1 * u1n ^ 2 * (1 - eps0);
h2 = h1 + u1n ^ 2 / 2 * (1 - eps0 ^ 2);
h2_ = tgas4(p2,r2);

fepso = h2 - h2_;

del = 1;
it = 0;
if(scplot)
    figure
    hold on
end

while err < del
	r2 = r1 / eps1;
	ee = 1 - eps1;
	B = atan((ee - sqrt(ee ^ 2 - 4 * eps1 * tan(delta) ^ 2))/(2 * eps1 * ...
		tan(delta)));
	u1n = v1 * sin(B);
	u2n = eps1 * u1n;
	p2 = p1 + r1 * u1n ^ 2 * (1 - eps1);
	h2 = h1 + u1n ^ 2 / 2 * (1 - eps1 ^ 2);

	[p, s, rho, e, a, h2_, T] = tgasM(2,p2,r2);

	feps = h2 - h2_;

	tmp = eps1;

	eps1 = eps1 - feps/((feps - fepso)/(eps1-eps0));

	eps0 = tmp;
	fepso = feps;
	del = abs(feps);
    
    if(scplot)
        plot(it, del)
    end
    it = it+1;
    if (it > 100 || imag(fepso))
        printf("Error: Could not converge to a valid solution\n")
        err = 1;
        break;
    end
end

T2 = T;
p21 = p2/p1;
T21 = T2/T1;
r21 = r2/r1;
h21 = h2/h1;

bd = B-delta;

downstream_state.v = sqrt((v1*cos(B))^2 + u2n^2);
downstream_state.p = p2;
downstream_state.r = r2;
downstream_state.a = a;
downstream_state.h = h2;
downstream_state.T = T; 
downstream_state.M = downstream_state.v/downstream_state.a;
