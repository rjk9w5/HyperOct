function [out] = equill_qdot_stag(fs_state, edge_state, wall_state, Rn)

	pe = edge_state.p;
	re = edge_state.r;
	Te = edge_state.T;
	mue = sutherlands(Te,1.786*10^-5, 110,288);

	Tw = wall_state.T;
	rw = wall_state.r;
	hw = wall_state.h
	muw = sutherlands(Tw,1.786*10^-5, 110,288);

	p1 = fs_state.p;

	due = (1/Rn)*sqrt((2*(pe-p1))/re);

	h0e = edge_state.h + edge_state.V^2/2
	hd=0;
	Le = 1.4;
	out = 0.763*(0.715)^(-0.6)*(re*mue)^(0.4)*(rw*muw)^(0.1)*sqrt(due)*(h0e-hw)*(1 + (Le^(0.52) - 1)*(hd/h0e));

	%out = 0.763*(0.715)^(-0.6)*(re*mue)^(0.5)*sqrt(due)*edge_state.cp*(Taw-Tw);