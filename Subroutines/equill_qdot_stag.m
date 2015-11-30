function [out] = equill_qdot_stag(fs_state, edge_state, wall_state, Rn)

	pe = edge_state.p;
	re = edge_state.r;
	Te = edge_state.T;
	mue = sutherlands(Te);

	Tw = wall_state.T;
	rw = wall_state.r;
	hw = wall_state.h;
	muw = sutherlands(Tw);

	p1 = fs_state.p;

	due = (1/Rn)*sqrt((2*(pe-p1))/re);

	haw = Taw_function(edge_state, sqrt(.715));
	hd=0;
	h0e=1;
	Le = 1.4;
	out = 0.76*(0.715)^(0.6)*(re*mue)^(0.5)*(rw*muw)^(0.1)*sqrt(due)*edge_state.cp*(haw-hw)*(1 + (Le^(0.52) - 1)*(hd/h0e));