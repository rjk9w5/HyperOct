function out = qdot_stag(fs_state, edge_state, Tw, Rn)

	pe = edge_state.p;
	re = edge_state.r;
	Te = edge_state.T;
	mue = sutherlands(Te,1.786*10^-5, 110,288);

	p1 = fs_state.p;

	due = (1/Rn)*sqrt((2*(pe-p1))/re);

	Taw = Taw_function(edge_state, sqrt(.715));

	out = 0.763*(0.715)^(-0.6)*(re*mue)^(0.5)*sqrt(due)*edge_state.cp*(Taw-Tw);