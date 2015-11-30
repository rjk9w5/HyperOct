function out = Taw_function(edge_state, r)
	Me = edge_state.M;
	Te = edge_state.T;
	gmma = edge_state.gma;

	out = Te*(1 + r*(gmma-1)/2*Me^2);