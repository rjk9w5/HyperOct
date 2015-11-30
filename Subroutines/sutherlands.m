function r = sutherlands(T)
	mref = 1.786*10^-5;
	Tref = 288;
	S = 110;

	r = mref*(T/Tref)^(3/2)*(Tref + S) / (T + S);