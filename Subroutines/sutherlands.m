function out = sutherlands(T, mref, S, Tref)
%	mref = 1.786*10^-5;
% 	Tref = 288;
% 	S = 110;

	out = mref*(T/Tref)^(3/2)*(Tref + S) / (T + S);