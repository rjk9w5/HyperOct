function [downstream_state, bd] = cpg_obs(upstream_state, delta, gmma)

p2p1 = @(M,B) 1 + 2 * gmma / (gmma + 1) * (M ^ 2 * sin(B)^2 - 1);

T2T1 = @(M,B) (2 * gmma * M ^ 2 * sin(B)^2 - (gmma - 1))*((gmma - 1) * M ^ 2 * sin(B) ^ 2  + 2 ) ...
	/ ((gmma + 1) ^ 2 * M ^ 2 * sin(B) ^ 2);

r2r1 = @(M,B) ((gmma + 1) * M ^ 2 * sin(B) ^ 2) / ((gmma - 1) * M ^ 2 * sin(B) ^ 2 + 2);

u2V1 = @(M,B) 1 - 2*(M^2*sin(B)^2 - 1) / ((gmma+1)*M^2);

v2V1 = @(M,B) 2*(M^2*sin(B)^2 - 1)*cos(B) / ((gmma+1)*M^2);

T1 = upstream_state.T;
p1 = upstream_state.p;
h1 = upstream_state.h;
r1 = upstream_state.r;
M1 = upstream_state.M;
v1 = upstream_state.v;
a1 = upstream_state.a;

B = T_B_M(delta,M1,gmma);
bd = B-delta;

downstream_state.T = T1*T2T1(M1,B); T2T1(M1,B);

downstream_state.r = r1*r2r1(M1,B); r2r1(M1,B);

downstream_state.v = sqrt(((r1/downstream_state.r)*v1*sin(B))^2 + (v1*cos(B))^2);

downstream_state.a = sqrt(gmma*286.9*downstream_state.T);

downstream_state.p = p1*p2p1(M1,B); p2p1(M1,B);

downstream_state.h = 1004.5*downstream_state.T;

downstream_state.M = downstream_state.v/downstream_state.a;
