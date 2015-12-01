fpr = @(T) sutherlands(T, 1.716*10^-5, 110.6,273.1)*1005.4/sutherlands(T,0.0241, 194,273.1);

flr = @(Pr) Pr^(1/2);

ftr = @(Pr) Pr^(1/3);


Pr0 = 0.715;
fsst = state_init('Pressure', 1171.87, 'Density', 1.80119*10^-2, 'Temperature', 226.65, 'Mach', 10);
% Laminar Case, Upper surface
% Results from P-M expansion from tables
Me = 12.23;
Te = fsst.T*(1 + .2*fsst.M^2)/(1 + .2*Me^2);
Pe = fsst.p*((1 + .2*fsst.M^2)/(1 + .2*Me^2))^(1.4/.4);
est = state_init('Mach', Me, 'Temperature', Te, 'Pressure', Pe, 'Velocity', Me*sqrt(1.4*287*Te));
wst = state_init('Temperature', 300);
r0 = flr(Pr0);
Taw = Taw_function(est, r0);

% Reference Temperature State
Tref = est.T + 0.5 * (wst.T - est.T) + 0.22 * (Taw - est.T);
refst = state_init('Temperature', Tref, 'Pressure', Pe, 'Density', Pe/(287*Tref));

Reref = refst.r*est.V*1/sutherlands(refst.T, 1.716*10^-5, 110.6,273.1);
Prref = fpr(refst.T);

twall = .5*refst.r*est.V^2*0.664/sqrt(Reref)
qdotwall = refst.r * est.V * est.cp * (Taw - wst.T) * 0.332 / sqrt(Reref) * Prref^(-2/3)/(1*10^4)

% Lower surface using turbulant

Me =  8.28609988;
Te = fsst.T*1.42547896;
Pe = fsst.p*3.02422697;
est = state_init('Mach', Me, 'Temperature', Te, 'Pressure', Pe, 'Velocity', Me*sqrt(1.4*287*Te));
wst = state_init('Temperature', 300);
r0 = ftr(Pr0);
Taw = Taw_function(est, r0);

% Reference Temperature State
Tref = est.T + 0.5 * (wst.T - est.T) + 0.22 * (Taw - est.T);
refst = state_init('Temperature', Tref, 'Pressure', Pe, 'Density', Pe/(287*Tref));

Reref = refst.r*est.V*1/sutherlands(refst.T, 1.716*10^-5, 110.6,273.1);
Prref = fpr(refst.T);

cfx = 0.455/(log(Reref)/log(10))^(2.52);
twall = .5*refst.r*est.V^2*cfx
qdotwall = .5*cfx*(Prref) ^ (-2/3) * refst.r * est.V * est.cp * (Taw - wst.T)/(1*10^4)