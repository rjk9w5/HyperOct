% Determine the FreeStream conditions
addpath('./Subroutines')
% flight conditions 
M = [7.828;7.831;7.938;7.938];
q_inf = [24.88;25.33;31.55;32.20]*1000;
h_alt = [34.48;34.31;33.05;32.89]*1000;

T = [235.594; 235.118; 231.590; 231.142];
P = [602.559; 617.613; 742.728; 760.484];
r = [0.00890991; 0.00915000; 0.0111724; 0.0114617];
u = M.*sqrt(1.4*287*T);
% check q_inf

% Find pressure from meassured q0 and M0, determine if Temp or density had largest error 
% by changing both and comparing to q_inf
p2 = 2/1.4.*q_inf./M.^2;
abs(p2-P)./p2;

r2 = p2./(287*T);
abs(.5.*r2.*u.^2 - q_inf)./q_inf;
T2 = p2./(287*r);
u2 = M.*sqrt(1.4*287*T2);
abs(.5.*(p2./(287*T2)).*u2.^2 - q_inf)./q_inf;

r3 = r2*0.9 + r*0.1;
abs(.5.*r3.*u.^2 - q_inf)./q_inf;
% For 
%   t1 -> T2*1.2
%   t2 -> T2*0.9
%   t3 -> T2*0.9 + T*0.1
%   t4 -> T2*0.8 + T*0.2
T3 = [T2(1)*1.2; T2(2)*0.9; T2(3)*0.9 + T(3)*0.1; T2(4)*0.8 + T(4)*0.2];
u2 = M.*sqrt(1.4*287*T3);
abs(.5.*(p2./(287*T3)).*u2.^2 - q_inf)./q_inf;
r2 = (p2./(287*T3));

