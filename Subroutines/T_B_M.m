% ----------------------------------------------------------------------------%
% Solve for the shock angle Beta from:
%
% tan(theta) = 2*cot(B)*((M^2*sin(B)^2-1) / (M^2*(gam+cos(2*B))+2))
%
% Inputs:
%		T 		-	deflection angle
%		M		- 	Free Stream Mach Number
%		gam 	-	specific heat ratio
%
% Outputs:
%		B 		-	shock angle
% ----------------------------------------------------------------------------%

function [Beta] = T_B_M(T,M,gam)

	if(nargin < 3)
		gam = 1.4;
	end
	func = @(B) 2*cot(B)*((M^2*sin(B)^2-1) / (M^2*(gam+cos(2*B))+2)) - tan(T);

	% Estimate beta for a given mach number and deflection angle
	B0 = (60/M+T)*pi/180;

	Beta = NewtonMethod(B0,func);
end
