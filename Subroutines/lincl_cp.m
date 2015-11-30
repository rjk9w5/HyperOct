function cp = lincl_cp(method, theta, state, gam)
% ----------------------------------------------------------------------------- 
% ----------------------------------------------------------------------------- 
% Inputs:
% 	method:	 
%		1 - CPG  
%		2 - HSP
%		3 - Newtionian Method
%       4 - Modified Newtonian Method w/ CPG
%       5 - Modified Newtonian Method w/ Equillibrium
%		6 - Small angles and HSL
%		7 - HSL
% 		8 - Equillibrium
%       9 - Tangent Wedge Method
%       10- Tangent Cone Method [Not Implemented]
% 	theta:	deflection angle in radians
%	gam:	specific heat ratio of free stream fluid
%	M1:		Free stream Mach number
% 	varargin:
% 		for CPG inputs are []
%		for HPS inputs are []
%		for Equillibrium inputs are []
% ----------------------------------------------------------------------------- 
% ----------------------------------------------------------------------------- 
% Outputs:
%	cp:		Pressure coefficeint approximation
% ----------------------------------------------------------------------------- 
% -----------------------------------------------------------------------------

	if (nargin < 4)
		gam = 1.4;
	end
    if(method ~= 3 && method ~= 6)
        M1 = state.M;
        fcp = @(p21) 2 / (gam * M1 ^ 2) * (p21 - 1);
        if(method ~= 4 && method ~= 5)
            B = T_B_M(theta,M1);
        end
    end
	
	switch(method)
	case 1	
		us_state = cpg_obs(state, theta, gam);
		prat = us_state.p/state.p;		
		cp = fcp(prat);
	case 2
		K = M1*theta;
		prat = 1 + (gam * (gam + 1)) / 4 * K ^ 2 + gam * K ^ 2 * sqrt(((gam + 1) / 4) ^ 2 + 1 / K ^ 2);
		cp = fcp(prat);
	case 3
		cp = 2 * sin(theta) ^ 2;
    case 4
        [p2p1,~,~,~] = ns(state, "CPG");
        cpmax = fcp(p2p1);
        cp = cpmax*sin(theta) ^ 2;
    case 5
        [p2p1,~,~,~] = ns(state, "Equillibrium");
        cpmax = fcp(p2p1);
        cp = cpmax*sin(theta) ^ 2;
	case 6
		cp = (gam + 1) * theta ^ 2;
	case 7
		cp = 4/(gam + 1)*sin(B) ^ 2;
	case 8
		us_state = equil_obs(state, theta);
		prat = us_state.p/state.p;
		cp = fcp(prat);
    case 9
        [prat, ~, ~, ~, ~] = obs(state, theta, "CPG");
        cp = fcp(prat);
	otherwise
		disp('Not a valid entry for method parameter!\n')
		disp('Please input a value from the follwing list:\n\t1: CPG\n\t2: HSP\n\t3: Equillibrium');
		method = input('Method: ')
	end
	
end
