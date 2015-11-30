function [p2p1, T2T1, r2r1, bd, downstream, B] = obs(upstream, flwdefl, method)

downstream = state_init("Units",upstream.u_flag);
bd=0;


if(nargin > 1)
    switch method
        case "CPG" % Calorically Perfect Gas
            M = upstream.M;
            gmma = upstream.gma;
            B = T_B_M(flwdefl,M,gmma);
            bd = B-flwdefl;
            
            p2p1 = 1 + 2 * gmma / (gmma + 1) * (M ^ 2 * sin(B)^2 - 1);

            T2T1 = (2 * gmma * M ^ 2 * sin(B)^2 - (gmma - 1)) * ... 
                ((gmma - 1) * M ^ 2 * sin(B) ^ 2  + 2 ) / ((gmma + 1) ...
                ^ 2 * M ^ 2 * sin(B) ^ 2);

            r2r1 = ((gmma + 1) * M ^ 2 * sin(B) ^ 2) / ((gmma - 1) * ...
                M ^ 2 * sin(B) ^ 2 + 2);
        case "Equillibrium" % Thermo-chemical Equillibrium
            [p2p1, r2r1, ~, T2T1, downstream_state, bd] = equil_obs(upstream_state, flwdefl)
        case "HSP" % Hypersonic Similarity Parameter
            fprintf("Error: Option not available")
        case "HL" % Hypersonic Limit
            fprintf("Error: Option not available")
        case "HL+SA" % Hypersonic Limit + Small Angles 
            fprintf("Error: Option not available")
        otherwise
            fprintf("Error: Invalid option \n")
            fprintf("Valid options are: \n")
            fprintf("1.\t\"CPG\"\n")
            fprintf("2.\t\"Equillibrium\"\n")
            fprintf("3.\t\HSP\"\n")
            fprintf("4.\t\"HSL\"\n")
            fprintf("5.\t\"HSL+S\"")
    end %switch
else
    fprintf("Error: No option specified\n")
    fprintf("Exiting...\n")
end %if


%{
state.T = -1;     % temperature [K]
state.p = -1;     % pressure [Pa]
state.r = -1;     % density [kg/m^3]
state.V = -1;     % velocity magnitude [m/s]
state.V_(1:3) = [0 0 0]; % velocity vector [m/s]
state.M = -1;     % Mach number
state.gma = -1;   % specific heat ratio
state.cp = -1;    % specific heat -> constant pressure [J/K]
state.cv = -1;    % specific heat -> constant volume [J/K]
state.e = -1;     % energy per unit mass [J/kg]
state.h = -1;     % enthalpy per unit mass [J/kg]
state.s = -1;     % entropy per unit mass [J/(K-kg)]
state.a = -1;     % speed of sound [m/s]
%}
