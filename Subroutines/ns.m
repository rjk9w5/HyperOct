function [p2p1, T2T1, r2r1, M, dstate] = ns(ustate, method)

% Shock angle is 90 degrees for a normal shock
B = pi/2;

dstate = state_init('Units', ustate.u_flag);

if(nargin > 1)
    switch method
        case 'CPG' % Calorically Perfect Gas
            M = ustate.M;
            gmma = ustate.gma;
            p2p1 = 1 + 2 * gmma / (gmma + 1) * (M ^ 2 * sin(B)^2 - 1);

            T2T1 = (2 * gmma * M ^ 2 * sin(B)^2 - (gmma - 1)) * ... 
                ((gmma - 1) * M ^ 2 * sin(B) ^ 2  + 2 ) / ((gmma + 1) ...
                ^ 2 * M ^ 2 * sin(B) ^ 2);

            r2r1 = ((gmma + 1) * M ^ 2 * sin(B) ^ 2) / ((gmma - 1) * ...
                M ^ 2 * sin(B) ^ 2 + 2);
            M = ((gmma - 1)*M ^2 + 2) / ((2*gmma*M^2) - (gmma - 1));
        case 'Equillibrium' % Thermo-Chemical Equillibrium
            [p2p1, r2r1, ~, T2T1, dstate] = equil_ns(ustate);
        case 'HSP' % Hypersonic Similarity Parameter
            fprintf('Error: Option not available')
        case 'HL' % Hypersonic Limit
            fprintf('Error: Option not available')
        case 'HL+SA' % Hypersonic Limit + Small Angles
            fprintf('Error: Option not available')
        otherwise
            fprintf('Error: Invalid option \n')
            fprintf('Valid options are: \n')
            fprintf('1.\t"CPG"\n')
            fprintf('2.\t"Equillibrium"\n')
            fprintf('3.\t"HSP"\n')
            fprintf('4.\t"HSL\"\n')
            fprintf('5.\t"HSL+S"')
    end %switch
else
    fprintf('ns: No option specified\n')
    fprintf('Exiting...\n')
end %if
