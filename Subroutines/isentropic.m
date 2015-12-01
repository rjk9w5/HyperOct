function [r] = isentropic(M, gma, ratio)

TT = 1 + (gma - 1) / 2 * M ^ 2;

switch ratio
    case "Temperature"
        r = TT;
    case "Pressure"
        r = TT ^ (gma / (gma - 1));
    case "Density"
        r = TT ^ (1/(gma + 1));
    otherwise
        error("insentropic: invalid option")
end
