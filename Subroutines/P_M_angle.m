function [nu] = P_M_angle(varargin)

if isstruct(varargin(1))
    state = cell2mat(varargin(1));
    M = state.M
    gam = state.gma
else
    M = cell2mat(varargin(1));
    if(nargin == 2)
        gam = cell2mat(varargin(2));
    else
        gam = 1.4;
    end
end

a = sqrt((gam + 1) / (gam - 1));

b = atan(sqrt(((gam - 1) / (gam + 1)) * (M ^ 2 - 1)));

c = atan(sqrt(M ^ 2 - 1));

nu = a * b - c;
