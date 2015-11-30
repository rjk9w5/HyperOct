% function [p, s, rho, e, a, h, T] = tgasM(group, value1, value2)
function [p, s, rho, e, a, h, T] = tgasM(group, value1, value2)
%
%--------------------------------------------------------------------------
%   Enter group = 1 for Energy, Density         |--> TGAS
%               = 2 for Pressure, Density       |--> TGAS
%               = 3 for Pressure, Entropy       |--> TGAS
%               = 4 for Enthalpy, Entropy       |--> airtbn
%               = 5 for Pressure, Temperature   |--> airtbn
%   (enter property values in the same order as written above)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%   Use Consistant SI units!!!!!!!!!!!
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%   Output  p   -> Pressure [Pa]
%           s   -> entropy [J/(kg*K)]
%           rho -> Density [kg/m^3]
%           e   -> energy [J/kg]
%           a   -> speed of sound [m/s]
%           h   -> enthalpy [J/kg]
%           T   -> temperature [K]
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%   Written by Aaron Erb
%   References:
%       "Simplified Curve Fits for the Thermodynamic Properties of 
%        Equilibrium Air"
%           S. Srinivasan, J.C. Tannehill, and K.J. Weilmuenster
%           NASA RP1181, August 1987
%       EQUAIR Applet
%           William Devenport
%           http://www.dept.aoe.vt.edu/~devenpor/tgas/
%--------------------------------------------------------------------------


load coefficients

%------------------------------
%   set sea level constants
%------------------------------
R = 287.06;
p0 = 101330;
s0 = 6779.2;
rho0 = 1.292;
e0 = 78408.4;
a0 = 331.3613;
h0 = 156873;
T0 = 273.15;

%--------------------------------------------------------------------------
%   Solve Initial Functions
%--------------------------------------------------------------------------
switch group
    %----------------------------------------------------------------------
    %   [Energy && Density]
    %----------------------------------------------------------------------
    case 1  
        e = value1;
        rho = value2;
        
        %-----------------------------------------------
        %   Pressure [p = p(e,rho)]
        %-----------------------------------------------
        p = pFUN(rho,rho0,e,T0,R,coefficients.A);
        
        %-----------------------------------------------
        %  Speed of Sound [a = a(e,rho)]
        %-----------------------------------------------
        a = aFUN(rho,rho0,e,T0,R,coefficients.A);
        
        %-----------------------------------------------
        %    Temperature [T = T(e,rho)]
        %-----------------------------------------------
        T = TFUN(rho,rho0,p,p0,T0,R,coefficients.B);
        
        %-----------------------------------------------
        %   Entropy [s = s(e,rho)]
        %-----------------------------------------------
        s = sFUN(rho,rho0,e,R,T0,coefficients.E);
        
        
    %----------------------------------------------------------------------
    %   [Pressure && Density]
    %----------------------------------------------------------------------        
    case 2
        p = value1;
        rho = value2;
        
        %-----------------------------------------------
        %   Enthalpy [h = h(p,rho)]
        %-----------------------------------------------
        h = hFUN(rho,rho0,p,p0,coefficients.C);
        
        %-----------------------------------------------
        %   Temperature [T = T(p,rho)]
        %-----------------------------------------------
        T = TFUN2(rho,rho0,p,p0,T0,R,coefficients.D);
        
       
    %----------------------------------------------------------------------
    %   [Pressure && Entropy]
    %----------------------------------------------------------------------        
    case 3
        p = value1;
        s = value2;
        
        %------------------------------------------------
        %   Density [rho = rho(p,s)]
        %------------------------------------------------
        rho = rhoFUN(s,R,p,p0,rho0,coefficients.F);
        
        %------------------------------------------------
        %   Energy [e = e(p,s)]
        %------------------------------------------------
        e = eFUN(s,R,p,p0,T0,coefficients.G);
        
        %------------------------------------------------
        %   Speed of Sound [a = a(p,s)]
        %------------------------------------------------
        a = aFUN2(s,R,p,p0,a0,coefficients.H);
        
        
    %----------------------------------------------------------------------
    %   [Enthalpy && Entropy]
    %----------------------------------------------------------------------
    case 4
        h = value1;
        s = value2;
        
        %------------------------------------------------------------------
        %   Iteratively Solves [p = p(h,s), e = e(h,s), rho = rho(h,s)]
        %------------------------------------------------------------------
        [p, e, rho] = perho(h,s,s0,p0,T0,rho0,h0,R,coefficients.G,coefficients.F);
        
    
    %----------------------------------------------------------------------
    %   [Pressure && Temperature]
    %----------------------------------------------------------------------    
    case 5
        p = value1;
        T = value2;
        
        %------------------------------------------------------------------
        %   Iteratively Solves [rho = rho(p,T)]
        %------------------------------------------------------------------
        rho = rhoFUN2(p,T,p0,T0,rho0,R,coefficients.D);
        
        
    otherwise
        fprintf('Choose a correct group...');
end



%--------------------------------------------------------------------------
%   Solve for Remaining Vaiables
%--------------------------------------------------------------------------
switch group
    case 1
        h = hFUN(rho,rho0,p,p0,coefficients.C);
        %s = sFUN(rho,rho0,e,R,T0,coefficients.E);
        
    case 2
        e = h - (p/rho);
        s = sFUN(rho,rho0,e,R,T0,coefficients.E);
        a = aFUN(rho,rho0,e,T0,R,coefficients.A);
        %a = aFUN2(s,R,p,p0,a0,coefficients.H);
        
    case 3
        T = TFUN(rho,rho0,p,p0,T0,R,coefficients.B);
        %T = TFUN2(rho,rho0,p,p0,T0,R,coefficients.D);
        h = hFUN(rho,rho0,p,p0,coefficients.C);
        
    case 4
        a = aFUN(rho,rho0,e,T0,R,coefficients.A);
        %a = aFUN2(s,R,p,p0,a0,coefficients.H);
        T = TFUN(rho,rho0,p,p0,T0,R,coefficients.B);
        %T = TFUN2(rho,rho0,p,p0,T0,R,coefficients.D);
        
    case 5
        h = hFUN(rho,rho0,p,p0,coefficients.C);
        e = h - (p/rho);
        s = sFUN(rho,rho0,e,R,T0,coefficients.E);
        a = aFUN(rho,rho0,e,T0,R,coefficients.A);
        %a = aFUN2(s,R,p,p0,a0,coefficients.H);
end
end



%--------------------------------------------------------------------------
%   Pressure [p = p(e,rho)]
%--------------------------------------------------------------------------
function p = pFUN(rho,rho0,e,T0,R,A)

Y = log10(rho/rho0);
Z = log10(e/(R*T0));

% mflag = 0;
% lflag = 0;
% kflag = 0;

% jflag = 0;
% Yhigh = 1;
% Ylow = 0;
% Ym = 0;
% rhosave = rho;
% phigh = 0;
% 
% if abs(Y + 4.5) < 0.025
%     iflag = 0;
%     rhosave = rho;
%     Ym = Y;
%     Y = -4.4750;
%     Yhigh = Y;
%     rho = 10^Y*rho0;
%     jflag = -1;
% elseif abs(Y + 0.5) < 0.005
%     iflag = 1;
%     rhosave = rho;
%     Ym = Y;
%     Y = -0.4950;
%     Yhigh = Y;
%     rho = 10^Y*rho0;
%     jflag = -1;
% else
%     iflag = -1;
% end

if -7.0 <= Y && Y <= -4.5
    layer = 1;
    if Z <= 0.65
        column = 1;
    end
    if 0.65 < Z && Z <= 1.50
        column = 2;
    end
    if 1.50 < Z && Z <= 2.20
        column = 3;
    end
    if 2.20 < Z && Z <= 3.05
        column = 4;
    end
    if 3.05 < Z && Z <= 3.40
        column = 5;
    end
    if 3.40 < Z
        column = 6;
    end
end
if -4.5 < Y && Y <= -0.5
    layer = 2;
    if Z <= 0.65
        column = 1;
    end
    if 0.65 < Z && Z <= 1.50
        column = 2;
    end
    if 1.50 < Z && Z <= 2.22
        column = 3;
    end
    if 2.22 < Z && Z <= 2.95
        column = 4;
    end
    if 2.95 < Z
        column = 5;
    end
end
% if -0.5 < Y && Y <= 3.0
if -0.5 < Y
    layer = 3;
    if Z <= 0.65
        column = 1;
    end
    if 0.65 < Z && Z <= 1.70
        column = 2;
    end
    if 1.70 < Z && Z <= 2.35
        column = 3;
    end
    if 2.35 < Z
        column = 4;
    end
end

a = A(:,column,layer);

% disp('a pressure')
% disp([column layer]);
% disp(a);

gamma = a(1) + a(2)*Y + a(3)*Z + a(4)*Y*Z + a(5)*Y^2 ....
        + a(6)*Z^2 + a(7)*Y^2*Z + a(8)*Y*Z^2 + a(9)*Y^3 ....
        + a(10)*Z^3 + (a(11) + a(12)*Y + a(13)*Z ....
        + a(14)*Y*Z + a(15)*Y^2 + a(16)*Z^2 + a(17)*Y^2*Z ....
        + a(18)*Y*Z^2 + a(19)*Y^3 + a(20)*Z^3)/(1 + a(25)*exp(a(21) ....
        + a(22)*Y + a(23)*Z + a(24)*Y*Z));

p = rho*e*(gamma - 1);

% if iflag == -1 || jflag == 1
%     return;
% end
% if jflag == 0
%     plow = p;
%     p = plow+(phigh-plow)/(Yhigh-Ylow)*(Ym-Ylow);
%     rho = rhosave;
% end  
end


%--------------------------------------------------------------------------
%   Speed of Sound [a = a(e,rho)]
%--------------------------------------------------------------------------
function a = aFUN(rho,rho0,e,T0,R,A)

Y = log10(rho/rho0);
Z = log10(e/(R*T0));

% jflag = 0;
% Yhigh = 1;
% Ylow = 0;
% Ym = 0;
% rhosave = rho;
% ahigh = 0;
% 
% if abs(Y + 4.5) < 0.025
%     iflag = 0;
%     jflag = -1;
%     rhosave = rho;
%     Ym = Y;
%     Y = -4.4750;
%     Yhigh = Y;
%     rho = 10^Y*rho0;
% elseif abs(Y + 0.5) < 0.005
%     iflag = 1;
%     jflag = -1;
%     rhosave = rho;
%     Ym = Y;
%     Y = -0.4950;
%     Yhigh = Y;
%     rho = 10^Y*rho0;
% else
%     iflag = -1;
% end

if -7.0 <= Y && Y <= -4.5
    layer = 1;
    if Z <= 0.65
        column = 1;
    end
    if 0.65 < Z && Z <= 1.50
        column = 2;
    end
    if 1.50 < Z && Z <= 2.20
        column = 3;
    end
    if 2.20 < Z && Z <= 3.05
        column = 4;
    end
    if 3.05 < Z && Z <= 3.40
        column = 5;
    end
    if 3.40 < Z
        column = 6;
    end
end
if -4.5 < Y && Y <= -0.5
    layer = 2;
    if Z <= 0.65
        column = 1;
    end
    if 0.65 < Z && Z <= 1.50
        column = 2;
    end
    if 1.50 < Z && Z <= 2.22
        column = 3;
    end
    if 2.22 < Z && Z <= 2.95
        column = 4;
    end
    if 2.95 < Z
        column = 5;
    end
end
%if -0.5 < Y && Y <= 3.0
if -0.5 < Y
    layer = 3;
    if Z <= 0.65
        column = 1;
    end
    if 0.65 < Z && Z <= 1.70
        column = 2;
    end
    if 1.70 < Z && Z <= 2.35
        column = 3;
    end
    if 2.35 < Z
        column = 4;
    end
end

a = A(:,column,layer);

% disp('a sound');
% disp([column layer]);
% disp(a);

gamma = a(1) + a(2)*Y + a(3)*Z + a(4)*Y*Z + a(5)*Y^2 ....
        + a(6)*Z^2 + a(7)*Y^2*Z + a(8)*Y*Z^2 + a(9)*Y^3 ....
        + a(10)*Z^3 + (a(11) + a(12)*Y + a(13)*Z ....
        + a(14)*Y*Z + a(15)*Y^2 + a(16)*Z^2 + a(17)*Y^2*Z ....
        + a(18)*Y*Z^2 + a(19)*Y^3 + a(20)*Z^3)/(1 + a(25)*exp(a(21) ....
        + a(22)*Y + a(23)*Z + a(24)*Y*Z));

delGAMMA_delY = a(2) + a(4)*Z + 2*a(5)*Y + 2*a(7)*Y*Z + a(8)*Z^2 + 3*a(9)*Y^2 ....
                + (a(12) + a(14)*Z + 2*a(15)*Y + 2*a(17)*Y*Z + a(18)*Z^2 + 3*a(19)*Y^2)/ ....
                (1 + a(25)*exp(a(21) + a(22)*Y + a(23)*Z + a(24)*Y*Z)) ....
                - a(25)*(a(11) + a(12)*Y + a(13)*Z + a(14)*Y*Z + a(15)*Y^2 + a(16)*Z^2 ....
                + a(17)*Y^2*Z + a(18)*Y*Z^2 + a(19)*Y^3 + a(20)*Z^3)*(a(22) + a(24)*Z)* ....
                (exp(a(21) + a(22)*Y + a(23)*Z + a(24)*Y*Z))/ ....
                (1 + a(25)*exp(a(21) + a(22)*Y + a(23)*Z + a(24)*Y*Z))^2;
delGAMMA_delLNrho_e = (1/log(10))*delGAMMA_delY;

delGAMMA_delZ = a(3) + a(4)*Y + 2*a(6)*Z + a(7)*Y^2 + 2*a(8)*Y*Z + 3*a(10)*Z^2 ....
              + (a(13) + a(14)*Y + 2*a(16)*Z + a(17)*Y^2 + 2*a(18)*Y*Z ....       % NASA Doc is wrong in this line [2a17*Y^2 -> a17*Y^2]
              + 3*a(20)*Z^2)/(1 + a(25)*exp(a(21) + a(22)*Y ....
              + a(23)*Z + a(24)*Y*Z)) - a(25)*(a(11) + a(12)*Y + a(13)*Z ....
              + a(14)*Y*Z + a(15)*Y^2 + a(16)*Z^2 + a(17)*Y^2*Z ....
              + a(18)*Y*Z^2 + a(19)*Y^3 + a(20)*Z^3)*(a(23) + a(24)*Y)* ....
                (exp(a(21) + a(22)*Y + a(23)*Z + a(24)*Y*Z))/ ....
                (1 + a(25)*exp(a(21) + a(22)*Y + a(23)*Z + a(24)*Y*Z))^2;
delGAMMA_delLNe_rho = (1/log(10))*delGAMMA_delZ;

a = sqrt(e*((gamma - 1)*(gamma + delGAMMA_delLNe_rho) + delGAMMA_delLNrho_e));

% if iflag == -1 || jflag == 1
%     return;
% end
% if jflag == 0
%     alow = a;
%     a = alow+(ahigh-alow)/(Yhigh-Ylow)*(Ym-Ylow);
%     rho = rhosave;
% end
end


%--------------------------------------------------------------------------
%   Temperature [T = T(e,rho)]
%--------------------------------------------------------------------------
function T = TFUN(rho,rho0,p,p0,T0,R,B)

Y = log10(rho/rho0);
X = log10(p/p0);
Z = X-Y;

if Z <= 0.25
    T = p/(rho*R);
else
    if -7.0 <= Y && Y <= -4.5
        layer = 1;
        if 0.25 < Z && Z <= 0.95
            column = 1;
        end
        if 0.95 < Z && Z <= 1.40
            column = 2;
        end
        if 1.40 < Z && Z <= 1.95
            column = 3;
        end
        if 1.95 < Z
            column = 4;
        end
    end
    if -4.5 < Y && Y <= -0.5
        layer = 2;
        if 0.25 < Z && Z <= 0.95
            column = 1;
        end
        if 0.95 < Z && Z <= 1.40
            column = 2;
        end
        if 1.40 < Z && Z <= 2.00
            column = 3;
        end
        if 2.00 < Z
            column = 4;
        end
    end
%     if -0.5 < Y && Y <= 3.0
    if -0.5 < Y
        layer = 3;
        if 0.25 < Z && Z <= 0.95
            column = 1;
        end
        if 0.95 < Z && Z <= 1.45
            column = 2;
        end
        if 1.45 < Z
            column = 3;
        end
    end
    
    b = B(:,column,layer);
    
%     disp('b');
%     disp([column layer]);
%     disp(b)
    
    T = T0*(10^(b(1) + b(2)*Y + b(3)*Z + b(4)*Y*Z + b(5)*Y^2 + b(6)*Z^2 ....
                + b(7)*Y^2*Z + b(8)*Y*Z^2 + b(9)*Y^3 + b(10)*Z^3 ....
                + (b(11) + b(12)*Y + b(13)*Z + b(14)*Y*Z ....
                + b(15)*Y^2 + b(16)*Z^2 + b(17)*Y^2*Z + b(18)*Y*Z^2 ....
                + b(19)*Y^3 + b(20)*Z^3)/(1 + exp(b(21) ....
                + b(22)*Y + b(23)*Z + b(24)*Y*Z))));
end
end


%--------------------------------------------------------------------------
%   Entropy [s = s(e,rho)]
%--------------------------------------------------------------------------
function s = sFUN(rho,rho0,e,R,T0,E)

Y = log10(rho/rho0);
Z = log10(e/(R*T0));

if Z <= 0.65
    s = 6779.2004 + (2.5*(Z - 0.4) - Y)*2.302585*R;
else
    if -7.0 <= Y && Y <= -4.5
        column = 1;
    end
    if -4.5 < Y && Y <= 0.5
        column = 2;
    end
    if 0.5 < Y && Y <= 3.0
        column = 3;
    end
    
    e = E(:,column);
    
%     disp('e')
%     disp(column);
%     disp(e);
    
    s = R*(e(1) + e(2)*Y + e(3)*Z + e(4)*Y*Z + e(5)*Y^2 + e(6)*Z^2 ....
           + e(7)*Y^2*Z + e(8)*Y*Z^2 + e(9)*Y^3 + e(10)*Z^3);
end
end


%--------------------------------------------------------------------------
%   Enthalpy [h = h(p,rho)]
%--------------------------------------------------------------------------
function h = hFUN(rho,rho0,p,p0,C)

Y = log10(rho/rho0);
X = log10(p/p0);
Z = X - Y;

if -7.0 <= Y && Y <= -4.5
    layer = 1;
    if Z <= 0.10
        column = 1;
    end
    if 0.10 < Z && Z <= 0.85
        column = 2;
    end
    if 0.85 < Z && Z <= 1.30
        column = 3;
    end
    if 1.30 < Z && Z <= 1.95
        column = 4;
    end
    if 1.95 < Z
        column = 5;
    end
end
if -4.5 < Y && Y <= -0.5
    layer = 2;
    if Z <= 0.10
        column = 1;
    end
    if 0.10 < Z && Z <= 0.95
        column = 2;
    end
    if 0.95 < Z && Z <= 1.50
        column = 3;
    end
    if 1.50 < Z && Z <= 2.00
        column = 4;
    end
    if 2.00 < Z
        column = 5;
    end
end
% if -0.5 < Y && Y <= 3.0
if -0.5 < Y
    layer = 3;
    if Z <= 0.10
        column = 1;
    end
    if 0.10 < Z && Z <= 1.05
        column = 2;
    end
    if 1.05 < Z && Z <= 1.60
        column = 3;
    end
    if 1.60 < Z
        column = 4;
    end
end

c = C(:,column,layer);

% disp('c');
% disp([column layer]);
% disp(c)

gamma = c(1) + c(2)*Y + c(3)*Z + c(4)*Y*Z + c(5)*Y^2 ....
        + c(6)*Z^2 + c(7)*Y^2*Z + c(8)*Y*Z^2 ....
        + c(9)*Y^3 + c(10)*Z^3 + (c(11) + c(12)*Y ....
        + c(13)*Z + c(14)*Y*Z + c(15)*Y^2 + c(16)*Z^2 ....
        + c(17)*Y^2*Z + c(18)*Y*Z^2 + c(19)*Y^3 ....
        + c(20)*Z^3)/(1 + c(25)*exp(c(21) + c(22)*Y .....
        + c(23)*Z + c(24)*Y*Z));

h = (p/rho)*(gamma/(gamma-1));
end


%--------------------------------------------------------------------------
%   Temperature [T = T(p,rho)]
%--------------------------------------------------------------------------
function T = TFUN2(rho,rho0,p,p0,T0,R,D)

Y = log10(rho/rho0);
X = log10(p/p0);
Z = X - Y;

if Z <= 0.25
    T = p/(rho*R);
else
    if -7.0 <= Y && Y <= -4.5
        layer = 1;
        if 0.25 < Z && Z <= 0.95
            column = 1;
        end
        if 0.95 < Z && Z <= 1.40
            column = 2;
        end
        if 1.40 < Z && Z <= 1.95
            column = 3;
        end
        if 1.95 < Z
            column = 4;
        end
    end
    if -4.5 < Y && Y <= -0.5
        layer = 2;
        if 0.25 < Z && Z <= 0.95
            column = 1;
        end
        if 0.95 < Z && Z <= 1.45
            column = 2;
        end
        if 1.45 < Z && Z <= 2.05
            column = 3;
        end
        if 2.05 < Z
            column = 4;
        end
    end
    if -0.5 < Y && Y <= 3.0
        layer = 3;
        if 0.25 < Z && Z <= 1.00
            column = 1;
        end
        if 1.00 < Z && Z <= 1.45
            column = 2;
        end
        if 1.45 < Z
            column = 3;
        end
    end
    
    d = D(:,column,layer);
    
%     disp('d');
%     disp([column layer]);
%     disp(d);
    
    T = T0*(10^(d(1) + d(2)*Y + d(3)*Z ....
                + d(4)*Y*Z + d(5)*Y^2 + d(6)*Z^2 + d(7)*Y^2*Z ....
                + d(8)*Y*Z^2 + d(9)*Y^3 + d(10)*Z^3 ....
                + (d(11) + d(12)*Y + d(13)*Z + d(14)*Y*Z ....
                + d(15)*Y^2 + d(16)*Z^2 + d(17)*Y^2*Z .....
                + d(18)*Y*Z^2 + d(19)*Y^3 + d(20)*Z^3)/(1 + exp(d(21) .....
                + d(22)*Y + d(23)*Z + d(24)*Y*Z))));
end
end


%--------------------------------------------------------------------------
%   Density [rho = rho(p,s)]
%--------------------------------------------------------------------------
function rho = rhoFUN(s,R,p,p0,rho0,F)

Y = log10(s/R);
X = log10(p/p0);
Z = X - Y;

if Y < 1.23
    rho = rho0*exp((log(p/p0)/1.4)-(s-s0)/(3.5*R));
else
    if 1.23 <= Y && Y < 1.42
        column = 1;
    end
    if 1.42 <= Y && Y < 1.592
        column = 2;
    end
    if 1.592 <= Y && Y < 1.70
        if Z <= (7.5269*Y-14.9366)
            column = 3;
        end
        if Z > (7.5269*Y-14.9366)
            column = 4;
        end
    end
    if 1.70 <= Y && Y < 1.180
        if Z <= (7.5269*Y-14.9366)
            column = 5;
        end
        if Z > (7.5269*Y-14.9366)
            column = 6;
        end
    end
    if 1.80 <= Y && Y < 1.90
        column = 7;
    end
    if 1.90 <= Y && Y < 2.00
        column =8;
    end
    if 2.00 <= Y && Y < 2.10
        column = 9;
    end
    if 2.10 <= Y
         column = 10;
    end
    f = F(:,column);
    
%     disp('f');
%     disp(column);
%     disp(f);
    
    rho = rho0*(10^(f(1) + f(2)*Y + f(3)*Z + f(4)*Y*Z + f(5)*Y^2 ....
                    + f(6)*Z^2 + f(7)*Y^2*Z + f(8)*Y*Z^2 ....
                    + f(9)*Y^3 + f(10)*Z^3 + (f(11) + f(12)*Y ....
                    + f(13)*Z + f(14)*Y*Z + f(15)*Y^2 + f(16)*Z^2 ....
                    + f(17)*Y^2*Z + f(18)*Y*Z^2 + f(19)*Y^3 ....
                    + f(20)*Z^3)/(1 + exp(f(21) ....
                    + f(22)*Y + f(23)*Z + f(24)*X + f(25)*Y^2))));
end
end


%--------------------------------------------------------------------------
%   Energy [e = e(p,s)]
%--------------------------------------------------------------------------
function e = eFUN(s,R,p,p0,T0,G)

Y = log10(s/R);
X = log10(p/p0);
Z = X - Y;

if Y < 1.23
    e = 2.5*e0*exp((log(p/p0)+(s-s0)/R)/3.5);
else
    if 1.23 <= Y && Y <= 1.40
        column = 1;
    end
    if 1.40 < Y && Y < 1.592
        column = 2;
    end
    if 1.592 <= Y && Y < 1.70
        if Z <= (-1.917*Y + 0.092)
            column = 3;
        end
        if Z > (-1.917*Y + 0.092)
            column = 4;
        end
    end
    if 1.70 <= Y && Y < 1.180
        if Z <= (-1.917*Y + 0.092)
            column = 5;
        end
        if Z > (-1.917*Y + 0.092)
            column = 6;
        end
    end
    if 1.80 <= Y && Y < 1.90
        column = 7;
    end
    if 1.90 <= Y && Y < 2.00
        column =8;
    end
    if 2.00 <= Y && Y < 2.10
        column = 9;
    end
    if 2.10 <= Y
         column = 10;
    end
    g = G(:,column);
    
%     disp('g');
%     disp(column);
%     disp(g);
    
    e = R*T0*(10^(g(1) + g(2)*Y + g(3)*Z + g(4)*Y*Z ....
                  + g(5)*Y^2 + g(6)*Z^2 + g(7)*Y^2*Z ....
                  + g(8)*Y*Z^2 + g(9)*Y^3 + g(10)*Z^3 ....
                  + (g(11) + g(12)*Y + g(13)*Z + g(14)*Y*Z ....
                  + g(15)*Y^2 + g(16)*Z^2 + g(17)*Y^2*Z ....
                  + g(18)*Y*Z^2 + g(19)*Y^3 + g(20)*Z^3)/(1 + exp(g(21) ....
                  + g(22)*Y + g(23)*Z + g(24)*X + g(25)*Y^2))));
end
end


%--------------------------------------------------------------------------
%   Speed of Sound [a = a(p,s)]
%--------------------------------------------------------------------------
function a = aFUN2(s,R,p,p0,a0,H)

Y = log10(s/R);
X = log10(p/p0);
Z = X - Y;

if Y < 1.23
    a = sqrt(exp(log(1.4*p0/rho0)+(log(p/p0)+(s-s0)/R)/3.5));
else
    if 1.23 <= Y && Y <= 1.40
        column = 1;
    end
    if 1.40 < Y && Y < 1.595
        column = 2;
    end
    if 1.595 <= Y && Y < 1.693
        if Z <= (-9.842*Y + 14.19)
            column = 3;
        end
        if Z > (-9.842*Y + 14.19)
            column = 4;
        end
    end
    if 1.693 <= Y && Y < 1.180
        if Z <= (-1.917*Y + 0.092)
            column = 5;
        end
        if Z > (-1.917*Y + 0.092)
            column = 6;
        end
    end
    if 1.80 <= Y && Y < 1.90
        column = 7;
    end
    if 1.90 <= Y && Y < 2.00
        column =8;
    end
    if 2.00 <= Y && Y < 2.10
        column = 9;
    end
    if 2.10 <= Y
         column = 10;
    end
    h = H(:,column);
    
%     disp('h');
%     disp(column);
%     disp(h);
    
    a = a0*(10^(h(1) + h(2)*Y + h(3)*Z + h(4)*Y*Z ....
                + h(5)*Y^2 + h(6)*Z^2 + h(7)*Y^2*Z + h(8)*Y*Z^2 ....
                + h(9)*Y^3 + h(10)*Z^3 + (h(11) + h(12)*Y ....
                + h(13)*Z + h(14)*Y*Z + h(15)*Y^2 ....
                + h(16)*Z^2 + h(17)*Y^2*Z + h(18)*Y*Z^2 ....
                + h(19)*Y^3 + h(20)*Z^3)/(1 + exp(h(21) ....
                + h(22)*Y + h(23)*Z + h(24)*X + h(25)*Y^2))));
end
end


%--------------------------------------------------------------------------
%   Iteratively Solves [p = p(h,s), e = e(h,s), rho = rho(h,s)]
%--------------------------------------------------------------------------
function [p, e, rho] = perho(h,s,s0,p0,T0,rho0,h0,R,G,F)

pp(50) = 0;
ff(50) = 0;

pp(1) = exp(-(s-s0)/R)*(h/h0)^5;
pp(2) = 1e-6*pp(1);

for i = 1:20        %   Bracket solution
    for j = 1:2
        p = pp(j)*p0;
        e = eFUN(s,R,p,p0,T0,G);
        rho = rhoFUN(s,R,p,p0,rho0,F);
        hn = e + p/rho;
        ff(j) = (hn - h)/h;
    end
    if ff(1)*ff(2) > 0.0
        pp(1) = 2*pp(1);
        pp(2) = (1/2)*pp(2);
    else
        break;
    end
end

for j = 1:50        %   Regular Falsi
    ppn = (pp(2)*ff(1)-pp(1)*ff(2))/(ff(1)-ff(2));
    p = ppn*p0;
    e = eFUN(s,R,p,p0,T0,G);
    rho = rhoFUN(s,R,p,p0,rho0,F);
    hn = e + p/rho;
    ffn = (hn - h)/h;
    affn = abs(ffn);
    if affn < 0.1
        break;
    end
    if ffn*ff(1) > 0.0
        pp(1) = ppn;
        ff(1) = ffn;
    else
        pp(2) = ppn;
        ff(2) = ffn;
    end
end

for j = 3:50        %   Secant Method
    pp(j) = pp(j-1) - ff(j-1)*(pp(j-1)-pp(j-2))/(ff(j-1)-ff(j-2));
    p = pp(j)*p0;
    e = eFUN(s,R,p,p0,T0,G);
    rho = rhoFUN(s,R,p,p0,rho0,F);
    hn = e + p/rho;
    ff(j) = (hn - h)/h;
    if abs(ff(j)) < 5e-4
        break;
    end
end
end


%--------------------------------------------------------------------------
%   Iteratively Solves [rho = rho(p,T)]
%--------------------------------------------------------------------------
function rho = rhoFUN2(p,T,p0,T0,rho0,R,D)

%gg(50) = 0;

rr(1) = p*(T0/T)/p0;
rr(2) = 0.9*rr(1);
Temp = T;

for ic = 1:20       %   Bracket Solution
    gg = zeros(2);
    for j = 1:2
        rho = rr(j)*rho0;
        T = TFUN2(rho,rho0,p,p0,T0,R,D);
        gg(j) = (T - Temp)/Temp;
    end
    if gg(1)*gg(2) > 0.0
        rr(1) = 2*rr(1);
        rr(2) = (1/2)*rr(2);
    else
        break;
    end
end

for j = 1:50        %   Regular Falsi
    rrn = (rr(2)*gg(1)-rr(1)*gg(2))/(gg(1)-gg(2));
    rho = rrn*rho0;
    T = TFUN2(rho,rho0,p,p0,T0,R,D);
    ggn = (T - Temp)/Temp;
    aggn = abs(ggn);
    if aggn < 0.1
        break;
    end
    if ggn*gg(1) > 0.0
        rr(1) = rrn;
        gg(1) = ggn;
    else
        rr(2) = rrn;
        gg(2) = ggn;
    end
end

rr = zeros(50);
for j = 3:50
    rr(j) = rr(j-1) - gg(j-1)*(rr(j-1)-rr(j-2))/(gg(j-1)-gg(j-2));
    rho = rr(j)*rho0;
    T = TFUN2(rho,rho0,p,p0,T0,R,D);
    gg(j) = (T - Temp)/Temp;
    if abs(gg(j)) < 1e-10
        break;
    end
end
%T = Temp;
end







