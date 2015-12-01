function [state] = state_init(varargin)
%===============================================================================%
%   Author: Ryan Krattiger                                                      %
%===============================================================================%
%   Inputs:                                                                     %
%       varargin -> Options List                                                %
%   Options:                                                                    %
%       'Mach'          ->                                                      %
%       'Temperature'   ->                                                      %
%       'TempTotal'     ->                                                      %
%       'Pressure'      ->                                                      %
%       'PresTotal'     ->                                                      %
%       'Density'       ->                                                      %
%       'DensTotal'     ->                                                      %
%       'Velocity'      ->                                                      %
%       'VelocityVec'   ->                                                      %
%       'SpecRatio'     ->                                                      %
%       'CP'            ->                                                      %
%       'CV'            ->                                                      %
%       'Viscosity'     ->                                                      %
%       'Energy'        ->                                                      %
%       'EnergyTotal'   ->                                                      %
%       'Enthalpy'      ->                                                      %
%       'EnthalpyTotal' ->                                                      %
%       'Entropy'       ->                                                      %
%       'SoS'           ->  Speed of sound                                      %
%       'Units'         ->  Flag for units: 1 -> Imperial, 0 -> SI              %
%   Default values for all varibles are in conditions of STP at sea level with  % 
%       zero velocity                                                           %
%-------------------------------------------------------------------------------%
%   Ouptus:                                                                     %
%       state:      structure of flow state variables                           %
%===============================================================================%

% Example: var = state_init('Mach', 4, 'Temperature', 600) 
%   sets var.M = 4 & var.T = 600, all other values remain default
if(nargin == 0)
    flag = 0;
end

state.T = 273;     % temperature [K]
state.Tt = 273;
state.p = 101325;     % pressure [Pa]
state.pt = 101325;
state.r = 1.225;     % density [kg/m^3]
state.rt = 1.225;
state.V = 0;     % velocity magnitude [m/s]
state.V_(1:3) = [0 0 0]; % velocity vector [m/s]
state.M = 0;     % Mach number
state.gma = 1.4;   % specific heat ratio
state.cp = 1004.5;    % specific heat -> constant pressure [J/K]
state.cv = state.cp/1.4;    % specific heat -> constant volume [J/K]
state.mu = 1.716*10^-5; % Viscosity coefficient
state.e = state.cv*state.T;     % energy per unit mass [J/kg]
state.et = state.e;
state.h = state.cp*state.T;     % enthalpy per unit mass [J/kg]
state.ht = state.h;
state.s = 0;     % entropy per unit mass [J/(K-kg)]
state.a = sqrt(state.gma*287*state.T);     % speed of sound [m/s]
state.u_flag = 0; %set a flag for the units

assert(mod(nargin,2)==0, 'state_init: Odd number of inputs not allowed')
for i = 1:2:nargin
    switch cell2mat(varargin(i))
        case 'Temperature'
            state.T = cell2mat(varargin(i+1));     % temperature [K]
        case 'TempTotal'
            state.Tt = cell2mat(varargin(i+1));
        case 'Pressure'
            state.p = cell2mat(varargin(i+1));     % pressure [Pa]
        case 'PresTotal'
            state.pt = cell2mat(varargin(i+1));
        case 'Density'
            state.r = cell2mat(varargin(i+1));     % density [kg/m^3]
        case 'DensTotal'
            state.rt = cell2mat(varargin(i+1));
        case 'Velocity'
            state.V = cell2mat(varargin(i+1));     % velocity magnitude [m/s]
        case 'VelocityVec'
            state.V_(1:3) = cell2mat(varargin(i+1)); % velocity vector [m/s]
        case 'Mach'
            state.M = cell2mat(varargin(i+1));     % Mach number
        case 'SpecRatio'
            state.gma = cell2mat(varargin(i+1));   % specific heat ratio
        case 'CP'
            state.cp = cell2mat(varargin(i+1));    % specific heat -> constant pressure [J/K]
        case 'CV'
            state.cv = cell2mat(varargin(i+1));    % specific heat -> constant volume [J/K]
        case 'Viscosity'
            state.mu = cell2mat(varargin(i+1));    % specific heat -> constant volume [J/K]
        case 'Energy'
            state.e = cell2mat(varargin(i+1));     % energy per unit mass [J/kg]
        case 'EnergyTotal'
            state.et = cell2mat(varargin(i+1));
        case 'Enthalpy'
            state.h = cell2mat(varargin(i+1));     % enthalpy per unit mass [J/kg]
        case 'EnthalpyTotal'
            state.ht = cell2mat(varargin(i+1));
        case 'Entropy'
            state.s = cell2mat(varargin(i+1));     % entropy per unit mass [J/(K-kg)]
        case 'SoS'
            state.a = cell2mat(varargin(i+1));     % speed of sound [m/s]
        case 'Units'
            state.u_flag = cell2mat(varargin(i+1)); %set a flag for the units
        otherwise
            error('state_init: Invalid Option')
    end
end
