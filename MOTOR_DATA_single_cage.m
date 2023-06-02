% *************************************************************
%     MOTOR_DATA
%*************************************************************
% This file contains the data of the machine.
%
% It is divided into sections:
%   * I/ paths
%
%   * stator
%     * winding
%
%   * rotor
%     * winding
%
%   * material
%
%   * input (e.g. optimization parameters)
%   * initialization
%   * computations dependent variables
%   * /O paths
%
% 2020/04/08
% *************************************************************************

%% I/ paths ---------------------------------------------------------------
if ~exist('Out2File_CC', 'var')
  filename = 'IM_single_cage';
else
  filename = 'IM_single_cage_CC';
end

%% General data -----------------------------------------------------------
mm = 1e-3;                     % millimeters
lc = mm;
% This variable defines the machine type
MachineType = 'IMsc';
p           = 1;                % number of pole pairs
nphase      = 3;                % number of phases
gap         = 0.3*mm;           % air-gap length
Lstk        = 200*mm;           % stack length
sim_poles   = 1;                % poles of the motor to analyze
sim_period  = 2*p/sim_poles;    % periodicity of the analysis

%% STATOR (m,rad) ---------------------------------------------------------
stator.name       = 's';         % nickname for the stator
stator.lc         = lc;
stator.p          = p;
stator.gap        = gap;
stator.airgap     = 2/5*gap;
stator.Lstk       = Lstk;
stator.K_fe       = 0.97;
stator.poles      = 2*p;          % number of stator poles
stator.sim_poles  = sim_poles;
stator.sim_period = sim_period;

% Geometry
stator.De         = 173*mm;      % [m] stator outer diameter
stator.Di         = 90*mm;       % [m] stator inner diameter

% Slot
stator.Q          = 24;          % number of stator slots
stator.hs         = 17*mm;       % [m], slot hight
stator.wt         = 5*mm;        % [m], tooth width
stator.wso        = 2.7*mm;      % [m], slot opening width
stator.hso        = 0.43*mm;     % [m], slot opening height
stator.hwed       = 0.77*mm;     % [m], wedge height
stator.slot.r     = 3.91*mm;     % wedge fillet radius
stator.slot.R     = 2.5*mm;      % slot-end fillet radius
stator.slot.label = 'circuit';   % 'material'/'circuit'

% Materials
stator.material             = 'isovac_330_35A';
stator.slot.material.name   = 'Cu';              % 'Cu'/'Al'/...
stator.slot.material.alphaT = 0.0039;            
stator.slot.material.sigma  = 56000000;          % [S/m]
stator.slot.material.rho    = 1.7857e-08;        % [Ohm.m]

stator.group = 1000; % stator group number

%% Stator Winding 
stator.winding.turns = 33;                       % number of turns per phase
stator.winding.nc    = 66;                       % number of series conductors per slot
stator.winding.npp   = 2;                        % number of parallel paths
stator.winding.yq    = 8;                        % coil pitch
stator.winding.m     = nphase;                   % number of phases
% stator.winding.Nlayer= 2;                      % single/double layer
stator.winding.kfill = 0.38;                     % slot fill factor
stator.winding.K_S   = load('K_slot_matrix.txt');

%% ROTOR (INNER)-----------------------------------------------------------
% Geometry
rotor.name       = 'r';                           % nickname for the rotor
rotor.lc         = lc;
rotor.p          = p;
rotor.gap        = gap;
rotor.airgap     = 2/5*gap;
rotor.Lstk       = Lstk;
rotor.K_fe       = 0.985;
rotor.poles      = 2*p;
rotor.sim_poles  = sim_poles;
rotor.sim_period = sim_period;
rotor.De         = stator.Di - 2*gap;             % [m] rotor outer diameter
rotor.Di         = 40*mm;                         % [m], rotor inner diameter
rotor.Dshaft     = rotor.Di;                      % [m], rotor shaft diameter
% Rotor slot parameters ---------------------------------------------------
% Slot
rotor.Q          = 18;                             % number of rotor slots
rotor.hs         = 20.2*mm - gap;                  % [m], slot hight
rotor.wt         = 6.5*mm;                         % [m], tooth width
rotor.wso        = 1*mm;                           % [m], slot opening width
rotor.hso        = 0.74*mm - gap;                  % [m], slot opening height
rotor.slot.r     = 3.78*mm;                        % wedge fillet radius
rotor.slot.R     = 1.28*mm;                        % slot-end fillet radius
rotor.slot.label = 'circuit';                      % 'circuit'/''
% rotor.slot.shape = 'trapz'; % 'trapz'/'?!circular'/'?!arc'
% Short Circuit Rings
rotor.ring.Di    = 52*mm;                           % [m]
rotor.ring.De    = 87.2*mm;                         % [m]
rotor.ring.D     = (rotor.ring.Di+rotor.ring.De)/2; % [m]
rotor.ring.W     = 6.5*mm;                          % [m]
rotor.ring.H     = (rotor.ring.De-rotor.ring.Di)/2; % [m]
rotor.ring.S     = rotor.ring.H*rotor.ring.W;       % [m^2]
% Materials
rotor.material             = 'isovac_330_35A';      % lamination name
rotor.slot.material.name   = 'Al';                  % 'Cu'/'Al'/...
rotor.slot.material.alphaT = 0.0039;                % [K^-1]
rotor.slot.material.sigma  = 35000000;              % [S/m]
rotor.slot.material.rho    = 2.8571e-08;            % [Ohm.m]
rotor.shaft.material       = 'Steel_1010';

rotor.group = 10;

%% Calculations of dependent variables ------------------------------------
stator.slot.angle = 360/stator.Q/2;                             % angle of the first slot in the drawing
stator.hbi        = abs(stator.De/2 - stator.Di/2 - stator.hs); % height of the stator yoke
rotor.slot.angle  = 360/rotor.Q/2;                              % angle of the first slot in the drawing
rotor.hbi         = abs(rotor.De/2 - rotor.Di/2 - rotor.hs);    % height of the statoyoke

stator.winding.ncs = stator.winding.nc/stator.winding.npp;      % number of series conductor per slot
stator.winding.t = gcd(stator.Q, stator.p);                     % stator winding periodicity

sim_period         = 2*p/sim_poles;                             % periodicity of the analysis
stator.sim_poles   = sim_poles;
stator.sim_period  = sim_period;
Qsim               = stator.Q/stator.sim_period;
QsimR              = rotor.Q/rotor.sim_period;
stator.Qsim        = Qsim;
rotor.Qsim         = QsimR;
rotor.sim_poles    = sim_poles;          %%%%% check!
rotor.sim_period   = sim_period;         %%%%% check!