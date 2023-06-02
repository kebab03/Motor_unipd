% -----------------------------------------------------------
%     SOLVING_CORE 
% -----------------------------------------------------------
% This script acts as the common ground for all the simulations to come.
%
% The requirements are:
%   A. Stator Isd and Isq currents
%   B. Rotor  Ird and Irq currents
%   C. MACHINE_DATA loaded
%
% It basically performs :
% two pre-processing operations:
%   1. currents computation
%   2. currents setting
% three simulation operations:
%   1. problem definition
%   2. saving of a "temp" file
%   3. analysis
% and other post-processing operations:
%   1. load the solution
%   2. computation of the slot area
%   3. computation of the flux linkages
%
%% PRE-PROCESSING ---------------------------------------------------------

% load the stator winding data

K_S   = stator.winding.K_S;
ncs_S = stator.winding.ncs;
QsimS = stator.Qsim;

% load rotor sinusoidal equivalent winding data

K_Rabc = rotor.winding.K_Rabc;
ncs_R  = rotor.winding.ncs;
QsimR  = rotor.Qsim;
    
% compute stator phase curents currents

[isa,isb,isc] = calc_dq2abc(Isd, Isq, 0);

% compute rotor phase currents

[ira,irb,irc] = calc_dq2abc(Ird, Irq, 0);
      
% set the stator slot current

Islot_S_vec = ones(1,QsimS);
for qs = 1:QsimS
    % the slot current is computed considering the phase currents 
    % and the stator winding matrix:
      Islot_S = ncs_S*(isa * K_S(qs,1) + isb * K_S(qs,2) + isc * K_S(qs,3));
    % here the current is set into the proper slot:
      mi_modifycircprop(['Islot', num2str(qs)], 1, Islot_S);
      Islot_S_vec(qs) = Islot_S;
end

% set the rotor slot current

Islot_R_vec = ones(1,QsimR);
for qr = 1:QsimR
    % the slot current is computed considering the phase currents 
    % and the stator winding matrix:
      Islot_R = ncs_R*(ira * K_Rabc(qr,1) + irb * K_Rabc(qr,2) + irc * K_Rabc(qr,3));
    % here the current is set into the proper slot:
      mi_modifycircprop(['IslotR', num2str(qr)], 1, Islot_R);
      Islot_R_vec(qr) = Islot_R;
end
    
%% SAVE AND ANALYSE -------------------------------------------------------

% make a folder where to save the "excited model", to leave clean the
% main one. The temp.fem model with currents is saved into the folder temp.

if ~exist('tempFolder', 'var') || ~exist('tempFN', 'var')
  tempFolder = 'temp';
  tempFN = 'temp';
  
  if ~exist('tempFolderMain','var')
    tempFolderMain = '.';
  end
  
  if ~exist([tempFolderMain, '/', tempFolder],'dir')
    mkdir([tempFolderMain, '/', tempFolder]);
  end
end

if ~exist('MeshMinAngle','var')
    % set the minimum angle for the mesh triangles:
      MeshMinAngle = 10;
end

% Use Previous solution by default for non linear iteration

if ~exist('UsePreviousSolution','var')
  UsePreviousSolution = 1;
end

if exist('UsePreviousSolution','var') && UsePreviousSolution
  if exist([tempFolder,'/',tempFN,'.ans'],'file') == 2
    mi_setprevious([tempFN,'.ans'], 0);
  end
end

% here the problem definition: frequnecy and some other variables. Then
% the simulation is run

mi_probdef(0,'meters','planar',1e-8,stator.Lstk*stator.K_fe,MeshMinAngle,0); % problem definitions
mi_saveas([tempFolderMain, '\\', tempFolder,'\\',tempFN,'.fem']); % save the file as temp
mi_analyze(1); % analyze the model

%% POST-PROCESSING ---------------------------------------------------------

post_processing;