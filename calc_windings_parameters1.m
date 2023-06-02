function [Lew, Rs, sigma_cage, sigma_cage_20, Kring, Sslot_S, Sslot_R, Lring, Rring] = calc_windings_parameters1(Tstator, Trotor, filename, stator, rotor, struct)

%% Computation of the stator and rotor slot cross-section area:

% here simulate the machine without any excitation, compute the stator
% and rotor slot area from the solution. Use the string filename to open
% the model.

openfemm(1)
opendocument([filename,'.fem'])
mi_analyze()
mi_loadsolution()

mo_groupselectblock(stator.group + 1); % first stator slot selection
Sslot_S = mo_blockintegral(5); % [m^2], stator slot area
mo_clearblock();

mo_groupselectblock(rotor.group + 1); % first rotor slot selection
Sslot_R = mo_blockintegral(5); % [m^2], rotor slot area
mo_clearblock();

closefemm()

%% End-winding Inductance:
mu_0 = 4*pi*1e-7;
p    = stator.p;
lstk = stator.Lstk;
Di   = stator.Di;
hs   = stator.hs;
Qs   = stator.Q;
ms   = stator.winding.m;
qs   = Qs/(ms*2*p);
ncs  = stator.winding.ncs;

Kfill = stator.winding.kfill;
lew   = 2.5*Di/p;

% The end-winding leakage inductance:
Lew   = mu_0*2*p*ncs^2*qs^2*lew*0.5; % [H]

%% Stator resistance
rho_Slots = stator.slot.material.rho*( 1 + stator.slot.material.alphaT*(Tstator-20) );
Rs        = rho_Slots*(ncs^2)*(2*qs*p)*(lstk + lew)/(Sslot_S*Kfill); % [Ohm], stator resistance

%% Rotor Aluminum Conductivity
sigma_cage_20 = rotor.slot.material.sigma/1e6;
rho_Slotr     = rotor.slot.material.rho*(1 + rotor.slot.material.alphaT*(Trotor-20));
sigma_cage    = 1/rho_Slotr/1e6;

%% Rotor Ring Coefficient
Qr     = rotor.Q;
D_ring = rotor.ring.D;
S_ring = rotor.ring.S;
Kring  = (2/pi)*(Qr/(2*p)^2)*(D_ring/lstk)*(Sslot_R/S_ring);

%% Rotor Ring Coefficient
W_ring       = rotor.ring.W;
H_ring       = rotor.ring.H;
lambda_a     = 0.365*log10(3*pi*(D_ring)/(4*(W_ring + H_ring)));
Lring        = pi*mu_0/Qr*D_ring*lambda_a;
Rring        = 1/sigma_cage*(pi*D_ring/Qr)/(S_ring*1e6);
end