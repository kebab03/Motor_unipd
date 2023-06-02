%% POST-PROCESSING

mi_loadsolution();

% post-processing operation:
%   1. load the solution
%   2. computation of the slot area
%   3. computation of the flux linkages
%   4. computation of the torque
%   5. computation of rotor Joule losses

%% Computation of the stator and rotor slot cross-section area

mo_groupselectblock(stator.group + 1); % first stator slot selection
Sslot_S = mo_blockintegral(5); % [m^2], stator slot area
mo_clearblock();

mo_groupselectblock(rotor.group + 1); % first rotor slot selection
Sslot_R = mo_blockintegral(5); % [m^2], rotor slot area
mo_clearblock();

%% Computation of the flux linkages

% Stator flux linkages
Az_S(QsimS,1) = 0; % Initialize the vector potential vector
for qs = 1:QsimS
    % !!! compute here the average integral of the vector potential in each
    % stator slot, separately and store them in the vector. Bare in mind to
    % multiply by the length and divide by the stator slot area.
    mo_groupselectblock(stator.group + qs);
    Az_S(qs,1) = mo_blockintegral(1)/Sslot_S;
    mo_clearblock();
end

% The moltiplication for Lstk is already taken into account (in FEMM!)
FluxABC_S = sim_period*ncs_S*K_S(1:QsimS,:)'*real(Az_S);
[fluxD_S, fluxQ_S] = calc_abc2dq(FluxABC_S(1), FluxABC_S(2), FluxABC_S(3), 0);


% Rotor flux linkages
Az_R(QsimR,1) = 0;
for qr = 1:QsimR
    % !!! compute here the average integral of the vector potential in each
    % rotor slot, separately and store them in the vector. Bare in mind to
    % multiply by the length and divide by the rotor slot area.
    mo_groupselectblock(rotor.group + qr);
    Az_R(qr,1) = mo_blockintegral(1)/Sslot_R;
    mo_clearblock();
end

% The moltiplication for Lstk is already taken into account (in FEMM!)
FluxABC_R = sim_period*ncs_R*K_Rabc(1:QsimR,:)'*real(Az_R);
[fluxD_R, fluxQ_R] = calc_abc2dq(FluxABC_R(1), FluxABC_R(2), FluxABC_R(3), 0);

%% Computation of the torque and radial force (selecting the rotor)

Tdq =  (3/2)*stator.p*(fluxD_S*Isq - fluxQ_S*Isd); % [N m]

%% Computation of rotor Joule Losses

for qr = 1:rotor.Qsim
    mo_groupselectblock(rotor.group + qr);
end
RotLosses = sim_period*mo_blockintegral(4);

% adjust considering the temperature and the ring factor
PjR = RotLosses*sigma_cage_20/sigma_cage*(1+rotor.winding.Kring); % [W]

%% Slip Angular Frequency

Wsl = stator.p*PjR/Tdq; % [rad/sec]
