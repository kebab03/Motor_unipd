% *************************************************************************
%         INDUCTION MOTOR ON-LOAD STATIC SIMULATION
% *************************************************************************
%
% This is the MAIN program
%
% This script allows the IM on-load permormance prediction, performing
% magneto-static analyses. The rotor is computed using the iterative
% procedure, based on the machine model in the Rotor Field Oriented (RFO)
% reference frame. This procedure is done in the file RFOA.m.
%
% Here you need:
%  1. MOTOR_DATA_single_cage.m loaded:
%  2. launch the function for the stator winding and rotor cage parameters
%     computation;
%  3. MOST IMPORTANT! Define stator d- and q-axis currents
%
%
% 2020/04/08
% *************************************************************************

clc; 
clear; 
close all;

% Reading motor data

  MOTOR_DATA_single_cage;

% Computation of the equivalent rotor winding

  rotor_winding;
  
% stacking coefficienct

  K_fe = stator.K_fe;

% Operating winding temperature

  Tstator = 120;        % stator temperature
  Trotor  = 150;        % rotor  temperature
  
%% rated data

  Vnom = 400; % [V] rms rated voltage !!! winding voltage !!!
  fs = 50; % [Hz] stator frequency
  Pnom = 7500; % [W] rated power
  RPM0 = 60*fs/stator.p; % [1/sec] no load speed

  Sest = Pnom/0.9/0.85; % [VA] estimated apparent power with eta = 0.9, cosphi = 0.85
  Iest = Sest/3/Vnom; % [A] rms estimated stator effective current
  Ws = 2*pi*fs; % [rad/sec] stator angular frequency

%% winding parameters

  [Lew, Rs, sigma_cage, sigma_cage_20, Kring, Sslot_S, Sslot_R, Lring, Rring] =...
    calc_windings_parameters1(Tstator, Trotor, filename, stator, rotor);
  stator.winding.R = Rs;
  stator.winding.Lew = Lew;
  rotor.winding.Kring = Kring;
  
%% Open FEMM

  openfemm(1)
  
%% No-Load test

% This test is meant to derive the voltage characteristic of the machine
% and determine the stator magnetizing current at the rated voltage

if ~exist('NLres.mat','file')
    
    % initialize variables
    Isd_vec = linspace(Iest*0.025*sqrt(2),Iest*sqrt(2),15);
    fluxD_S_vec = ones(size(Isd_vec));
    fluxD_R_vec = ones(size(Isd_vec));
    Ls_vec = ones(size(Isd_vec));
    Lm_vec = ones(size(Isd_vec));
    Vsq_vec = ones(size(Isd_vec));
    
    for isdIdx = 1:length(Isd_vec)
        
        % stator currents
        Isd = Isd_vec(isdIdx);
        Isq = 0;
        
        % rotor currents
        Ird = 0;
        Irq = 0;
        
        % perform no-load simulation
        RFOA;
        Vsq = Ws*fluxD_S;
        
        % fill results vectors
        fluxD_S_vec(isdIdx) = fluxD_S; % [V s] stator d-axis flux linkage
        fluxD_R_vec(isdIdx) = fluxD_R; % [V s] rotor d-axis flux linkage
        Ls_vec(isdIdx) = Ls; % [H] stator synchronous inductance
        Lm_vec(isdIdx) = Lm; % [H] mutual synchronous inductance
        Vsq_vec(isdIdx) = Vsq; % [V] stator phase peak voltage
    end
    save('NLres.mat')
else
   load('NLres.mat')   
end

% stator voltage
figure(1), grid on; hold on; box on;
plot(Isd_vec/sqrt(2),Vsq_vec/sqrt(2),'.-','linewidth',1.5,'markersize',17.5)
xlabel('magnetizing current I_{sd} [A] rms')
ylabel('phase voltage V_{sq} [V] rms')

% magnetizing inductances
figure(2), grid on; hold on; box on;
plot(Isd_vec/sqrt(2),Ls_vec*1000,'.-','linewidth',1.5,'markersize',17.5)
plot(Isd_vec/sqrt(2),Lm_vec*1000,'.-','linewidth',1.5,'markersize',17.5)
legend('stator synchronous','mutual synchronous')
xlabel('magnetizing current I_{sd} [A] rms')
ylabel('magnetizing inductances [mH]')

% stator & rotor d-axis flux linkage
figure(3), grid on; hold on; box on;
plot(Isd_vec/sqrt(2),fluxD_S_vec,'.-','linewidth',1.5,'markersize',17.5)
plot(Isd_vec/sqrt(2),fluxD_R_vec,'.-','linewidth',1.5,'markersize',17.5)
legend('stator flux','rotor flux')
xlabel('magnetizing current I_{sd} [A] rms')
ylabel('d-axis flux \lambda_{sd}, \lambda_{sq} [V s]')

% compute the rated magnetizing current
% consider a voltage drop of 5% in the stator impedance
Imu_nom = interp1(Vsq_vec/sqrt(2), Isd_vec/sqrt(2), Vnom*0.95); % [A] rms

figure(1), hold on;
plot([0 1]*max(Isd_vec/sqrt(2)), [1 1]*Vnom*0.95, 'k', 'linewidth', 1)
text(0.1*max(Isd_vec/sqrt(2)), 1.05*Vnom*0.95, ['V_{nom} = ',num2str(Vnom*0.95)])

plot([1 1]*Imu_nom, [0 1]*max(Vsq_vec/sqrt(2)), 'k', 'linewidth', 1)
text(1.05*Imu_nom, 0.05*max(Vsq_vec/sqrt(2)), ['I_{\mu} = ',num2str(Imu_nom)])

%% RFOA on-load simulation

if ~exist('OLres.mat','file')

% stator rated magnetizing current
  Isd = Imu_nom*sqrt(2); % [A] peak
  
% initialize the torque current vector
  Isq_vec = linspace(0.25*Iest*sqrt(2),3*Iest*sqrt(2),10);
  
% initialize results vectors
  Is_vec = ones(size(Isq_vec));     % stator current
  Irq_FO_vec = ones(size(Isq_vec)); % rotor current
  Lt_vec = ones(size(Isq_vec));     % transient inductance
  Lphi_vec = ones(size(Isq_vec));   % magnetizing inductance
  Tdq_vec = ones(size(Isq_vec));    % torque
  Wsl_vec = ones(size(Isq_vec));    % slip angular frequency
  Pmecc_vec = ones(size(Isq_vec));  % mechalical power
  slip_vec = ones(size(Isq_vec));   % slip
  PJS_vec = ones(size(Isq_vec));    % stator Joule losses
  PJR_vec = ones(size(Isq_vec));    % rotor Joule losses
  
  for isqIdx = 1:length(Isq_vec)
      
      % stator  current
        Isq = Isq_vec(isqIdx); % [A] peak
        Is = hypot(Isd,Isq); % [A] peak 
        Is_vec(isqIdx) = Is; % [A] peak
       
      % initialize rotor currents
        Ird = 0;
        Irq = -Isq;
        
      % launch Rotor Field Oriented Analysis
        RFOA;
        PjS = (3/2)*Rs*Is^2; % [W] stator Joule losses
        slip = Wsl/Ws; % [-] slip
        Pmecc = Tdq*Ws/stator.p*(1-slip); % [W] mech power
        
      % store results
        Irq_FO_vec(isqIdx) = Irq_FO; % [A] peak
        Lt_vec(isqIdx) = Lt; % [H]
        Lphi_vec(isqIdx) = Lphi; % [H]
        Tdq_vec(isqIdx) = Tdq; % [N m]
        Wsl_vec(isqIdx) = Wsl; % [rad/sec]
        Pmecc_vec(isqIdx) = Pmecc; % [W]
        slip_vec(isqIdx) = slip; % [-]
        PJS_vec(isqIdx) = PjS; % [W]
        PJR_vec(isqIdx) = PjR; % [W]
  end
  closefemm()
  save('OLres.mat');
else
    load('OLres.mat');
end

% rated values
slip_nom = interp1(Pmecc_vec, slip_vec, Pnom); % [-] rated slip
T_nom = interp1(Pmecc_vec, Tdq_vec, Pnom); % [N m] rated torque
Is_nom = interp1(Pmecc_vec, Is_vec/sqrt(2), Pnom); % [A] rms rated current

% mechanical power
figure(4), grid on; hold on, box on;
plot(slip_vec, Pmecc_vec, '.-', 'linewidth', 1.5, 'markersize', 17.5)
plot([0 1]*max(slip_vec), [1 1]*Pnom, 'k', 'linewidth', 0.75)
plot([1 1]*slip_nom, [0 1]*max(Pmecc_vec), 'k', 'linewidth', 0.75)
text(0.025*max(slip_vec),1.05*Pnom,['P_{nom} = ',num2str(Pnom),' [W]'])
text(1.05*slip_nom,0.025*max(Pmecc_vec),['slip_{nom} = ',num2str(slip_nom)])
xlabel('slip')
ylabel('mechanical power [W]')

% torque
figure(5), grid on; hold on, box on;
plot(slip_vec, Tdq_vec, '.-', 'linewidth', 1.5, 'markersize', 17.5)
plot([0 1]*max(slip_vec), [1 1]*T_nom, 'k', 'linewidth', 0.75)
plot([1 1]*slip_nom, [0 1]*max(Tdq_vec), 'k', 'linewidth', 0.75)
text(0.025*max(slip_vec),1.05*T_nom,['T_{nom} = ',num2str(T_nom),' [N m]'])
text(1.05*slip_nom,0.025*max(Tdq_vec),['slip_{nom} = ',num2str(slip_nom)])
xlabel('slip')
ylabel('torque [N m]')

% stator current
figure(6), grid on; hold on, box on;
plot(slip_vec, Is_vec/sqrt(2), '.-', 'linewidth', 1.5, 'markersize', 17.5)
plot([0 1]*max(slip_vec), [1 1]*Is_nom, 'k', 'linewidth', 0.75)
plot([1 1]*slip_nom, [0 1]*max(Is_vec/sqrt(2)), 'k', 'linewidth', 0.75)
text(0.025*max(slip_vec),1.05*Is_nom,['I_{snom} = ',num2str(Is_nom),' [A] rms'])
text(1.05*slip_nom,0.025*max(Is_vec),['slip_{nom} = ',num2str(slip_nom)])
xlabel('slip')
ylabel('stator current [A] rms')

