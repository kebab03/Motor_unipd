% This script allows the rotor flux orientation to be achieved.

rapp1 = Isq/Irq;
rapp2 = Irq/Isq;
if abs(rapp1) > 1.15 || abs(rapp2) < 1/1.15
    Isq = -Irq;
end

% Initialize the rotor q-axis flux:
if Irq == 0
    opendocument([filename, '.fem']);
    solving_core;
    Ls = fluxD_S/Isd;
    Lm = fluxD_R/Isd;
else
    fluxQ_R = 1;
    while abs(fluxQ_R) > 1e-4
        
        %%%%%%%%%%%%%%%% FIRST SIMULATION %%%%%%%%%%%%%%%%
        
        opendocument([filename, '.fem']);
        solving_core;
        x1 = fluxQ_R;
        y1 = Irq;
        
        %%%%%%%%%%%%%%%% SECOND SIMULATION %%%%%%%%%%%%%%%%
        
        if Isd == 0
            Irq = y1-fluxQ_R;
        else
            Irq = y1-fluxQ_R/fluxD_R*Isd;
        end
        opendocument([filename, '.fem']);
        solving_core;        
        x2 = fluxQ_R;
        y2 = Irq;
                
        %%%%%%%%%%%%%%%% THIRD SIMULATION %%%%%%%%%%%%%%%%
        
        Irq  = (x2*y1-x1*y2)/(x2-x1);
        opendocument([filename, '.fem']);
        solving_core;
        mi_close;
        
        x3 = fluxQ_R;
        y3 = Irq;
        
        fluxQ_R
    end
    Irq_FO = y3;
    Irq = Irq_FO;
      
    n  = -Irq_FO/Isq;
    Ls = fluxD_S/Isd;
    Lt = fluxQ_S/Isq;
    Lphi = Ls-Lt;
    Lm = fluxD_R/Isd;
    Lr = Lm^2/Lphi;
end