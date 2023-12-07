function dydt = ode_system(t, y)
    % Parameters
    gZ = 0.1;
    gA = 0.01;
    gB = 0.01;
    gA1 = 0.001;
    dZ = 0.00001;
    dA = 0.00001;
    dB = 0.00001;
    dA1 = 0.00001;
    pZA = 0.2;
    pZB = 0.05;
    pAA1 = 0.2;
    uA1 = 0.1;
    uB = 0.1;
    kA1 = 0.1;
    kB = 0.1;
    qAA = 0.3;
    qB = 0.3;
    r = 0.5;
    
    % Variables
    A1 = y(1);
    B = y(2);
    tA1 = y(3);
    tB = y(4);
    
    % ODE system
    dydt = zeros(4, 1);
    dydt(1) = (gA1 - dA1) + uA1*(tA1)/(kA1+tA1)*B - uB*(tB)/(kB+tB)*A1;;
    dydt(2) = (gB - dB) - uA1*(tA1)/(kA1+tA1)*B + uB*(tB)/(kB+tB)*A1;
    dydt(3) = -qAA * tA1;
    dydt(4) = -qB * tB;
end