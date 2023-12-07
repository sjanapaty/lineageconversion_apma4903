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
    
    % Variables
    Z = y(1);
    A = y(2);
    A1 = y(3);
    B = y(4);
    tA1 = y(5);
    tB = y(6);
    
    % ODE system
    dydt = zeros(6, 1);
    dydt(1) = (gZ - dZ) - (pZA + pZB) * Z;
    dydt(2) = (gA - dA) + pZA*Z - pAA1*A; 
    dydt(3) = (gA1 - dA1) + pAA1*A + uA1*(tA1)/(kA1+tA1)*B - uB*(tB)/(kB+tB)*A1;
    dydt(4) = (gB - dB) +  pZB*Z - uA1*(tA1)/(kA1+tA1)*B + uB*(tB)/(kB+tB)*A1;
    dydt(5) = -qAA * tA1;
    dydt(6) = -qB * tB;
end