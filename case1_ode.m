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
    
    % Variables
    Z = y(1);
    A = y(2);
    A1 = y(3);
    B = y(4);
    
    % ODE system
    dydt = zeros(4, 1);
    dydt(1) = (gZ - dZ) - (pZA + pZB) * Z;
    dydt(2) = (gA - dA) + pZA*Z - pAA1*A; 
    dydt(3) = (gA1 - dA1) + pAA1*A;
    dydt(4) = (gB - dB) + pZB*Z;
end