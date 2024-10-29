function [znext] = stepRK4_Burn_v3(t, z, dt, index, OrbitRadii)

%Runge Kutta 4 method that returns the next Burn state vector, given inputs of
%time(vec), current state(matrix) and time step


%% Laying out coefficients A,B,C,D for averaging

dz = stateDeriv_Burn_v3(t, z, index, OrbitRadii);

A = dt * dz;


% Adjusting Coefficients based on stateDeriv inputs
B = dt * stateDeriv_Burn_v3( (t + dt/2), (z + A/2), index, OrbitRadii);


C = dt * stateDeriv_Burn_v3( (t + dt/2), (z + B/2), index, OrbitRadii);


D = dt * stateDeriv_Burn_v3( (t + dt), (z + C), index, OrbitRadii);

% Combining into next state
znext = z +(A + 2*B + 2*C + D)/6;



