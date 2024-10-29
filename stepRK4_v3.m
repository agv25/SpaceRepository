function [znext] = stepRK4_v3(t, z, dt)

%Runge Kutta 4 method that returns the next state vector, given inputs of
%time(vec), current state(matrix) and time step


%% Laying out coefficients A,B,C,D for averaging

dz = stateDeriv_v3(t, z);

A = dt * dz;


% Adjusting Coefficients based on stateDeriv inputs
B = dt * stateDeriv_v3(t + dt/2, z + A/2);


C = dt * stateDeriv_v3(t + dt/2, z + B/2);


D = dt * stateDeriv_v3(t + dt, z + C);

%% Combining into next state
znext = z +(A + 2*B + 2*C + D)/6;



