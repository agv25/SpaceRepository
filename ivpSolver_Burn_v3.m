function [t, z] = ivpSolver_Burn_v3(t0, z0, dt, tend, index, OrbitRadii)

% Given initial Values and a time range, the function will return
% a matrix of Burn state vectors Z, as well as a time vector t
% Solves Initial Value Problem/Boundary problem


% Setting Time as a Vector for while loop.

t(1) = t0;

%% Assigning x to one page, y to another page. Burn state variables will get put
% into columns on each page.

z(:,1,1) = z0(:,:,1);
z(:,1,2) = z0(:,:,2);

%% Iterating updates for time until tend by time step dt
n = 1;
while t(n)<= tend
    % Increment time vector by one time step
    t(n+1) = t(n) + dt;
    
    % Apply RK4 method for one time step
    z(:, n+1, :) = stepRK4_Burn_v3(t(n), z(:, n, :), dt, index, OrbitRadii);
    
    
    n = n+1;
end