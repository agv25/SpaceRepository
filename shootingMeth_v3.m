function [finalAlpha] = shootingMeth_v3(alpha0, alpha1)
% function to determine the alpha value for spacecraft insertion height
% above the surface of Venus. Implementing Shooting Method over a limited
% range.

%% Controlling range of alphas due to limitation of Linear Shooting Method
%alpha vector holding the current two alpha values being examined
alpha = [alpha0 ; alpha1];

% Sorting in ascending order
alpha = sort(alpha);

% Applying range limits.
if(alpha(2)> 1.33e6)
    alpha(2) = 1.33e6;
elseif(alpha(1)<1.328e6)
    alpha(1) = 1.328e6;
end

%%  Initiliasing variables

% m - wanted value of H (apoapsis distance to surface of Venus)
Hwanted = 1200e3;

% Initial State Vector components to be used

x0 = 10000e3; % m - initial x position of Spacecraft
vx0 = -11000; % m/s - initial x velocity of Spacecraft


% Setting up time parameters
timeParam

% Radius of Venus (m)
R = 6051.8e3;

% State vectors in x direction
z0(:,:,1) = [x0 ; vx0];




%plotting Venus if wanting to plot each iteration trajectory
plotCircle
hold on;

%% Use IVP solver to find solution for 1st guess

% Initial altitude in y axis of spacecraft(first guess)
y01 = R + alpha(1);

% initial state vectors in y
z0(:,:,2) = [y01 ; 0];

% Solve IVP for this guess

[t, z] = ivpSolver_v3(t0, z0, dt, tend);

%Finding Apoapsis point of trajectory

radius = hypot(z(1,:,1), z(1,:,2));

adjustedRadii1 = radius(2000:end)-R;


H1 = max(adjustedRadii1);



%% While Loop to continuously improve alpha shot taken
epsilon = 0;

while ( (epsilon < 0.9999) || (epsilon > 1.0001) )
    
    
    % Running for second guess
    y02 = R + alpha(2); % m
    
    % state vector in y direction
    z0(:,:,2) = [y02 ; 0];
    
    
    % Solve IVP for this guess
    
    [t, z] = ivpSolver_v3(t0, z0, dt, tend);
    
    
    %Finding Apoapsis point of trajectory
    
    radius = hypot(z(1,:,1), z(1,:,2));
    
    adjustedRadii = radius(2000:end)-R;
    
    H2 = max(adjustedRadii);
    
    %% Error Ratio based on linear relationship between H and alpha
    
    ratio = (Hwanted - H1) / (H2 - H1);
    
    
    
    % Assigning new alpha based on error ratio and previous guess
    alphaNew = (alpha(2)-alpha(1))*abs(ratio)+alpha(1);
    
    % Shifting values in order to keep the two most up to date in loop.
    alpha(1) = alpha(2);
    alpha(2) = abs(alphaNew);
    
    H1 = H2;
    
    
    
    
    
    
    
    % Determining error ratio to tell loop to keep iterating or stop.
    
    epsilon = H1/Hwanted;
    
    % To plot trajectory shots
    hold on;
    figure(1)
    title('Shooting Method Trajectory Plot');
    hold on
    plot( z(1,:,1), z(1,:,2), 'k');
    legend('Shooting Method Trajectories');
    ylabel('y Co-ordinate (m)');
    xlabel('x Co-ordinate (m)');
    axis equal;
    drawnow;
    
end

% Output final value alpha

finalAlpha = alpha(2);



