%% Assignment 2: Aerocapture of Spacecraft around Venus. Modelling using Runge-Kutta 4.
% Modelling the system using RK4.

%% Setting Up Variables

% Accessing Planet(Venus) information.
planet
%Accessing time Parameters
timeParam


% Shooting Method Initial Guesses

alpha1 = 1333000; % initial guess (m)
alpha0 = 1328000; % second initial guess (m)  

% Using 2 guesses, optimised Alpha is the output for 1.2e6 = H.
finalAlpha = shootingMeth_v3(alpha0, alpha1);

%% Setting up Parameters for ivpSolver

x0 = 10000e3;           % m - initial x position of Spacecraft
vx0 = -11000;           % m/s - initial x velocity of Spacecraft

y0 = R + finalAlpha;    % m - initial y position of Spacecraft
vy0 = 0;                % m/s - initial y velocity of Spacecraft

% Combining initial State Vectors in a single 3D array.
z0(:,:,1) = [x0 ; vx0];
z0(:,:,2) = [y0 ; vy0];



%% Running ivpSolver for z output matrix of aerocapture trajectory

[t, z] = ivpSolver_v3(t0, z0, dt, tend);


%% Finding H (apoapsis) & h (perigee)

% Taking the radius of the spacecraft from the origin of the planet (m)
radius = hypot(z(1,:,1), z(1,:,2));

% adjustedRadii is a vector of altitudes after 1st part of trajectory (m)
adjustedRadii = radius(2000:end)-R;

% Copy of state matrix in order to clip trajectory approaching to planet
zCopy = [z(1,:,1) ; z(1,:,2)];
zCopy(:, 1:2000) = [];

% Taking the maximum altitude of the spacecraft as the apoapsis time index
% with respect to trajectory time steps
[H, index] = max(adjustedRadii);

% Creating a normal line from the apoapsis to the origin
apoapsis = [0 zCopy(1, index) ; 0 zCopy(2,index)];

% Mirror process to find perigee, but during approach to planet
[h, index2] = min(radius(1:2000)-R);
perigee = [0 z(1, index2, 1) ; 0 z(1, index2, 2)];

% Printing H and h to command window
Hprint = ['H = ', num2str(H)];
disp(Hprint);
hprint = ['h = ', num2str(h)];
disp(hprint)   

%% Creating combined magnitude vectors for Velocity and Acceleration for plotting

% Differentiating velocity to find acceleration at each time step
ddzX = diff(z(2,:,1)) / dt;
ddzY = diff(z(2,:,2)) / dt;

% Taking the magnitude of the sum of Acceleration components
combAcc = hypot(ddzX, ddzY);

% Taking the magnitude of the sum of Velocity components
combVel = hypot(z(2,:,1),z(2,:,2));


%% Plot information on basic trajectory


% Plotting on top of plot from Shooting Method plot with finalAlpha
figure(1)
hold on
plotCircle;                                    % plot Venus
plot(z(1,:,1),z(1,:,2),'g');                   % plot position x&y
ylabel('y Co-ordinate (m)');
xlabel('x Co-ordinate (m)');
legend('Venus','Shooting Method Trajectories');
axis equal;

% Plot of final aerocapture trajectory with detail
figure(2)
title('Final Alpha Trajectory Plot');
hold on
plotCircle;                                    % plot Venus
plot(z(1,:,1),z(1,:,2),'g');                   % plot position x&y
plot (apoapsis(1,:),apoapsis(2,:),'b');        % plot Apoapsis
plot (perigee(1,:),perigee(2,:),'m');          % plot Perigee
ylabel('y Co-ordinate (m)');    
xlabel('x Co-ordinate (m)');
legend('Venus','Final Alpha Trajectory', 'Normal to H', 'Normal to h');
axis equal;

% Plotting Resultant Acceleration over time
figure(3)
title('Acceleration vs Time Plot')
hold on
plot(t(2:end), combAcc, 'r');          
ylabel('Resultant Acceleration (m/s^2)');
xlabel ('Time (s)');

% Plotting Resultant Velocity over time
figure(4)
title('Velocity vs Time Plot');
hold on
plot(t, combVel, 'b');                  
ylabel('Resultant Velocity (m/s)');
xlabel('Time (s)');
axis auto;

axis equal;


%% Calculating trajectory with a burn at the point of apoapsis
% Purpose of burn is to place spacecraft in circular orbit of radius H

% Setting up vector of orbital radii for Hohmann transfer
% Based on Apoapsis and Perigee
OrbitRadii = [(h+R), (H+R) ];

% Calculating new trajectory
[t, zBurn] = ivpSolver_Burn_v3(t0, z0, dt, tend, index, OrbitRadii);

% New radius vector for the newly calculated trajectory
radiusBurn = hypot(zBurn(1,:,1), zBurn(1,:,2));

%% Creating combined magnitude vectors for Burn Velocity & Burn Acceleration
% Can be used to plot Burn trajectory characteristics

% Differentiating Burn velocity to find Burn acceleration at each time step
ddzXburn = diff(zBurn(2,:,1))/ dt;
ddzYburn = diff(zBurn(2,:,1))/ dt;

% Taking the magnitude of the sum of Acceleration components
combAccBurn = hypot(ddzXburn, ddzYburn);

% Taking the magnitude of the sum of Velocity components
combVelBurn = hypot(zBurn(2,:,1),zBurn(2,:,2));

%% Plot both trajectories in a 3D environment

% To orbit around equator, setting Z-axis positions and velocities to 0.
z(:,:,3) = 0;
zBurn(:,:,3) = 0;

% Making a unit sphere and scaling it to match Venus's Radius
[a, b, c] = sphere(100);
a = a*R;
b = b*R;
c = c*R;

% Making a unit sphere and scaling it to match the outer edge of the
% atmosphere
[aAtm, bAtm, cAtm] = sphere(100);
aAtm = aAtm*(R+190e3);
bAtm = bAtm*(R+190e3);
cAtm = cAtm*(R+190e3);


% Extracting Basic Trajectory Positions
xOut = z(1,:,1);
yOut = z(1,:,2);
zOut = z(1,:,3);

% Extracting Burn Trajectory Positions
xBurnOut = zBurn(1,:,1);
yBurnOut = zBurn(1,:,2);
zBurnOut = zBurn(1,:,3);


figure(5)
set(figure(5), 'color', '#AAAAAA')  % Grey background for better contrast
set(gca, 'color', '#AAAAAA')        % Grey background for better contrast
grid on                             % Grid to help with depth perception

xlabel('X Co-ordinate (m)');
ylabel('Y Co-ordinate (m)');
zlabel('Z Co-ordinate (m)');

hold on



%Plotting Venus and changing sphere attributes for better visualisation
q = surf(a,b,c, 'EdgeColor', 'none');
axis equal
colormap('copper')

shading interp

%Plotting Venus's atmosphere as 
atmosphere = surf(aAtm, bAtm, cAtm, 'EdgeColor', '#d9d9d9');
set(atmosphere, 'FaceAlpha', 0.2)
shading interp



% Plotting the trajectories

plot3(xBurnOut,yBurnOut,zBurnOut, 'b-', 'Linewidth', 0.6)

plot3(xOut, yOut, zOut, 'r--', 'Linewidth', 0.6)

legend('Venus','Venus Atmosphere', 'Trajectory with Burn','Trajectory without Burn');