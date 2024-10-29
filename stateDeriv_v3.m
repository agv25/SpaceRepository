function [dz] = stateDeriv_v3(t, z)
% calculates the state derivative for aerocapture of Spacecraft
% around Venus. Ignoring Drag for now.

% Takes inputs of z which includes the position and velocity in x and y
% and an input of time
% Outputs dz vector which includes the derivative state (acceleration=dz2)

%% Defining variables:
planet
Spacecraft

% Altitude of spacecraft (m)
r = (sqrt(z(1,end,1)^2 + z(1,end,2)^2))-6051.8e3;

% Atmospheric density kg/m^3;
rho = profileVenus(r); 


%% Defining angle of spacecraft characteristics relative to Venus

% Defining angle of Spacecraft's position relative to x axis
theta = atan2( z(1,end,2), z(1,end,1) )+ (2*pi)*(z(1,end,2)<0);

% Defining angle of Spacecraft's velocity relative to x axis
phi = atan2( z(2,end,2), z(2,end,1) ) + (2*pi)*(z(2,end,2)<0);



combVel = hypot(z(2,end,1),z(2,end,2));


%% Describing direction dependent ODEs

xacc = (cos(theta)*((-G*Mv)/(z(1,end,1)^2 + z(1,end,2)^2)))+(cos(phi-pi)*(.5*A*Cd*rho*combVel^2)/Ms);
yacc = (sin(theta)*((-G*Mv)/(z(1,end,1)^2 + z(1,end,2)^2)))+(sin(phi-pi)*(.5*A*Cd*rho*combVel^2)/Ms);


%% Defining state derivatives, with [x y] convention

dz1x = z(2, end, 1);  % velocity x


dz2x = xacc;        %acceleration x


dz1y = z(2, end, 2);  % velocity y


dz2y = yacc;        % acceleration y


%% dz is a matrix of single column, 2 pages(representing 2 axes) holding the
% state derivatives

dz(:,:,1) = [dz1x ; dz2x];
dz(:,:,2) = [dz1y ; dz2y];




