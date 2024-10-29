function [dz] = stateDeriv_Burn_v3(t, z, index, OrbitRadii)
% calculates the state derivative for aerocapture of Spacecraft
% around Venus. Ignoring Drag for now.

% Takes inputs of z which includes the position and velocity in x and y
% and an input of time
% Outputs dz vector which includes the derivative state (acceleration=dz2)

%% Defining variables:
planet
Spacecraft


% Gravitational parameter of Venus
mu = G*Mv;

% Altitude of spacecraft (m)
r = (sqrt(z(1,end,1)^2 + z(1,end,2)^2))-6051.8e3;

% Atmospheric density (kg/m^3)
rho = profileVenus(r);

%% Extracting Orbit Radii:
r1 = OrbitRadii(1);
r2 = OrbitRadii(2);

%% Defining angle of spacecraft characteristics relative to Venus

% Defining angle of Spacecraft's position relative to x axis
theta = atan2( z(1,end,2), z(1,end,1) )+ (2*pi)*(z(1,end,2)<0);

% Defining angle of Spacecraft's velocity relative to x axis
phi = atan2( z(2,end,2), z(2,end,1) ) + (2*pi)*(z(2,end,2)<0);

% Defining magnitude of Velocity
combVel = hypot(z(2,end,1),z(2,end,2));



%% Defining Thrust force to go into circular orbit

vCircular = sqrt(mu/r2);
vAlpha = vCircular*(1-sqrt((2*r1)/(r1+r2)));

%Thrust burn for single time step. Time being when reached apoapsis.
if( t ~= (index+2000))
    
    aBurn = 0;
    
else
    
    aBurn = (vAlpha)/0.32;     %Applying a correction factor to vAlpha
    
end

%% Describing direction dependent ODEs

xacc = (cos(theta)*((-G*Mv)/(z(1,end,1)^2 + z(1,end,2)^2)))+(cos(phi-pi)*(.5*A*Cd*rho*combVel^2)/Ms)+(cos((2*pi)-phi)*aBurn);
yacc = (sin(theta)*((-G*Mv)/(z(1,end,1)^2 + z(1,end,2)^2)))+(sin(phi-pi)*(.5*A*Cd*rho*combVel^2)/Ms)+(sin((2*pi)-phi)*aBurn);


%% Defining state derivatives, with [x y] convention

dz1x = z(2, end, 1);  % velocity x


dz2x = xacc;        %acceleration x


dz1y = z(2, end, 2);  % velocity y


dz2y = yacc;        % acceleration y


%% dz is a matrix of single column, 2 pages(representing 2 axes) holding the
% state derivatives

dz(:,:,1) = [dz1x ; dz2x];
dz(:,:,2) = [dz1y ; dz2y];




