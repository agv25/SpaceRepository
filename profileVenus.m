function [rho] = profileVenus(z)
%% profileVenus  Atmospheric density profile model for Venus 
%
% [RHO] = profileVenus(Z) outputs a vector RHO of modelled density values 
% at the corresponding altitude values in the input vector Z. Input 
% altitudes must be specified in meters and output densities are given in 
% kg/m^3.
% 
% For altitudes 0-100km, empirical data taken from:
% V.I. Moroz, "The atmosphere of Venus", Space Science Reviews,
% vol. 29, no. 1, 1981, pp. 3-127
% 
% For altitudes 130-190km, model data taken from:
% https://sci.esa.int/web/venus-express/-/57736-density-profiles-of-venus-polar-atmosphere
% 
% Alan Hunter, University of Bath, Oct 2020.


%% Atmospheric model data for Venus

% Altitude samples
z0 = [...
    1e-3
    10
    20
    30
    40
    50
    60
    70
    75
    80
    85
    90
    95
    100
    130
    140
    150
    160
    170
    180
    190
    ] * 1e3; % m

% Density values
rho0 = [...
    64.8
    37.3
    20.1
    9.95
    4.28
    1.53
    4.41e-1
    7.68e-2
    2.90e-2
    1.03e-2
    3.50e-3
    1.07e-3
    3.03e-4
    8.25e-5
    1e-8
    1e-9
    7e-11
    1e-11
    3e-12
    1e-12
    5e-13
    ]; % kg/m^3


%% Interpolate data at requested altitude(s)

% Interpolate data logarithmically for accuracy
z0log = log(z0);
rho0log = log(rho0);

zlog = log(z);

% Handle negative altitudes
zlog(z<=0) = 0;

% Interpolate
rholog = interp1(z0log,rho0log,zlog,'linear',-Inf);
rho = exp(rholog);

% Assign infinite density for negative altitudes
rho(z<=0) = Inf;