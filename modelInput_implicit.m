function modelInput_implicit()
% Read model parameters and stock them in global variable

global  parSPM;

% Time parameters - dynamic model
parSPM.T  = 0.5e6.*365;                                                     % Model duration (day)
parSPM.dt = 2000.*365;                                                       % Time step (day)
parSPM.t  = [parSPM.dt:parSPM.dt:parSPM.T];                                 % Time vector
parSPM.nt = numel(parSPM.t);                                                % Number of time steps

% Space parameters
parSPM.dx = 50;                                                            % Number of model points
parSPM.L  = 10000;                                                          % Model horizontal length

% Uplift and precipitation rates
parSPM.U  = 0.01./365;                                                      % Uplift rate (m/day)
parSPM.P  = 5./365;                                                         % Precipitation (m/day)

% Stream Power erosion law
parSPM.K  = 1e-6;                                                           % Stream power efficiency for rivers K.A^m.S^n.
parSPM.m = 0.5;                                                            % Stream power area exponent for rivers K.A^m.S^n.