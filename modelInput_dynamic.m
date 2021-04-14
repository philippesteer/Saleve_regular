function modelInput_dynamic()
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
parSPM.m1 = 0.5;                                                            % Stream power area exponent for rivers K.A^m.S^n.
parSPM.m2 = 0.24;                                                           % Stream power area exponent for colluvial valleys K.A^m.S^n.
parSPM.m3 = 0;                                                              % Stream power area exponent for hillslopes K.A^m.S^n.
parSPM.n  = 1;                                                              % Stream power slope exponent K.A^m.S^n (THE SOLUTION WORKS ONLY FOR N=1 (could be changed))
parSPM.Qc1 = 0;                                                             % River to colluvial valleys discharge transition
parSPM.Qc2 = 0;                                                             % Colluvial valleys to rivers discharge transition
% parSPM.Sc=30;parSPM.Qc2 = mean(mean((parSPM.U./(parSPM.K.*parSPM.Qc1.^(parSPM.m1-parSPM.m2))).^(1/parSPM.m2).*tand(parSPM.Sc).^-(parSPM.n./parSPM.m2))); % Colluvial valleys to rivers discharge transition

