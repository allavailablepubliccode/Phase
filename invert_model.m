function [LAP3, LAP4] = invert_model(Y, U, n, pE, pC, prec, varparam)
% outputs LAP3 and LAP4 are model inversions using first and second-order
% equations of motion, respectively

% Y is the input data
% U is the driving input
% n is the number of gradients
% pE are model parameters (prior means)
% pC are prior variabce
% prec are precisions (noise)
% varparam is the scaling of prior variance

% 1st order, with cross-connections

% initial states
x.z = zeros(n,1);

% observation function
g   = @(x,v,P) P.s*x.z;

% equation of motion (Eq. (1) in text)
f   = @(x,v,P) P.A*x.z + P.C*v;

% model parameters
DEM.M(1).pE        = pE;

% prior variace
DEM.M(1).pC        = diag(spm_vec(pC))*varparam;

% initial conditions
DEM.M(1).x         = x;

% observation function
DEM.M(1).g         = g;

% equation of motion
DEM.M(1).f         = f;

% precision of observation noise
DEM.M(1).V         = exp(prec);

% precision of state noise
DEM.M(1).W         = exp(prec);

% initial causes
DEM.M(2).v         = 0;

% precision of exogenous causes
DEM.M(2).V         = exp(prec);

% driving input
DEM.U              = U;

% data
DEM.Y              = Y;

% LAP1 = spm_LAP(DEM);

% second order, with cross-connections
clear x
x.y     = zeros(n,1);
x.z     = zeros(n,1);
g       = @(x,v,P) P.s*x.y;
% Eq. (2) in text
f       = @(x,v,P) [x.z; P.A*x.y + P.C*v];

DEM.M(1).x         = x;
DEM.M(1).g         = g;
DEM.M(1).f         = f;

% LAP2 = spm_LAP(DEM);

% 1st order, no cross-connections
pC.A = pC.A .* eye(n);

clear x
x.z = zeros(n,1);
g   = @(x,v,P) P.s*x.z;
f   = @(x,v,P) P.A*x.z + P.C*v;

DEM.M(1).pC        = diag(spm_vec(pC))*varparam;
DEM.M(1).x         = x;
DEM.M(1).g         = g;
DEM.M(1).f         = f;

LAP3 = spm_LAP(DEM);

% second order, no cross-connections
clear x
x.y     = zeros(n,1);
x.z     = zeros(n,1);
g       = @(x,v,P) P.s*x.y;
f       = @(x,v,P) [x.z; P.A*x.y + P.C*v];

DEM.M(1).x         = x;
DEM.M(1).g         = g;
DEM.M(1).f         = f;

LAP4 = spm_LAP(DEM);

