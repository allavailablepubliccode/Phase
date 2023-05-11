function [LAP1, LAP2, LAP3, LAP4] = invert_model(Y, U, n, pE, pC, prec)

% 1st order, with cross-connections
x.z = zeros(n,1);
g   = @(x,v,P) P.s*x.z;
f   = @(x,v,P) P.A*x.z + P.C*v;

DEM.M(1).pE        = pE;
DEM.M(1).pC        = diag(spm_vec(pC));
DEM.M(1).x         = x;
DEM.M(1).g         = g;
DEM.M(1).f         = f;
DEM.M(1).V         = exp(prec);
DEM.M(1).W         = exp(prec);
DEM.M(2).v         = 0;
DEM.M(2).V         = exp(prec);
DEM.U              = U;
DEM.Y              = Y;

LAP1 = spm_LAP(DEM);

% second order, with cross-connections
clear x
x.y     = zeros(n,1);
x.z     = zeros(n,1);
g       = @(x,v,P) P.s*x.y;
f       = @(x,v,P) [x.z; P.A*x.y + P.C*v];

DEM.M(1).x         = x;
DEM.M(1).g         = g;
DEM.M(1).f         = f;

LAP2 = spm_LAP(DEM);

% 1st order, no cross-connections
pC.A = pC.A .* eye(n);

clear x
x.z = zeros(n,1);
g   = @(x,v,P) P.s*x.z;
f   = @(x,v,P) P.A*x.z + P.C*v;

DEM.M(1).pC        = diag(spm_vec(pC));
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

