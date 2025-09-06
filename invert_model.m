function [LAP3, LAP4] = invert_model(Y, U, n, pE, pC, prec, varparam)
%==========================================================================
% invert_model.m
%
% Description:
%   Performs model inversion for synthetic dynamic causal models using
%   both first- and second-order equations of motion, with and without
%   cross-connections. Returns Laplace-approximated posterior estimates.
%
% Inputs:
%   Y         - Observed data [n x T]
%   U         - External (driving) input [n x T]
%   n         - Number of sources (nodes)
%   pE        - Prior expectations of model parameters (struct)
%   pC        - Prior covariances of model parameters (struct)
%   prec      - Precision of noise (scalar, log precision)
%   varparam  - Scaling factor for prior variance
%
% Outputs:
%   LAP3      - Laplace estimate (first-order, no cross-connections)
%   LAP4      - Laplace estimate (second-order, no cross-connections)
%
% Requires:
%   - SPM12 (specifically, spm_LAP and spm_vec)
%
%==========================================================================

%==========================================================================
% First-order model with cross-connections (Eq. 1 in paper)
%==========================================================================
x.z = zeros(n,1);                          % Initial latent states
g    = @(x,v,P) P.s * x.z;                 % Observation function
f    = @(x,v,P) P.A * x.z + P.C * v;       % Evolution function

% Build DEM model structure
DEM.M(1).x   = x;
DEM.M(1).f   = f;
DEM.M(1).g   = g;
DEM.M(1).pE  = pE;
DEM.M(1).pC  = diag(spm_vec(pC)) * varparam;
DEM.M(1).V   = exp(prec);                 % Observation noise precision
DEM.M(1).W   = exp(prec);                 % State noise precision

DEM.M(2).v   = 0;                         % Initial exogenous causes
DEM.M(2).V   = exp(prec);                % Cause noise precision

DEM.U       = U;                         % External input
DEM.Y       = Y;                         % Observed data

% This model is skipped from output, could be retained as LAP1 if needed
% LAP1 = spm_LAP(DEM);

%==========================================================================
% Second-order model with cross-connections (Eq. 2 in paper)
%==========================================================================
x.y = zeros(n,1);                          % Position-like states
x.z = zeros(n,1);                          % Velocity-like states
g    = @(x,v,P) P.s * x.y;                 % Observe only position
f    = @(x,v,P) [x.z; P.A * x.y + P.C * v];% Second-order dynamics

DEM.M(1).x  = x;
DEM.M(1).f  = f;
DEM.M(1).g  = g;

% This model is also skipped from output, could be retained as LAP2
% LAP2 = spm_LAP(DEM);

%==========================================================================
% First-order model WITHOUT cross-connections
%==========================================================================
pC.A = pC.A .* eye(n);                    % Remove off-diagonal prior covariances

x.z = zeros(n,1);
g   = @(x,v,P) P.s * x.z;
f   = @(x,v,P) P.A * x.z + P.C * v;

DEM.M(1).x  = x;
DEM.M(1).f  = f;
DEM.M(1).g  = g;
DEM.M(1).pC = diag(spm_vec(pC)) * varparam;

LAP3 = spm_LAP(DEM);                     % Perform Laplace inversion

%==========================================================================
% Second-order model WITHOUT cross-connections
%==========================================================================
x.y = zeros(n,1);
x.z = zeros(n,1);
g   = @(x,v,P) P.s * x.y;
f   = @(x,v,P) [x.z; P.A * x.y + P.C * v];

DEM.M(1).x = x;
DEM.M(1).f = f;
DEM.M(1).g = g;

LAP4 = spm_LAP(DEM);                     % Perform Laplace inversion
