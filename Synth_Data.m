%==========================================================================
% generate_synthetic_DCM_datasets.m
%
% Description:
%   This script generates synthetic datasets based on first-order and
%   second-order dynamic causal models (DCMs) using the SPM12 DEM toolbox.
%
%   Each iteration creates randomized parameters for internal and external
%   connections, generates sinusoidal inputs, simulates dynamics, and
%   saves the resulting dataset after model inversion.
%
% Requirements:
%   - SPM12 installed and on MATLAB path
%   - 'invert_model' function defined in MATLAB path
%
% Output:
%   Saves one .mat file per dataset with the inferred parameters for both
%   1st- and 2nd-order models.
%
%==========================================================================

clear; clc; close all;

%-------------------------------
% User Configuration Parameters
%-------------------------------

dir_out = '';           % Output directory for saved .mat files (default = current folder)
num_datasets = 1000;    % Number of synthetic datasets to generate
n = 3;                  % Number of observed sources
N = 405;                % Number of time points per dataset

% Loop through and generate each dataset
for ii = 1:num_datasets

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: Generate 1st-Order System
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define initial state vector for each node
    clear x
    x.z = zeros(n,1);    % Latent neural states

    % Define observation function: maps latent states to observed outputs
    g = @(x,v,P) x.z;

    % Define evolution function: linear first-order dynamics
    f = @(x,v,P) P.A * x.z + P.C * v;

    % Initialize DCM parameters (P)
    clear P
    P.A = zeros(n);                        % Internal (self) connections
    r = 0.3 + (0.7 - 0.3) * rand(1,n);     % Random negative self-connections
    P.A(logical(eye(n))) = -r;            % Set diagonals of A

    P.C = 0.3 + (0.7 - 0.3) * rand(1,n);   % External connections
    P.C = P.C .* eye(n);                  % Diagonal input connections only

    % Define sinusoidal input U with random amplitude, frequency, phase
    U = rand(n,1) .* sin(rand(n,1) .* (1:N) + rand);

    % Construct SPM DEM model structure (1st order)
    clear M
    M(1).x  = x;            % Initial states
    M(1).f  = f;            % Evolution function
    M(1).g  = g;            % Observation function
    M(1).pE = P;            % Parameter estimates
    M(1).V  = exp(4)*ones(n,1);  % Observation noise precision
    M(1).W  = exp(4)*ones(n,1);  % State noise precision

    % Second-level model (exogenous causes)
    M(2).v = zeros(n,1);
    M(2).V = exp(4)*ones(n,1);

    % Generate synthetic data using SPM’s DEM engine
    DEM_1 = spm_DEM_generate(M, U, P);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: Generate 2nd-Order System
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear x
    x.y = zeros(n,1);   % Position-like state
    x.z = zeros(n,1);   % Velocity-like state

    % Observation: only observe position
    g = @(x,v,P) x.y;

    % Evolution: second-order dynamics as 2 coupled 1st-order equations
    f = @(x,v,P) [x.z; P.A * x.y + P.C * v];

    % Reuse same parameters and input
    M(1).x = x;
    M(1).f = f;
    M(1).g = g;

    DEM_2 = spm_DEM_generate(M, U, P);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 3: Invert Models to Estimate Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define prior expectations (mean)
    clear pE
    pE.A = zeros(n);
    pE.C = zeros(n,1);
    pE.s = 1;

    % Define prior covariances (uncertainty)
    clear pC
    pC.A = eye(n);
    pC.C = ones(n,1);
    pC.s = 1;

    % No inputs during inversion
    U_null = zeros(1,N);

    % Variational precision (inverse noise)
    prec = 4;

    % Variational mode (1 = Laplace, 0 = fixed)
    varparam = 1;

    % Invert 1st- and 2nd-order models
    [LAP_g1_1, LAP_g1_2] = invert_model(DEM_1.Y, U_null, n, pE, pC, prec, varparam);
    [LAP_g2_1, LAP_g2_2] = invert_model(DEM_2.Y, U_null, n, pE, pC, prec, varparam);

    % Save results
    save(fullfile(dir_out, sprintf('%04d.mat', ii)), ...
         'LAP_g1_1','LAP_g1_2','LAP_g2_1','LAP_g2_2');

end
