clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate 1st order synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% enter output directory:
dir_out = '';

% number of datasets to generate
for ii = 1:1000

    % number of sources
    n       = 3;

    % initial states
    clear x
    x.z     = zeros(n,1);

    % observation function
    g       = @(x,v,P) x.z;

    % equation of motion
    f       = @(x,v,P) P.A*x.z + P.C*v;

    % number of time points
    N       = 405;

    % randomized internal self-connections
    clear P
    P.A     = zeros(n);
    r = 0.3 + (0.7-0.3)*rand(1,n);
    P.A(logical(eye(size(P.A)))) = -r;

    % randomized external connections
    P.C     = 0.3 + (0.7-0.3)*rand(1,n);
    P.C     = P.C .* eye(n);

    % randomized sinusoidal driving inputs
    U       = rand(n,1).*sin(rand(n,1).*(1:N) + rand);

    % first level state space model
    %--------------------------------------------------------------------------
    clear M
    M(1).x  = x;                             % initial states
    M(1).f  = f;                             % equations of motion
    M(1).g  = g;                             % observation mapping
    M(1).pE = P;                             % model parameters
    M(1).V  = exp(4)*ones(n,1);                       % precision of observation noise
    M(1).W  = exp(4)*ones(n,1);                       % precision of state noise

    % second level  causes or exogenous forcing term
    %--------------------------------------------------------------------------
    M(2).v  = zeros(n,1);                             % initial causes
    M(2).V  = exp(4)*ones(n,1);                       % precision of exogenous causes

    DEM_1     = spm_DEM_generate(M,U,P);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate 2nd order synthetic data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear x
    x.y     = zeros(n,1);
    x.z     = zeros(n,1);
    g       = @(x,v,P) x.y;
    f       = @(x,v,P) [x.z; P.A*x.y + P.C*v];

    M(1).x  = x;                             % initial states
    M(1).f  = f;                             % equations of motion
    M(1).g  = g;                             % observation mapping

    DEM_2     = spm_DEM_generate(M,U,P);

    clear pE
    pE.A = zeros(n);
    pE.C = zeros(n,1);
    pE.s = 1;

    clear pC
    pC.A = eye(n);
    pC.C = ones(n,1);
    pC.s = 1;

    U = zeros(1,N);

    prec = 4;

    varparam = 1;

    [LAP_g1_1, LAP_g1_2] = invert_model(DEM_1.Y, U, n, pE, pC, prec, varparam);
    [LAP_g2_1, LAP_g2_2] = invert_model(DEM_2.Y, U, n, pE, pC, prec, varparam);

    save([dir_out num2str(ii) '.mat'],'LAP_g1_1','LAP_g1_2','LAP_g2_1','LAP_g2_2')

end