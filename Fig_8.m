% generate Fig. 8
clear;clc;close all;

% specify states and parameters
%==========================================================================
N           = 110;      % 10x max. no. of free params.

% model states (initial conditions)
%--------------------------------------------------------------------------
x.y         = 1;
x.z         = 1;

% Intrinsic connectivity
%--------------------------------------------------------------------------
P.A         = -1;

% observation function (to generate timeseries)
%--------------------------------------------------------------------------
g           = @(x,v,P) P.y*x.y + P.z*x.z;

% equations of motion (1st or 2nd order)
%--------------------------------------------------------------------------
f           = @(x,v,P) [x.z; P.A*(P.y*x.y + P.z*x.z)];

% causes or exogenous input
%--------------------------------------------------------------------------
U           = zeros(1,N);

% first level state space model
%--------------------------------------------------------------------------
M(1).x      = x;
M(1).f      = f;
M(1).g      = g;

% second level (causes or exogenous forcing term)
%--------------------------------------------------------------------------
M(2).v      = 0;

switchord   = [0 1; 1 0];
dt          = [1/32 1/4];
errs        = [24 12];

model{1}    = '1 comp, ord 1';
model{2}    = '1 comp, ord 2';

% create data with known parameters (P)
%==========================================================================
k           = 0;
f           = figure;
f.Position  = [650 650 650 650];
for cc      = 1:numel(model)
    k       = k + 1;

    disp(['Generating synthetic data, ' model{cc}])

    % model parameters
    %----------------------------------------------------------------------
    P.y         = switchord(cc,1);
    P.z         = switchord(cc,2);

    M(1).pE     = P;
    M(1).E.dt   = dt(cc);

    M(1).V      = exp(errs(cc));
    M(1).W      = exp(errs(cc));
    M(2).V      = exp(errs(cc));

    D           = spm_DEM_generate(M,U,P);

    subplot(4,4,k)
    plot(D.Y')
    title([model{cc} ', tser'])
    axis tight
    drawnow

    k           = k + 1;
    subplot(4,4,k)

    if cc == 1
        plot(D.Y',zeros(size(D.Y')),'.')
        grid on
        title([model{cc} ', ph sp'])
        axis tight
        drawnow
    else
        plot(D.Y(1:end-1)',diff(D.Y),'.')
        grid on
        title([model{cc} ', ph sp'])
        axis tight
        drawnow
    end

    % Now try to recover model parameters from data features
    %======================================================================

    % initialization of priors over parameters
    %----------------------------------------------------------------------
    pE.A        = -1/4;     % prior parameters
    pE.y        = 0;
    pE.z        = 0;

    pC.A        = 1;       % prior variance
    pC.y        = 1;
    pC.z        = 1;

    % place new model states and priors in generative model
    %----------------------------------------------------------------------
    D.M(1).pE   = pE;
    D.M(1).pC   = diag(spm_vec(pC));

    disp(['mod inv, ~2 mins, ' model{cc}])

    % Inversion using generalised filtering
    %======================================================================
    LAP         = spm_LAP(D);

    % use Bayesian model reduction to test different hypotheses
    %======================================================================

    % apply precise shrinkage priors to of diagonal coupling elements
    %----------------------------------------------------------------------
    clear PC
    PC{1}     = pC; PC{1}.y = 0;
    PC{2}     = pC; PC{2}.z = 0;

    %  evaluate the evidence for these new models or prior constraints
    %----------------------------------------------------------------------
    qE          = LAP.qP.P{1};
    qC          = LAP.qP.C;
    pE2         = LAP.M(1).pE;
    pC2         = LAP.M(1).pC;
    clear F
    for m       = 1:numel(PC)
        rC      = diag(spm_vec(PC{m}));
        F(m)    = spm_log_evidence(qE,qC,pE2,pC2,pE2,rC);
    end

    % report marginal log likelihood or evidence
    %----------------------------------------------------------------------
    F           = F - min(F);

    k = k + 1;
    subplot(4,4,k)
    bar(F,'c')
    title([model{cc} ', log ev.'])
    axis square,
    box off

    k = k + 1;
    subplot(4,4,k)
    bar(spm_softmax(F(:)),'c')
    title([model{cc} ', prob.'])
    axis square,
    box off,
    drawnow

end

x.y             = ones(3,1);
x.z             = ones(3,1);
M(1).x          = x;

A(:,:,1)        = [-2 -1/2 -1/2;
    1/4 -1 -1/4;
    1/2  1/2 -1/2];

A(:,:,2)        = [-1/8 -1/32 -1/32;
    1/32 -1/8 -1/32;
    1/32  1/32 -1/8];

M(1).V          = exp(20);
M(1).W          = exp(20);
M(2).V          = exp(20);

dt              = [1/32 1/2];

model{1}    = '3 comp, ord 1';
model{2}    = '3 comp, ord 2';

for cc          = 1:numel(model)
    k           = k + 1;

    disp(['Generating synthetic data, ' model{cc}])

    % model parameters
    %----------------------------------------------------------------------
    P.A         = A(:,:,cc);
    P.y         = switchord(cc,1);
    P.z         = switchord(cc,2);

    M(1).pE     = P;
    M(1).E.dt   = dt(cc);

    D           = spm_DEM_generate(M,U,P);

    subplot(4,4,k)
    plot(D.Y')
    title([model{cc} ', tser'])
    axis tight
    drawnow

    k           = k + 1;
    subplot(4,4,k)

    if cc == 1
        plot3(D.Y(1,:)',D.Y(2,:)',D.Y(3,:)','.')
        grid on
        view(-16,28)
        title([model{cc} ', ph sp'])
        axis tight
        drawnow
    else
        title([model{cc} ', ph sp, 6-D'])
        drawnow
    end

    % Now try to recover model parameters from data features
    %==========================================================================

    % initialization of priors over parameters
    %--------------------------------------------------------------------------
    pE.A        = -eye(3)/4;     % prior parameters
    pE.y        = 0;
    pE.z        = 0;

    pC.A        = ones(3);       % prior variance
    pC.y        = 1;
    pC.z        = 1;

    % place new model states and priors in generative model
    %--------------------------------------------------------------------------
    D.M(1).pE   = pE;
    D.M(1).pC   = diag(spm_vec(pC));

    disp(['Inverting model, ' model{cc}])

    % Inversion using generalised filtering
    %=========================================================================
    LAP         = spm_LAP(D);

    % use Bayesian model reduction to test different hypotheses
    %==========================================================================

    % apply precise shrinkage priors to of diagonal coupling elements
    %--------------------------------------------------------------------------
    clear PC
    PC{1}     = pC; PC{1}.y = 0;
    PC{2}     = pC; PC{2}.z = 0;

    %  evaluate the evidence for these new models or prior constraints
    %--------------------------------------------------------------------------
    qE          = LAP.qP.P{1};
    qC          = LAP.qP.C;
    pE2         = LAP.M(1).pE;
    pC2         = LAP.M(1).pC;
    clear F
    for m       = 1:numel(PC)
        rC      = diag(spm_vec(PC{m}));
        F(m)    = spm_log_evidence(qE,qC,pE2,pC2,pE2,rC);
    end

    % report marginal log likelihood or evidence
    %--------------------------------------------------------------------------
    F           = F - min(F);

    k = k + 1;
    subplot(4,4,k)
    bar(F,'c')
    title([model{cc} ', log ev.'])
    axis square,
    box off

    k = k + 1;
    subplot(4,4,k)
    bar(spm_softmax(F(:)),'c')
    title([model{cc} ', prob.'])
    axis square,
    box off,
    drawnow

end

