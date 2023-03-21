% Generate Fig. 9
clear;clc;close all;

% directory where model inversions will be saved
outdir = '';

% directory where temporary averages will be saved
middir = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory where HCP data is saved
% datadir = '';
% subs   = {'rest_hcp160123_ts'...
%     'rest_hcp154734_ts'...
%     'rest_hcp155637_ts'...
%     'rest_hcp156637_ts'...
%     'rest_hcp159340_ts'...
%     'rest_hcp161731_ts'...
%     'rest_hcp162733_ts'...
%     'rest_hcp163129_ts'...
%     'rest_hcp176542_ts'...
%     'rest_hcp178950_ts'...
%     'rest_hcp188347_ts'...
%     'rest_hcp189450_ts'...
%     'rest_hcp190031_ts'...
%     'rest_hcp192540_ts'...
%     'rest_hcp196750_ts'...
%     'rest_hcp198451_ts'...
%     'rest_hcp199655_ts'...
%     'rest_hcp201111_ts'...
%     'rest_hcp208226_ts'...
%     'rest_hcp211417_ts'...
%     'rest_hcp211720_ts'...
%     'rest_hcp212321_ts'...
%     'rest_hcp214423_ts'...
%     'rest_hcp221319_ts'...
%     'rest_hcp239944_ts'...
%     'rest_hcp245333_ts'...
%     'rest_hcp280739_ts'...
%     'rest_hcp298051_ts'...
%     'rest_hcp366446_ts'...
%     'rest_hcp397760_ts'...
%     'rest_hcp414229_ts'...
%     'rest_hcp499566_ts'...
%     'rest_hcp654754_ts'...
%     'rest_hcp672756_ts'...
%     'rest_hcp751348_ts'...
%     'rest_hcp756055_ts'...
%     'rest_hcp792564_ts'...
%     'rest_hcp856766_ts'...
%     'rest_hcp857263_ts'...
%     'rest_hcp899885_ts'};
%
% % save HCP data
% z      = [];
% for ii = 1:numel(subs)
%     disp(['saving HCP data, ' num2str(round(ii*100/numel(subs))) '% complete'])
%     clear zt
%     zt = load([datadir subs{ii}]);
%     z  = cat(3,z,zt);
% end
% save('HCPdat.mat','z')
%
% % save HCP data correlated with functional gradients
% for jj = [1 3]
% directory where gradients are saved
%     gr = importdata(['/Margulies/' num2str(jj)]);
%
%     z      = [];
%     for ii = 1:numel(subs)
%         disp(['saving gradient correlated data, ' num2str(round(ii*100/numel(subs))) '% complete'])
%         clear zt
%         zt = load([datadir subs{ii}]);
%         zt = corr(zt',gr);
%         z  = cat(3,z,zt);
%     end
%
%     save(['HCPdat_grads_' num2str(jj) '.mat'],'z')
%
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform model inversions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for jj = [1 3]
%
%     clear z
%     load(['HCPdat_grads_' num2str(jj) '.mat'],'z')
%
%     c = -1;
%     for ii = 1:size(z,3)
%         disp(['inverting model, gradient number' num2str(jj) ', ' num2str(round(ii*100/size(z,3))) '% complete'])
%         c = c*-1;
%         if c == 1
%             dat = z(1:600,:,ii)';
%         else
%             dat = z(601:1200,:,ii)';
%         end
%         LAP = modelinv(dat);
%         save([outdir 'subj_' num2str(ii) '_gr' num2str(jj) '.mat'],'LAP')
%     end
%
%     c = 1;
%     for ii = 1:size(z,3)
%         disp(['inverting model, alternative split, gradient number' num2str(jj) ', ' num2str(round(ii*100/size(z,3))) '% complete'])
%         c = c*-1;
%         if c == 1
%             dat = z(1:600,:,ii)';
%         else
%             dat = z(601:1200,:,ii)';
%         end
%         LAP = modelinv(dat);
%         save([outdir 'subj_' num2str(ii) '_gr' num2str(jj) '_alt.mat'],'LAP')
%     end
%
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate figure 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figselect = [12 1; 24 1; 3 3; 10 3];

figure
kk = 1;
for aa = 1:size(figselect,1)

    clear LAP
    load([outdir 'subj_' num2str(figselect(aa,1)) '_gr' num2str(figselect(aa,2)) '_alt.mat'],'LAP')
    subplot(4,4,kk)
    plot(LAP.Y','k')

    kk = kk + 1;
    subplot(4,4,kk)
    if aa == 1
        plot(LAP.Y',ones(size(LAP.Y')),'k.')
    else if aa == 2
            plot(LAP.Y(1:end-1)',diff(LAP.Y'),'k.')
    else if aa == 3
            plot3(LAP.Y(1,:),LAP.Y(2,:),LAP.Y(3,:),'k.')
    end
    end
    end

    clear pC
    pC.A      = ones(figselect(aa,2));
    pC.y      = 1;
    pC.z      = 1;

    clear PC
    PC{1}     = pC; PC{1}.y = 0;
    PC{2}     = pC; PC{2}.z = 0;

    qE               = LAP.qP.P{1,1};
    qC               = LAP.qP.C;
    pE               = LAP.M(1).pE;
    pC               = LAP.M(1).pC;
    clear F
    for m            = 1:numel(PC)
        rC           = diag(spm_vec(PC{m}));
        F(m,1)       = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    end
    F                = F - min(F);

    kk = kk + 1;
    subplot(4,4,kk)
    bar(F,'c')
    title('Log evid')
    axis tight, box off

    kk = kk + 1;
    subplot(4,4,kk)
    bar(spm_softmax(F(:)),'c')
    title('Prob')
    axis tight, box off

    kk = kk + 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform model averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
kk = 1;
for jj = [1 3]

    mkdir(middir)
    rmdir(middir,'s')
    mkdir(middir)

    clear LAPs
    c = 0;
    for ii = 1:40
        clear LAP

        toload = [outdir 'subj_' num2str(ii) '_gr' num2str(jj) '.mat'];

        if exist(toload)

            clear LAP
            load(toload,'LAP')

            if isfield(LAP,'qP')
                c = c + 1;
                LAPs(c) = LAP;
            end

        end
    end

    c = 0;
    for ii = 1:numel(LAPs)
        if ~isempty(LAPs(ii).qP)
            c        = c + 1;
            clear DCM
            DCM.Ep   = LAPs(ii).qP.P{1,1};
            DCM.Cp   = LAPs(ii).qP.C;
            DCM.F    = LAPs(ii).F(end);
            DCM.M.pC = LAPs(ii).M(1).pC;
            DCM.M.pE = LAPs(ii).M(1).pE;
            save([middir num2str(c) '.mat'],'DCM')
        end
    end

    clear names
    for ii           = 1:c
        names{ii}    = [middir num2str(ii) '.mat'];
    end

    DCMav            = spm_dcm_average(names,'kav');

    clear pC
    pC.A      = ones(jj);
    pC.y      = 1;
    pC.z      = 1;

    clear PC
    PC{1}     = pC; PC{1}.y = 0;
    PC{2}     = pC; PC{2}.z = 0;

    qE               = DCMav.Ep;
    qC               = DCMav.Cp;
    pE               = DCMav.M.pE;
    pC               = DCMav.M.pC;
    clear F
    for m            = 1:numel(PC)
        rC           = diag(spm_vec(PC{m}));
        F(m,1)       = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    end
    F                = F - min(F);

    subplot(2,2,kk)
    bar(F,'c')
    title('Log evid')
    axis tight, box off

    kk = kk + 1;

    subplot(2,2,kk)
    bar(spm_softmax(F(:)),'c')
    title('Prob')
    axis tight, box off

    rmdir(middir,'s')

    kk = kk + 1;

end

function LAP = modelinv(dat)

n             = size(dat,1);

x.y           = ones(n,1);
x.z           = ones(n,1);

P.A           = -eye(n)/4;
P.y           = 0;
P.z           = 0;

pC            = spm_unvec(ones(spm_length(P),1),P);

g             = @(x,v,P)  P.y*x.y + P.z*x.z;

f             = @(x,v,P) [x.z; P.A*(P.y*x.y + P.z*x.z)];

D.M(1).x        = x;
D.M(1).pE       = P;
D.M(1).pC       = diag(spm_vec(pC));
D.M(1).f        = f;
D.M(1).g        = g;
D.M(1).V        = exp(4);
D.M(1).W        = exp(4);

D.M(2).v        = 0;
D.M(2).V        = exp(4);

D.Y             = dat;

LAP             = spm_LAP(D);

end
