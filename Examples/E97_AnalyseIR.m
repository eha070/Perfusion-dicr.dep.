% ------------------------------------------------------------------------- 
%                               Description 
% 
% 
% 
% 
% 
%                                          (c)Constantin Heck, 09-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 



clc;
close all;
    
%which results to show?
showFlowMaps       = 0;
showMultipleCurves = 0; %remember to setup indices idxiD and idxjD
showSingleCurve    = 0;
playground         = 1;
renewCalculation   = 0;


%setup area where to run the deconvolution
idxi  = (1:64);
idxj  = (1:64);


%setup oscillation index OI
OI = .5; %play with it!



%% load data

% settings
prm               = settings;
prm.aiftype       = 'delta';
basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
basenameFlow      = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim');
foldername        = './results/';

%setup paths
pathloadFlow = ['synt-createflowTPFA-' basenameFlow '.mat'];
pathload = [foldername,'synt-createindicatorpde-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
% pathload = [foldername,'synt-createindicatorpde-' basenameindicator
% '.mat']; %full dataset

%load concentration curves
msg = ['Loading ' pathload];
disp(msg);
D   = load(pathload);

%setup C etc
Cmat     = squeeze(D.Cmat);
AIF      = D.aifval;
timeline = D.prm.timeline;


%setup CBF
msg = ['Loading ' pathloadFlow];
disp(msg);
E = load(pathloadFlow);


%perfmat
qmat     = E.qmat;
CBF      = E.perfmat; 
CBF(1,1) = 2*CBF(1,1);
CBF(end,end) = 2*CBF(end,end);


%averaging
% qmat = E.qmat;
% qcc = perfusion1c.convertFlowStagToCC(qmat);
% CBF = sqrt(qcc{1}.^2 + qcc{2}.^2)/2;




%% prepare data
if renewCalculation

    %prepare sizes
    m       = [size(Cmat,1),size(Cmat,2)];
    n       = prod(m); %number of datapoints
    k       = size(Cmat,3); %number of timepoints
    mk      = [m,k]; %...


    %setup area where to run the deconvolution

    %downsmaple the data
    numi        = numel(idxi)*numel(idxj);

    %initialize variables to store results
    CBFest = zeros(m);
    Iest   = zeros([m,k]);
    Cest   = zeros([m,k]);

    %get deconvolution matrix
    deltaT  = timeline(2)-timeline(1);
    A       = perfusion1c.getLinearConvolutionMatrix(AIF,deltaT);
    fprintf('Starting SVD...');
    tic; [U,S,V] = svd(A);
    fprintf('...done. Elapsed time: %1.3fs\n',toc);

    %do the job voxelwise
    h = waitbar(0);
    for i = idxi;
        for j = idxj;


            %do the deconvolution
            Cij = squeeze(Cmat(i,j,:));
            [F,Irec,Crec] = perfusion1c.linearDeconvolution(Cij,timeline,OI,U,S,V);
            Irec          = A\Cij(:);
            F             = max(Irec);
            Crec          = A*Irec;

            %store results
            CBFest(i,j) = F;
            Iest(i,j,:) = Irec;
            Cest(i,j,:) = Crec;

            %fill the waitbar
            perc = (numel(idxj)*(i-min(idxi)) + j-min(idxj)+1 )/numi;
            waitbar(perc,h,sprintf('Calculating...%2.2f%%...',perc*100));


        end
    end

    close(h);
end


return;
%% results for a single voxel

if showSingleCurve
    
    %position of single curve
    pos = [35,35];

    s = @(v) squeeze(v(:));

    figure(3);clf;
    set(3,'name','Comparison on single voxel');

    subplot(1,2,1);
    plot(timeline,s(Cmat(pos(1),pos(2),:)),timeline,s(Cest(pos(1),pos(2),:)),'LineWidth',3);
    title('Ctrue (blue) and reconstructed C (red) at position pos')
    
    subplot(1,2,2);
    plot(s(Iest(pos(1),pos(2),:)),'LineWidth',3);
    title('estimated I at position pos')
end








%% PLAYGROUND!



if playground
    
    
    % setup parameters
    k = numel(timeline);
    pathnameGamma = './results/synt-createindicatorpde-phi-flat-K-flat-dim-64-aif-gamma-T-90-red-100.mat';
    E        = load(pathnameGamma);    
    
    %different values for indices
%     idx = [1,1];
%     idx = [35,50]; p0 = [1.8e-4;3.2;1];
%     idx = [60,40]; p0 = [5e-5;6;1];
    idx = randi([1,64],1,2)

    % get CDelta, IDelta
    CDelta = squeeze(Cmat(idx(1),idx(2),:));
    IDelta = squeeze(Cmat(idx(1),idx(2),:));


%     %analyse pOpt
%     p0   = [1.8e-4;2.2;1];
%     pOpt = FitGaussian(timeline,CDelta,p0,'showIni',false);
%     pOpt
%     pOpt(2)^2
    
    
    % get CGamma
    %do the deconvolution: Get deconvolution matrix 
    deltaT  = timeline(2)-timeline(1);
    A       = perfusion1c.getLinearConvolutionMatrix(E.aifval,deltaT);
    [U,S,V] = svd(A);

    %do the job voxelwise
    CGamma = squeeze(E.Cmat(idx(1),idx(2),:));
    [F,IGamma,Crec] = perfusion1c.linearDeconvolution(CGamma,timeline,OI,U,S,V);


    %get scaling factor
    sfac = max(IGamma)/max(IDelta)
    
    %show results
    figure(1);clf;
    plot(timeline,IGamma,timeline,IDelta)
    legend('IGamma','IDelta');
   
    norm(IGamma-IDelta)/norm(IGamma);
    
    % compare the integrals;
    Av = spdiags(ones(k,2),[0,1],k-1,k);
    e  = timeline(2:end)-timeline(1:end-1);
    e  = e*Av;
    e  = e(:);
    
    

end