% ------------------------------------------------------------------------- 
%                    E02_CircularDeconvolutionConvolution
% 
% Circular deconvolution on one of the synthetic flow models.
% 
% 
%                                          (c)Constantin Heck, 19-Nov-2014 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clc;
close all;



%which results to show?
showFlowMaps       = 1;
showMultipleCurves = 1; %remember to setup indices idxiD and idxjD
showSingleCurve    = 0;
writeImage         = 0;
saveData           = 0;
saveSingleCurve    = 0;


%setup area where to run the deconvolution
idxi  = (32:40);
idxj  = (32:40);


%setup oscillation index OI
OI = .005; %probably needs to be tuned



%% setup area where to run the deconvolution

%downsmaple the data
numi        = numel(idxi)*numel(idxj);


%initialize variables to store results
CBFest = zeros(m);
Iest   = zeros([m,2*k]);
Cest   = zeros([m,2*k]);

%get deconvolution matrix
deltaT  = timeline(2)-timeline(1);
A       = perfusion1c.getCircularConvolutionMatrix(aifval,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);


%do the job voxelwise
fh = waitbar(0);
for i = idxi;
    for j = idxj;
                    
       
        %do the deconvolution
        Cij = squeeze(Cmat(i,j,:));
        [F,Irec,Crec] = perfusion1c.circularDeconvolution(Cij,timeline,OI,U,S,V);
        
        %store results
        CBFest(i,j) = F;
        Iest(i,j,:) = Irec;
        Cest(i,j,:) = Crec;
        
        %fill the waitbar
        perc = (numel(idxj)*(i-min(idxi)) + j-min(idxj)+1 )/numi;
        waitbar(perc,fh,sprintf('Calculating...%2.2f%%...',perc*100));

        
    end
end

close(fh);




%% show results


if showFlowMaps
    
    
    %setup "segmented" Flow CBFseg
    seg = zeros(m);
    seg(idxi,idxj) = 1;
    CBFseg = CBF.*seg;

    %setup scaling
    tmp  = 100*60*[reshape(CBFseg(CBFseg~=0),[],1);reshape(CBFest(CBFest~=0),[],1)];
    cmin = min(tmp);
    cmax = max(tmp);

    %setup 
    REMap = abs(CBFest-CBFseg)./CBFseg;

    %show flow
    figure(1);clf;
    colormap jet(512);
    set(gcf,'name','Results for Convolution Model')


    subplot(1,3,1);
    imagesc(CBFseg*100*60);
    caxis([cmin,cmax]);
    axis image;
    title('True flow (ml/min/100ml)')

    subplot(1,3,2);
    imagesc(CBFest*100*60);
    caxis([cmin,cmax]);
    axis image;
    title('Estimated flow (ml/min/100ml)')

    %show division
    subplot(1,3,3);
    imagesc(REMap);
    caxis([0,1]);
    axis image;
    title('Relative Error in Flows |Fest-Ftrue|/Ftrue')

end


%% show some curves


if showMultipleCurves

    %setup curves to display
    idxiD       = (32:5:40);
    idxjD       = (32:5:40);    
    numi        = numel(idxiD)*numel(idxjD);
    timelineExt = linspace(0,2*timeline(end),2*k);




    Ctr = Cmat(idxiD,idxjD,:);
    Ctr = reshape(Ctr,numi,[]);
    Ctr = Ctr';

    Crec = Cest(idxiD,idxjD,:);
    Crec = reshape(Crec,numi,[]);
    Crec = Crec';

    Irec = Iest(idxiD,idxjD,:);
    Irec = reshape(Irec,numi,[]);
    Irec = Irec';


    REF = (CBFest(idxiD,idxjD)-CBF(idxiD,idxjD))./CBF(idxiD,idxjD);
    REF = mean(REF(:));


    figure(2);clf;

    subplot(1,2,1);
    plot(timelineExt,Irec,'-r','lineWidth',1);
    title(sprintf('Reconstructed I, Average RE in Flow=%1.4f',REF))

    subplot(1,2,2);
    plot(timeline,Ctr,'-b',timelineExt,Crec,'-r','lineWidth',1);
    title(sprintf('True C'))
    
end

%% results for a single voxel

if showSingleCurve
    
    %position of single curve
    pos = [32,32];

    s = @(v) squeeze(v(:));

    figure(3);clf;
    set(3,'name','Comparison on single voxel');

    subplot(1,2,1);
    plot(timelineLow,s(Clow(pos(1),pos(2),:)),timelineHigh,s(Cest(pos(1),pos(2),:)));
    title('Ctrue (blue) and reconstructed C (red) at position pos')
    
    subplot(1,2,2);
    plot(timelineHigh,s(Iest(pos(1),pos(2),:)));
    title('estimated I at position pos')
end









return;
%% playground

% idx = [32,32];
idx = randi([1,64],1,2);

% i1 = IestDeltaPhi(idx(1),idx(2),1,1:k);
i1 = IestGamma(idx(1),idx(2),1,:);
i2 = IestDeltaB(idx(1),idx(2),1,:);
% i2 = IestUHigh(idx(1),idx(2),1,:);

i1  = i1(:);
i2  = i2(:);

%scaling
s   = max(i2)/max(i1)
i1s = s*i1;


% plot(timeline',i1,timeline',i2,timeline',i1s)
% title(sprintf('(i,j)=(%i,%i), RE=%1.4f, s=%1.4e',idx(1),idx(2),norm(i1s-i2)/norm(i2),s));
% legend('i1','i2','i1s')









