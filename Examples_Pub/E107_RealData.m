% ------------------------------------------------------------------------- 
%                               CTP-Analysis 
% 
% 
% 
% 
% 
%                                          (c)Constantin Heck, 07-Aug-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


clc;
close all;
clear;



%% settings
%{
%Parameters for 128x128x80
dataset  = 'D2';
m        = [128,128,80];
k        = 24;

maskMode = 'head';  %mask where to do the deconvolution
tRes     = 5;          %time resolution  in seconds after interpolation
sd       = .5;           %sd for prior smoothing
fsize    = [3,3,3];   
thres    = .04;         %global threshold for svd
memorySave = false;

%}



%parameters for 512x512x320
% {
dataset  = 'D2';
m        = [512,512,320];
k        = 24;

maskMode = 'awesome';  %mask where to do the deconvolution
tRes     = 5;          %time resolution after interpolation
sd       = 1.5;          %sd for prior smoothing
fsize    = [5,5,5];   
thres    = .04;        %global threshold for svd
memorySave = true;
%}

fpath = '/Volumes/Macintosh_home/check/Documents/data/CTP-Matlab/';
fname = sprintf('%s_%ix%ix%ix%iint16.nii',dataset,m(1),m(2),m(3),k);





%% load and prepare data (caution, might take a while)
% {

fprintf('Loading data...');tic;

%setup main variables
n  = prod(m);
mk = [m,k];

%generate filename

%load data
[data,omega,m] = loadnii_lars([fpath,fname]);
data = double(data.img);

% load timeline
fname = sprintf('%s_timeline.txt',dataset);
fid   = fopen([fpath,fname],'r');
timeline = fscanf(fid,'%f');
fclose(fid);

% load mask for aif
fname = sprintf('%s_maskAif_%ix%ix%i.mat',dataset,m(1),m(2),m(3));
load([fpath,fname]);
maskAif = int16(maskAif);

fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);





data(:,:,:,(end-2:end)) = [];
timeline(end-2:end) = [];
k = numel(timeline);
mk = [m,k];


%% setup global correction factor for CBF and CBV. See
%1) doi:10.1097/00004647-199907000-00001
%2) Technical Aspects of Perfusion-Weighted Imaging (Wu 2005)
%proportionality constants for conversion are supposed to be equal.
%
% Then C = sfac*CBF*(R(x)aif); 
% and CBV = 1/sfac*intintegral(C)/integral(aif)


HctSmall = .25; %small vessel hematocrit
HctLarge = .45; %large vessel hematocrit
rho      = 1.04; %brain density in g/mL

sfac     = rho*(1-HctSmall)/(1-HctLarge);





%% smooth input data

if sd>0
    fprintf('Smoothing data...');tic;
    D = data;
    for i = 1:k
        D(:,:,:,i) = smooth3(data(:,:,:,i),'gaussian',fsize,sd);
        fprintf('.%1.0f%%.',i/k*100);
    end
    fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);
elseif sd == 0;
    D = data;
else
    error('sd does not make sense');
end

if memorySave
    clear('data');
end



%% generate mask

fprintf('Setting up mask...');tic;
switch maskMode
    case 'head'
        D1   = squeeze(D(:,:,:,1));
        mask = (D1>-150) & (D1<1000);
        
        if memorySave
            clear('D1');
        end
    case 'smallBox'
        
        %setup image center and small width/ big width
        c = round(m/2);
        w = round(m*2.5/100);
        idxi = c(1)-w(1) + (1:2*w(1));
        idxj = c(2)-w(2) + (1:2*w(2));
        idxk = c(3)-w(3) + (1:2*w(3));
        
        
        mask = zeros(m);
        mask(idxi,idxj,idxk) = 1;

    case 'bigBox'
        
        %setup image center and small width/ big width
        c = round(m/2);
        w = round(m*10/100);
        idxi = c(1)-w(1) + (1:2*w(1));
        idxj = c(2)-w(2) + (1:2*w(2));
        idxk = c(3)-w(3) + (1:2*w(3));
        
        
        mask = zeros(m);
        mask(idxi,idxj,idxk) = 1;
        
    case 'awesome'
        
        fname = sprintf('%s_%ix%ix%iint16-segm.nii',dataset,m(1),m(2),m(3));
        tmp = loadnii_lars([fpath,fname]);
        mask = tmp.img;
        clear('tmp');
        
    otherwise
        error('maskMode unknown');
end
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);




%% show something of the input data
%{
scrollView(data,omega,m,3,'name','Original Data','fig',1);
scrollView(D,omega,m,3,'name','Smoothed data','fig',2);
scrollView(D(:,:,:,1),omega,m,3,'mask',maskAif,'name','Smoothed data with AIF','fig',3);
scrollView(D(:,:,:,1),omega,m,3,'mask',mask,'name','Smoothed data with mask','fig',4);
return;
%}

%% get uptake curves and convert to concentrations
%reshape D to image-first
D       = reshape(D,n,k);


%get uptake curves
fprintf('Getting uptake curves...');tic;
idxMask = (mask(:)~=0);
tC      = D(idxMask,:);
ncurve  = size(tC,1);
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);

%convert to concentrations
fprintf('Calculating concentrations...');tic;
C      = bsxfun(@minus,tC,tC(:,1));
C(C<0) = 0;
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);




%% get aif
fprintf('Calculating aif curve...');tic;
idxAif     = (maskAif(:)~=0);
tAif       = D(idxAif,:);
tAif       = mean(tAif,1);
aif        = tAif-tAif(1);
aif(aif<0) = 0;
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);

%reshape D to standard sizes
D = reshape(D,mk);

if memorySave
    clear('tC','tAif','D');
end


% figure(1);clf; plot(timeline,aif);

%% interpolation to regulary spaced timeline




fprintf('Starting interpolation....');tic;

%new timeline
timelineNew = (0:tRes:timeline(end));
kNew        = numel(timelineNew);
T           = perfusion1c.getLinearInterMatrix1D(timeline,timelineNew);

%setup interpolated curves
aifNew = T*aif(:);
CNew   = C*T';

fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);


if memorySave
    clear('C','aif');
end


%% do deconvolution analysis

[F,phi,Irec,Crec,para] = perfusion1c.fastPerfusionAnalysis(CNew,aifNew,timelineNew,thres);


%% calculate CBF and CBV


CBF = zeros(m);
CBF(idxMask) = F*100*60/sfac;


CBV = zeros(m);
CBV(idxMask) = phi*100/sfac;


if memorySave
    clear('F','phi');
end


%% plot nDisp curves from Irec, Crec, C

nDisp = 1; %number of curves to plot

%setup indices and a new, longer timeline for display
idxRand      = randperm(ncurve,nDisp);
timelineNewL = linspace(0,2*timeline(end),2*kNew);

figure(1);clf;

subplot(1,3,1);
plot(timelineNewL,Irec(idxRand,:)');
title('IRec');

subplot(1,3,2);
plot(timelineNewL,Crec(idxRand,:)','r',timelineNew,CNew(idxRand,:)','b');
title('CRec and C');


subplot(1,3,3);
plot(timelineNew,aifNew);
title('AIF');



figure(4);clf;
set(4,'name',sprintf('Mean Curves. CBF(Cmean)=%1.3f, mean(CBF)=%1.3f',para.Fmean*100*60/sfac,para.meanF*100*60/sfac));
subplot(1,3,1);
plot(timelineNewL,para.Imean);
title('IRec');

subplot(1,3,2);
plot(timelineNewL,para.CmeanR,'r',timelineNewL,para.Cmean,'b');
title('CRec and C');


subplot(1,3,3);
plot(timelineNew,aifNew);
title('AIF');



%% show CBF and CBV

% 
scrollView(CBF,omega,m,3,'name','CBF in ml/min/100ml','fig',2);
scrollView(CBV,omega,m,3,'name','CBV in percent','fig',3);
%}

return;
%% save everything as matlab file
%{
fname = sprintf('%s_%ix%ix%i_mask-%s_sd-%1.1f_thres-%1.2f.mat',dataset,m(1),m(2),m(3),maskMode,sd,thres);
sp
save([fpath,fname],'-v7.3');
%}

%% save CBF for mevislab

fname = sprintf('%s_PK_%ix%ix%i_mask-%s_sd-%1.1f_thres-%1.2f.nii',dataset,m(1),m(2),m(3),maskMode,sd,thres);

PKParameters = cat(4,CBF,CBV);
savenii_lars(PKParameters,omega,m,[fpath,fname],'First: CBF in ml/min/100ml, Second: CBV in percent');


return;
%% graphics paper


%plot Cmean vs CMean R
figure(5);clf;
plot(timelineNewL,para.Cmean,'b',timelineNewL,para.CmeanR,'r');
xlabel('time [s]');
ylabel('rel. concentration');
matlab2tikz('real_meanC.tikz','width','\fwd');


%% plot AIF
figure(5);clf;
plot(timelineNew,aifNew);
xlabel('time [s]');
ylabel('rel. concentration');
matlab2tikz('real_AIF.tikz','width','\fwd');


%% plot CBF
scrollView(CBF,omega,m,3,'cmap','parula(512)','scale','slices');
caxis(ca);

%%

Si = squeeze(CBF(:,:,160));

figure(5);clf;
colormap('parula(512)');
imagesc(Si);
axis image;
axis off;
% caxis([0,500]);
caxis([0,150]);
colorbar('south','color','white')
export_fig('real_axial160.pdf','-transparent',5);
