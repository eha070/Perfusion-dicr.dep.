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





%% settings

%{
%OK parameters for 128x128
maskMode = 'head';  %mask where to do the deconvolution
tRes     = .5;          %time resolution after interpolation
sd       = 1;           %sd for prior smoothing
fsize    = [3,3,3];   
thres    = .04;         %global threshold for svd
%}



%OK parameters for 512x512
maskMode = 'smallBox';      %mask where to do the deconvolution
tRes     = .5;          %time resolution after interpolation
sd       = 1;           %sd for prior smoothing
fsize    = [5,5,5];   
thres    = .04;         %global threshold for svd
memorySave = true;



%% load and prepare data (caution, might take a while)
% {

clearvars -except maskMode tRes sd fsize thres

fprintf('Loading data...');tic;

%setup main variables
dataset  = 'D2';
% m        = [128,128,80];
m        = [512,512,320];
k        = 24;

n  = prod(m);
mk = [m,k];

%generate filename
fname = sprintf('%s_%ix%ix%ix%iint16.nii',dataset,m(1),m(2),m(3),k);
fpath = '/Volumes/Macintosh_home/check/Documents/data/CTP-Matlab/';

%load data
[data,omega,m] = loadnii_lars([fpath,fname]);
data = double(data.img);

% load timeline
fname = sprintf('%s_timeline.txt',dataset);
fid   = fopen([fpath,fname],'r');
timeline = fscanf(fid,'%f');
fclose(fid);



% load aif
fname = sprintf('%s_maskAif_%ix%ix%i.mat',dataset,m(1),m(2),m(3));
load([fpath,fname]);

fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);




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


fprintf('Smoothing data...');tic;
D = data;
for i = 1:k
    D(:,:,:,i) = smooth3(data(:,:,:,i),'gaussian',fsize,sd);
    fprintf('.%1.0f%%.',i/k*100);
end
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);

if memorySave
    clear('data');
end



%}


%% generate mask


switch maskMode
    case 'head'
        D1   = squeeze(D(:,:,:,1));
        mask = (D1>-150);
        
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
        
    otherwise
        error('maskMode unknown');
end
        


% scrollView(D(:,:,:,1),omega,m,3,'mask',mask);



%% get uptake curves and convert to concentrations

%convert to concentrations
DConc   = reshape(D,n,k);
DConc   = bsxfun(@minus,DConc,DConc(:,1));
DConc(DConc<0) = 0;


idxMask = (mask(:)~=0);
C       = DConc(idxMask,:);
ncurve  = size(C,1);

% scrollView(DConc,omega,m,3);
% return;



%% get aif
idxAIF = (maskAif(:)~=0);
aif    = DConc(idxAIF,:);
aif    = mean(aif,1);

if memorySave
    clear('DConc');
end

%% interpolation to regulary spaced timeline

disp('Starting interpolation....');tic;

%new timeline
timelineNew = (0:tRes:timeline(end));
kNew        = numel(timelineNew);

%setup aiflc
aifNew = interp1(timeline,aif,timelineNew,'linear');

%setup C
CNew   = interp1(timeline,C',timelineNew,'linear');
CNew   = CNew';

fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);


if memorySave
    clear('C');
end


%% do deconvolution analysis

[F,phi,Irec,Crec] = perfusion1c.fastPerfusionAnalysis(CNew,aifNew,timelineNew,thres);


%% show CBF and CBV


CBF = zeros(m);
CBF(idxMask) = F*100*60/sfac;
scrollView(CBF,omega,m,3,'name','CBF in ml/min/100ml','fig',1);


CBV = zeros(m);
CBV(idxMask) = phi*100/sfac;
scrollView(CBV,omega,m,3,'name','CBV in percent','fig',2);


if memorySave
    clear('F','phi');
end


%% plot nDisp curves from Irec, Crec, C

nDisp = 1; %number of curves to plot


%setup indices and a new, longer timeline for display
idxRand      = randperm(ncurve,nDisp);
timelineNewL = linspace(0,2*timeline(end),2*kNew);

figure(2);clf;

subplot(1,3,1);
plot(timelineNewL,Irec(idxRand,:)');
title('IRec');

subplot(1,3,2);
plot(timelineNewL,Crec(idxRand,:)','r',timelineNew,CNew(idxRand,:)','b');
title('CRec and C');


subplot(1,3,3);
plot(timelineNew,aifNew);
title('AIF');


%% save results
%{
fpath = '/Volumes/Macintosh_home/check/Documents/data/CTP-Matlab/';
fname = sprintf('%s_%ix%ix%i_mask-%s_sd-%1.1f_thres-%1.2f.mat',dataset,m(1),m(2),m(3),maskMode,sd,thres);
sp
save([fpath,fname],'-v7.3');
%}

%% save CBF for mevislab

fpath = '/Volumes/Macintosh_home/check/Documents/data/CTP-Matlab/';
fname = sprintf('%s_PK_%ix%ix%i_mask-%s_sd-%1.1f_thres-%1.2f.nii',dataset,m(1),m(2),m(3),maskMode,sd,thres);

PKParameters = cat(4,CBF,CBV);
savenii_lars(PKParameters,omega,m,[fpath,fname],'First: CBF in ml/min/100ml, Second: CBV in percent');
