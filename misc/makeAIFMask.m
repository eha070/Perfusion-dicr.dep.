% ------------------------------------------------------------------------- 
%                               makeAIFMask 
% 
% 
% 
% 
% 
%                                          (c)Constantin Heck, 09-Aug-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 



%% show volume i
i = 10;
D = squeeze(data(:,:,:,i));
scrollView(D,omegaM(m),m,3,'scale','slices');



%% show slice
i = 39;
close all;
clc;

Di = squeeze(D(:,:,i));
mi = m(1:2);

figure(1);clf;
viewImage2Dsc(Di(:),omegaM(mi),mi,'colormap','gray(512)');
setGrayValueWindowFcnt

%% setup 2D mask

%set color axis
ca = [-44.7822 174.704];

%z slices to create the mask
z  = [38,39,40,41];



close all;
mi = m(1:2);

W = false(m);
for i = z;

    Di = squeeze(D(:,:,i));
    

    figure(1);clf; colormap gray(512);
    imagesc(Di);
    axis image;
    caxis(ca);


    h = imfreehand;
    W(:,:,i) = W(:,:,i) + h.createMask;

    close(1);

end    
   


scrollView(D,omegaM(m),m,3,'mask',W);


%% save data

fpath   = '~/Documents/data/CTP-Matlab/';
fname   = 'D2_maskAif_128x128x80.mat';

maskAif = W;

save([fpath,fname],'maskAif')