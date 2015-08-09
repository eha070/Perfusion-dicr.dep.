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
i = 7;
D = squeeze(data(:,:,:,i));
scrollView(D,omega,m,3,'scale','slices');



%% show slice
i = 151;
close all;
clc;

Di = squeeze(D(:,:,i));
mi = m(1:2);

figure(1);clf;
viewImage2Dsc(Di(:),omegaM(mi),mi,'colormap','gray(512)');
setGrayValueWindowFcnt

%% setup 2D mask

%z slices to create the mask
z = [151,151,152,152];



close all;
mi = m(1:2);

W = false(m);
for i = z;

    Di = squeeze(D(:,:,i));
    

    figure(1);clf; colormap gray(512);
    ca = [-111.288 325.34];
    imagesc(Di);
    axis image;
    caxis(ca);


    h = imfreehand;
    W(:,:,i) = W(:,:,i) + h.createMask;

    close(1);

end    
   


scrollView(W,omega,m,3);