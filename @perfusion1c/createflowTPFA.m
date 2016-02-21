function [] = createflowTPFA


% settings
[prm,Fmat] = settings;
dim = prm.dim;


% Given diagonal permeability tensor K as ones everywhere
% value of K in the range of 1e-6
if isequal(prm.Kopt,'flat')
    % K = 5e-6*ones(dim);
    % n = 50*ones(dim);
    % r = 5*1e-3*ones(dim);
    % K = (1/24)*pi*n.*r.^4;
    % K = 5e-6*ones(dim);
    
    % Units in mm^2
    K = 5e-6*ones(dim);
    Kmat(:,:,:,1) = K;
    Kmat(:,:,:,2) = K;
    Kmat(:,:,:,3) = K;
    
elseif isequal(prm.Kopt,'rand')
    
    % % Or we give random values
    K = exp(5*smooth3(smooth3(randn(3,dim(1),dim(2),dim(3)))));
    K1(:,:,:,1) = squeeze(K(1,:,:,:));
    K1(:,:,:,2) = squeeze(K(2,:,:,:));
    K1(:,:,:,3) = squeeze(K(3,:,:,:));
    Kmat = K1;

%     K = 10*rand([dim,3]);
%     filt = fspecial('gaussian',3,1);
%     K(:,:,:,1) = imfilter(K(:,:,:,1),filt);
%     K(:,:,:,2) = imfilter(K(:,:,:,2),filt);
%     K(:,:,:,3) = imfilter(K(:,:,:,3),filt);
%     K = exp(K);
elseif isequal(prm.Kopt,'sin')

    [x,y] = ndgrid(1:dim(1),1:dim(2));
    a = 10*(sin(0.05*x-5) + cos(0.05*y-5));
    a = scale(a,0.1,30);
    Kmat(:,:,:,1) = a;
    Kmat(:,:,:,2) = a;
    Kmat(:,:,:,3) = a;
    meanKmat = mean(Kmat(:));
    Kmat = Kmat/meanKmat;
    Kmat = Kmat*0.1*1e-6;
%    show(Kmat(:,:,:,1),1)
%    pause
elseif isequal(prm.Kopt,'brain')
    
    % compute the distance from the source and sink
%     X = (linspace(1,35,dim(1)/2))';
%     X = repmat(X,1,dim(2));
%     d = X;
%     d = [d;flipud(d)];
    x = (1:35)';
    %n = (linspace(70,130,dim(1)/2))';
    %n = ((130-30)/35)*x + 30;
    n = 29+exp(0.13*x);
    n = repmat(n,1,dim(2));    
    n = [n;flipud(n)];
    
    % radius or arteries
%    r = (linspace(100,5,dim(1)/2))';
    r = 100-exp(0.13*x);
    r = r/1000;
    %r = r - 100;
    r = repmat(r,1,dim(2));    
    r = [r;flipud(r)];
    %r = r.^0.25;
    
%     
%     % diameter of arteries
%     r0 = 0.1;
%     % map of diameter
%     % r = 0.90.^(0.95*(d-1))*r0;
%     r = 0.95.^(0.95*(d-1))*r0;

    
    K = (1/24)*pi*n.*r.^4;
   % multiply by a gradient in x-direction 
   % (by EAH  24.11.14)
   % Structure of the gradient:
   %
   % 1 .x .x .x  0  
   % 1 .x .x .x .x
   % 1  1  1  1  1  
   % 1 .x .x .x .x
   % 1 .x .x .x  0
   
   y1=linspace(0.1,1,dim(2))';
   y2=linspace(0.1,1,dim(2)/2)';
   G=(1-bsxfun(@times,y1',y2));
   G=([flipud(G);G]);
   
   K = K.*G;   
   K(1,end) = K(1,end-1);
   K(end,end) = K(end,end-1);
  
   Kmat(:,:,:,1) = K;
   Kmat(:,:,:,2) = K;
   Kmat(:,:,:,3) = K;
   
end;

% make a porosity map
if isequal(prm.phiopt,'flat')
    phimat = 0.05*ones(dim);
elseif isequal(prm.phiopt,'sin')
    [x,y] = ndgrid(1:dim(1),1:dim(2));
    phimat = sin(0.1*x) + cos(0.1*y);
    phimat = scale(phimat,0,1);
elseif isequal(prm.phiopt,'brain')
    phimat = (pi.*r.*r.*n)/(prmin.h(1)*prmin.h(2));
end;

% find the pressure field and flux
% [pmat,Fmat] = perfusion1c.syntTPFAsimp(Kmat,Fmat,prm);
[pmat,qmat] = perfusion1c.syntTPFA(Kmat,Fmat,prm);
% Fmat{1} = Fmat{1}/100;
% Fmat{2} = Fmat{2}/100;

% scale the pressure to minimum zero pressure. Does not change the problem
pmat = pmat - min(pmat(:));
msg = ['Maximum (scaled) pressure (kPa): ' num2str(max(pmat(:)))];
disp(msg);

% show(Fmat{1},1);colorbar
% show(Fmat{2},2);colorbar
qin = sum(abs(Fmat(:)))/2;
msg = ['Flux in: ' num2str(qin)];
disp(msg);

mid = round(dim(1)/2);
a = qmat{1}(mid,:);
qslice = sum(a);
msg = ['Flux x slice: ' num2str(qslice)];
disp(msg);

mid = round(dim(2)/2);
a = qmat{2}(:,mid);
qslice = sum(a);
msg = ['Flux y slice: ' num2str(qslice)];
disp(msg);

% Fmat{1}
% Fmat{2}
% pause

makelenmat = 1;
if makelenmat == 1
    msg = ['Making and saving lenmat'];
    disp(msg);
    lenmat = perfusion1c.arclength(qmat,Fmat,prm.h);
    save lenmat.mat lenmat prm;
end;

% find the streamlines
D = load('lenmat.mat');
if ~isequal(D.prm.dim,prm.dim) || ~isequal(D.prm.physdim,prm.physdim)
    msg = ['Making and saving lenmat'];
    disp(msg);
    lenmat = perfusion1c.arclength(qmat,Fmat,prm.h);
    save lenmat.mat lenmat prm;
else
    msg = ['Loading lenmat'];
    disp(msg);
    lenmat = D.lenmat;
end;
% convert flux to perfusion
[perfmat,perfmatn] = perfusion1c.flux2perf(qmat,lenmat,Fmat,prm.h);

% test the inflow and the estimated perfusion values to see that they agree
% with each other
totperf = sum(perfmat(:))*prm.voxelvol;
msg = ['Total (unnormalized) perfusion (m^3/s): ' num2str(totperf)];
disp(msg);
msg = ['Mean (normalized) perfusion (ml/min/100ml): ' num2str(mean(perfmatn(:)))];
disp(msg);

% show(pmat,3);colorbar
% show(perfmat,4);colorbar
% show(perfmatn,5);colorbar


savpaper = 0;
savdata = 1;
basename = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);
if savpaper
    h = figure(1);
    imagesc(Kmat(:,:,:,1));colormap(gray);axis off;colorbar;axis image;
    pathsave = ['figs' '/synt-' mfilename '-' basename '-K.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

    h = figure(2);
    imagesc(pmat);axis off;colormap(gray);colorbar;axis image;
    pathsave = ['figs' '/synt-' mfilename '-' basename '-p.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

    h = figure(3);
    imagesc(phimat);axis off;colorbar;axis image;colormap(gray);
    pathsave = ['figs' '/synt-' mfilename '-' basename '-phi.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

    h = figure(4);
    imagesc(qmat{1});colormap(gray);axis off;colorbar;axis image;
    brighten(0.7);
    pathsave = ['figs' '/synt-' mfilename '-' basename '-fluxx.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

    h = figure(5);
    imagesc(qmat{2});colormap(gray);axis off;colorbar;axis image;
    brighten(0.7);
    pathsave = ['figs' '/synt-' mfilename '-' basename '-fluxy.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

    h = figure(6);
    imagesc(perfmat);colormap(gray);axis off;colorbar;axis image;
    brighten(0.7);
    pathsave = ['figs' '/synt-' mfilename '-' basename '-perf.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

    h = figure(7);
    imagesc(perfmatn,[0,300]);colormap(gray);axis off;colorbar;axis image;
    brighten(0.7);
    pathsave = ['figs' '/synt-' mfilename '-' basename '-perfn.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')
    
    h = figure(8);
    imagesc(lenmat);colormap(gray);axis off;colorbar;axis image;
    pathsave = ['figs' '/synt-' mfilename '-' basename '-lenmat.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(h,pathsave,'-deps')

end;

if savdata
    pathsave = ['results/synt-' mfilename '-' basename '.mat'];
    msg = ['Saving ' pathsave];
    disp(msg);
    save(pathsave,'pmat','qmat','perfmat','perfmatn','phimat','Fmat','prm','Kmat','-v7.3')
end;

