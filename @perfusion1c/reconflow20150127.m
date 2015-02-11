function [qmat,perfmat,perfmatn,phimat] = reconflow(prmin)
% [qmat,perfmat,perfmatn,phimat] = perfusion1c.reconflow(0,'','syntconv','full','classic');


% settings
[prm,Fmat] = settings;
prm = perfusion1c.mergestruct(prm,prmin);
perfusion1c.printstructscreen(prm);

% basename
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);

% indicator basename
basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);

    
if isequal(prm.reconflowdata,'synt')

    % load true field
    pathload = ['synt-createflowTPFA-' basenameflow '.mat'];
    msg = ['Loading ' pathload];
    disp(msg);
    D = load(pathload);
    perfmatnt = D.perfmatn;
    phimatt = D.phimat;

    pathload = [prm.reconflowdata '-createindicatorpde-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
    D = load(pathload);
    msg = ['Loading ' pathload];
    disp(msg);
    Cmat = D.Cmat;
    timeline = D.prm.timeline;
    aifval = D.aifval;
    h = prm.h;
    
elseif isequal(prm.reconflowdata,'syntconv')
    
    % load true field
    pathload = ['synt-createflowTPFA-' basenameflow '.mat'];
    msg = ['Loading ' pathload];
    disp(msg);
    D = load(pathload);
    perfmatnt = D.perfmatn;
    phimatt = D.phimat;

    pathload = [prm.reconflowdata '-createindicatorconv-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
    D = load(pathload);
    msg = ['Loading ' pathload];
    disp(msg);
    Cmat = D.Cmat;
    timeline = D.prm.timeline;
    aifval = D.aifval;
    h = prm.h;
    
elseif isequal(prm.reconflowdata,'brain')

    D = load('ep_perf_mri_20030919_shifted.mat');
    Cmat = D.im;
    h = D.h;
    h(3) = h(2);%??
    clear D;
    dim = size(Cmat);
    dim3 = dim(1:3);
    ntime = dim(4);
    dt = 1.44;
    timeline = (0:ntime-1)'*dt;
    a = mean(Cmat,4);
    roi.organ = a > 150;
    load ball3;se = getball(ball,3,3);
    roi.organ = imopen(roi.organ,se);
    roi.organ = bwkeep(roi.organ,1,6);
    roi.organ = imfill(roi.organ,'holes');
%         roi.organ = ones(dim3);
    D = load('ep_perf_mri_20030919_shifted-aif.mat');
    roi.aif = D.aifroi;
    TE = 0.6; % s? M? der v?re. 0.6ms er ekstremt kort
    r = 4.1;%s, taken from DCE-MRI sequence using Dotarem?
    a = zeros(dim);
    eps = 0.01;
    % convert to concentration time curves
    % Using equation 11.7a in Kennan, Jager, T2 and T2*-w DCE-MRI:
    % blood perfusion and Volume estimation using Bolus Tracking, page
    % 371 (part of the BMED360 literature)
    for j = 1 : ntime
        a(:,:,:,j) = log(Cmat(:,:,:,j)./(Cmat(:,:,:,1) + eps));
    end;
    a = -a/(TE*r);
    % no negative concentrations
    a(a < 0) = 0;
    Cmat = a;
    % v = im2vec4D(Cmat,roi.aif);
    phiin = 1;
end;

% compute porosity
phimat = perfusion1c.porosity(Cmat,aifval,timeline);

% NaN and inf in prm.reconflowdata???
Cmat(isinf(Cmat)) = 0;
Cmat(isnan(Cmat)) = 0;

dim = size(Cmat);
dim3 = dim(1:3);
if isequal(prm.reconflowroitype,'full')
    roi.organ = ones(dim3);
end;
% run perfusion prm.reconflowcode
if isequal(prm.reconflowcode,'minder')
    [qmat,perfmat,perfmatn] = perfusion1cminder(Cmat,timeline,roi.aif,phiin,roi.organ,h,prm.reconflowalpha,prm.reconflowregularization);
elseif isequal(prm.reconflowcode,'linear')
    [qmat,perfmat,perfmatn] = perfusion1c(Cmat,timeline,roi.aif,phiin,roi.organ,h,prm.reconflowalpha);
elseif isequal(prm.reconflowcode,'classic')
    [qmat,perfmat,perfmatn] = perfusion1c.reconflowclassic(Cmat,timeline,roi.organ,aifval,prm);
end;

savpaper = 1;
savdata = 1;

msg = ['Average reconstructed perfusion: ' num2str(mean(perfmatn(:)))];
disp(msg);
msg = ['True, average reconstructed perfusion: ' num2str(mean(perfmatnt(:)))];
disp(msg);

basenamerecon = perfusion1c.providenamerecon(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T,prm.reconflowdeconvmethod,prm.stepred);

if savpaper == 1
    figure(1)
    val = 100*(perfmatnt - perfmatn)./perfmatnt;
    imagesc(val);colormap(gray);axis off;colorbar;axis image;
    brighten(-0.5)
    pathsave = ['figs/' prm.reconflowdata  '-' mfilename  '-method-'  prm.reconflowcode  '-' basenamerecon '-perfn'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(gcf,pathsave,'-deps')

    try
        h2 = figure(2);
        imagesc(D.delay);colormap(gray);axis off;colorbar;axis image;
        pathsave = ['figs/' prm.reconflowdata  '-' mfilename  '-method-'  prm.reconflowcode  '-' basenamerecon  '-delay'];
        msg = ['Saving ' pathsave];
        disp(msg);
        print(gcf,pathsave,'-deps')
    catch
        close(h2);
    end;
    figure(3)
    val = 100*(phimat - phimatt)./phimatt;
    imagesc(val);colormap(gray);axis off;colorbar;axis image;
    brighten(-0.5)
    pathsave = ['figs/' prm.reconflowdata '-' mfilename '-method-'  prm.reconflowcode  '-' basenamerecon  '-porosity'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(gcf,pathsave,'-deps')

end;

if savdata == 1
    % save results
    pathsave = [prm.reconflowdata '-' mfilename  '-method-'  prm.reconflowcode  '-' basenamerecon  '.mat'];
    msg = ['Saving ' pathsave];
    disp(msg);
    save(pathsave,'roi','h','timeline','qmat','perfmat','perfmatn','phimat');
end;
    
