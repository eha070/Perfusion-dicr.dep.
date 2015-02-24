function [] = createindicatorpde(prmin)


% settings
[prm,Fmat] = perfusion1c.settings;
prm = perfusion1c.mergestruct(prm,prmin);

% basename
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);

pathload = ['results/synt-createflowTPFA-' basenameflow '.mat'];
msg = ['Loading ' pathload];
disp(msg);
D = load(pathload);

% create aif
aifval = perfusion1c.createaif(prm.aiftype,prm.timeline);

% run forward problem    
Cmat = perfusion1c.syntforwprobpde(D.phimat,D.qmat,D.Fmat,aifval,prm);



savpaper = 1;
savdata = 1;
savmov = 0;
basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
if savmov
    
    % save the concentration maps
    step = 1/prm.dt;

    % delete figs/synthetic-forwprob/*
    pathsave = ['figs' '/' 'synt-' mfilename '-' basenameindicator '.avi'];
    mat2movie(Cmat,pathsave,1,10,lim,step)
    
end;

if savpaper
    
    lim = [0,prm.scaling*max(aifval)];
    % make panel
    start = 0;
    starti = find(prm.timeline > start,1);
    n = round(numel(prm.timeline)/2) - starti;
    step = floor(n/16);
    c = starti;
    timeplot = zeros(16,1);
    C = 0;
    for i = 1 : 4
        for j = 1 : 4            
            C = C + 1;
            A{i,j} = Cmat(:,:,:,c);
            clim{i,j} = lim*D.phimat(1,1);
            timeplot(C) = prm.timeline(c);
            c = c + step;
        end;
    end;
    timeplot
    H = perfusion1c.panelstruct(A,0.01,700,'lim',clim);
    pathsave = ['figs' '/' 'synt-' mfilename '-'  basenameindicator '-panel.eps'];
    
    print(H,pathsave,'-deps')
    
    % makemovieraw(numberim,1,1,pathsave);
end

if savdata
    pathsave = ['results/synt-' mfilename '-' basenameindicator '.mat'];
    msg = ['Saving ' pathsave];
    disp(msg);
    save(pathsave,'aifval','Cmat','prm','-v7.3')
end;






