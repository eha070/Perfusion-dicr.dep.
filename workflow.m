%%
% create the flow
perfusion1c.createflowTPFA;

%%

% create the timelapse of indicator
clear prmin;
aiftype = {'parker','gamma'};
% aiftype = {'gamma'};
for i = 1 : numel(aiftype)
    prmin.aiftype = aiftype{i};
    perfusion1c.createindicatorpde(prmin);
end;

%%

% Create the convolution flow
clear prmin;
aiftype = {'parker','gamma'};
% aiftype = {'gamma'};
for i = 1 : numel(aiftype)
    prmin.aiftype = aiftype{i};

    % create indicator using convolution
    perfusion1c.createindicatorconv(prmin);

end;


%%

% reduce time steps

% parameters
[prm,temp] = settings;

clear prmin;
aiftype = {'parker','gamma'};
% aiftype = {'gamma'};
stepred = [25,100,500];
for i = 1 : numel(aiftype)
    for j = 1 : numel(stepred)
        prm.aiftype = aiftype{i};
        
        basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
        pathload = ['results/synt-createindicatorpde' '-' basenameindicator '.mat'];                
        perfusion1c.reducetimesteps(pathload,stepred(j));
                
        basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
        pathload = ['results/syntconv-createindicatorconv' '-' basenameindicator '.mat'];                
        perfusion1c.reducetimesteps(pathload,stepred(j));
        
    end;
end;

return;
%%

%
% Recognize the flow using classical methods
%

% the options to run over
aiftype = {'parker','gamma'};
% deconvmethod = {'circSVD'};
stepred = [25,100,500];
datatype = {'synt','syntconv'};


clear prmin;
for i = 1 : numel(aiftype)
    for j = 1 : numel(deconvmethod)
        for k = 1 : numel(stepred)
            for l = 1 : numel(datatype)
                prmin.aiftype = aiftype{i};
                prmin.reconflowdeconvmethod = deconvmethod{j};
                prmin.stepred = stepred(k);
                prmin.reconflowdata = datatype{l};
                perfusion1c.reconflow(prmin);
            end;
        end;
    end;
end;

%%

%
% Load results for making table
%

% the optionts to run over
aiftype = {'parker','gamma'};
deconvmethod = {'circSVD'};
stepred= [100,500];
datatype = {'synt','syntconv'};

% load true perfusion
[prm,Fmat] = settings;
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);
pathload = ['results/synt-createflowTPFA-' basenameflow '.mat'];
msg = ['Loading true solution ' pathload];
disp(msg);
true = load(pathload);


clear A;
clear header
for l = 1 : numel(datatype)
    c = 0;
    for i = 1 : numel(aiftype)
        for j = 1 : numel(deconvmethod)
            for k = 1 : numel(stepred)
            
                c = c + 1;
                [prm,Fmat] = settings;
                
                prm.aiftype = aiftype{i};
                prm.reconflowdeconvmethod = deconvmethod{j};
                prm.stepred = stepred(k);

                basenamerecon = perfusion1c.providenamerecon(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T,prm.reconflowdeconvmethod,prm.stepred);
                pathload = ['results/' datatype{l} '-reconflow'  '-method-'  prm.reconflowcode  '-' basenamerecon  '.mat'];
                msg = ['Loading ' pathload];
                disp(msg);
                recon = load(pathload);

                val = 100*(recon.perfmatn - true.perfmatn)./true.perfmatn;
                val = abs(val);
                A.(datatype{l})(1,c) = mean(val(:));
                    
                val = 100*(recon.phimat - true.phimat)./true.phimat;
                val = abs(val);
                A.(datatype{l})(2,c) = mean(val(:));
                header{1,c} = [prm.aiftype '-' prm.reconflowdeconvmethod '-' num2str(prm.stepred)];
            end;
        end;
    end;
end;

p = 30;
for l = 1 : numel(datatype)
    A.(datatype{l}) = perfusion1c.mat2celldirect(A.(datatype{l}));
    desc = {'';'Perfusion';'Porosity'};
    pathsave = ['stat-' datatype{l} '.txt'];
    msg = ['Saving ' pathsave];
    disp(msg);
    fid = fopen(pathsave,'wt');
    B = [header;A.(datatype{l})];
    B = [desc,B];
    printcell(fid,B,p);
    fclose(fid);
    pathsave = ['stat-' datatype{l} '.tex'];
    msg = ['Saving ' pathsave];
    disp(msg);    
    perfusion1c.cell2tex(B,pathsave,2);
end;

