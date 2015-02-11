function [] = reducetimesteps(pathload,stepred)



% basename indicator file
msg = ['Loading ' pathload];
disp(msg);
E = load(pathload);

Cmat = E.Cmat(:,:,:,1:stepred:end);
aifval = E.aifval(1:stepred:end);
prm.timeline = E.prm.timeline(1:stepred:end);
prm.dt = prm.timeline(2)-prm.timeline(1);

msg = ['Timestep after reduction: ' num2str(prm.dt) 's' ];
disp(msg);
[a,b,c] = fileparts(pathload);
pathsave = [a '/' b '-red-' num2str(stepred) c];
clear a b c;
disp(['Saving ' pathsave]);
save(pathsave,'Cmat','aifval','prm');

% 
% % settings
% [prm,Fmat] = settings;
% prm = perfusion1c.mergestruct(prm,prmin);
% 
% step = prm.stepred;
% 
% % basename indicator file
% basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
% pathload = ['results/synt-createindicatorpde-' basenameindicator '.mat'];
% msg = ['Loading ' pathload];
% disp(msg);
% E = load(pathload);
% 
% Cmat = E.Cmat(:,:,:,1:step:end);
% aifval = E.aifval(1:step:end);
% prm.timeline = E.prm.timeline(1:step:end);
% prm.dt = prm.timeline(2)-prm.timeline(1);
% 
% msg = ['Timestep after reduction: ' num2str(prm.dt) 's' ];
% disp(msg);
% [a,b,c] = fileparts(pathload);
% pathsave = [a '/' b '-red-' num2str(step) c];
% clear a b c;
% disp(['Saving ' pathsave]);
% save(pathsave,'Cmat','aifval','prm');