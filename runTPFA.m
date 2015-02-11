nx = 32;
ny = 32;
for i = 1 : 5
    nx = nx*2;
    ny = ny*2;
    
    p = [10,10,1];
    
    % Gimse code
    Grid.Nx=nx; Grid.hx=p(1)/Grid.Nx;
    Grid.Ny=ny; Grid.hy=p(2)/Grid.Ny;
    Grid.Nz=1; Grid.hz=p(3)/Grid.Nz;
    Grid.K=ones(3,Grid.Nx,Grid.Ny);
    N=Grid.Nx*Grid.Ny*Grid.Nz; 
    q=zeros(N,1); 
    q([1 N])=[1 -1];
    physdim = prod(p);
    voxelvol = Grid.hx*Grid.hy*Grid.hz;
    % q = q/voxelvol;
    
    'NEW'

    [P,V]=TPFAgimse(Grid,Grid.K,q,voxelvol);
    P = P - min(P(:));
    [nx, ny]
    mid = nx/2;
    a = V.x(mid,:);
    sum(a)
    

    % My code
    dim = [Grid.Nx,Grid.Ny,Grid.Nz];
    Kmat = zeros([dim,3]);
    Kmat(:,:,:,1) = squeeze(Grid.K(1,:,:,:));
    Kmat(:,:,:,2) = squeeze(Grid.K(2,:,:,:));
    Kmat(:,:,:,3) = squeeze(Grid.K(3,:,:,:));
    Qmat = reshape(q,Grid.Nx,Grid.Ny);
    prmin.h = [Grid.hx,Grid.hy,Grid.hz];
    prmin.dim = dim;
    prmin.mu = 1;
    [pmat,qmat] = perfusion1c.syntTPFA(Kmat,Qmat,prmin);
    
    % [P,V]=TPFAoriginal(Grid,Grid.K,q);
    pmat = pmat - min(pmat(:));
    [nx, ny]
    mid = nx/2;
    a = qmat{1}(mid,:);
    sum(a)
    pause
    
end;
