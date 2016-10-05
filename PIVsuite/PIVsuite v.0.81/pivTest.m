example_01_Image_pair_simple;

X = pivData.X;
Y = pivData.Y;
U = pivData.U;
V = pivData.V;

U = inpaint_nans(U);
V = inpaint_nans(V);

ppp.iaSizeX = 32;
ppp.iaSizeY = 32;
ppp.iaStepX = 16;
ppp.iaStepY = 16;

IAinfo = piv2GetIAinfo(size(imread(im1)),ppp);

tic;
[Xia,Yia,Uia,Via] = arrayfun(@(XXX) piv2InterpolateDisplacement(XXX,X,Y,U,V),IAinfo,'UniformOutput',false);
toc


Xi = cell(size(IAinfo,1),size(IAinfo,2));
Yi = cell(size(IAinfo,1),size(IAinfo,2));
Ui = cell(size(IAinfo,1),size(IAinfo,2));
Vi = cell(size(IAinfo,1),size(IAinfo,2));

tic;
for kx = 1:size(IAinfo,2);
    for ky=1:size(IAinfo,1);
        [Xi{ky,kx},Yi{ky,kx},Ui{ky,kx},Vi{ky,kx}]=piv2InterpolateDisplacement(IAinfo(ky,kx),X,Y,U,V);
    end;
end;
toc


tic;
% gpuIAinfo = gpuArray(IAinfo);
gpuX = gpuArray(X);
gpuY = gpuArray(Y);
gpuU = gpuArray(U);
gpuV = gpuArray(V);
[gpuXia,gpuYia,gpuUia,gpuVia] = arrayfun(@(XXX) piv2InterpolateDisplacement(XXX,X,Y,U,V),IAinfo,'UniformOutput',false);
Xia = gather(Xia);
Yia = gather(Yia);
Uia = gather(Uia);
Via = gather(Via);
toc

