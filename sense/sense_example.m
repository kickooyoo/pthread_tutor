% sense_example.m

% make some image dimensions
nx = 128;
ny = 128;
nz = 128;
ncoil = 32;

x = randn(nx,ny,nz) + 1i*randn(nx,ny,nz);
smap = randn(nx,ny,nz,ncoil) + 1i*randn(nx,ny,nz,ncoil);

% do the matlab method
tic;
y1 = zeros(nx,ny,nz,ncoil);

for i=1:ncoil
    y1(:,:,:,i) = fftshift(fftn(ifftshift(x.*smap(:,:,:,i))));
end
toc;

% do the mex method
tic;
y2 = reshape(sense_example_mex(x, smap, int32(nx), int32(ny), int32(nz), ...
    int32(ncoil), int32(feature('numcores')), int32(0)), [nx ny nz ncoil]);
toc;

printf('normdiff is %f',norm(y2(:)-y1(:)));
