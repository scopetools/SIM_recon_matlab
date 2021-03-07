function [OTF_scaled]=resampleOTF(inputOTF,pxl_dim_PSF,inputData,pxl_dim_data,nphases,norders,norientations)
[ny_PSF,nx_PSF,nz_PSF,~,~] = size(inputOTF);
dk_PSF=1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF);

[ny_data,nx_data,nimgs_data] = size(inputData);
nz_data=nimgs_data/(nphases*norientations);
dk_data=1./([ny_data,nx_data,nz_data].*pxl_dim_data);

OTF_yy=[-ceil((ny_PSF-1)/2):floor((ny_PSF-1)/2)]*dk_PSF(1);
OTF_xx=[-ceil((nx_PSF-1)/2):floor((nx_PSF-1)/2)]*dk_PSF(2);
OTF_zz=[-ceil((nz_PSF-1)/2):floor((nz_PSF-1)/2)]*dk_PSF(3);

map_yy=[-ceil((ny_data-1)/2):floor((ny_data-1)/2)]*dk_data(1);
map_xx=[-ceil((nx_data-1)/2):floor((nx_data-1)/2)]*dk_data(2);
map_zz=[-ceil((nz_data-1)/2):floor((nz_data-1)/2)]*dk_data(3);

[OTF_xx_arr,OTF_yy_arr,OTF_zz_arr]=meshgrid(OTF_xx,OTF_yy,OTF_zz);
[map_xx_arr,map_yy_arr,map_zz_arr]=meshgrid(map_xx,map_yy,map_zz);

[OutputSizeX, OutputSizeY, OutputSizeZ] = size(map_xx_arr);                      % output sizes
OTF_scaled=complex(zeros(OutputSizeX, OutputSizeY, OutputSizeZ, norders, norientations)); % preallocate the full array for performance   
tic
fprintf('OTF interpolating');
for jj=1:norientations
    for ii=1:norders
        OTF_scaled(:,:,:,ii,jj) = interp3(OTF_xx_arr,OTF_yy_arr,OTF_zz_arr,inputOTF(:,:,:,ii,jj),map_xx_arr,map_yy_arr,map_zz_arr,'cubic',0+0i);
        fprintf('.');
    end
end
elapsed = toc;
fprintf('Resampled OTF to image (%d x %d x %d) at %d orientations, %d orders. Took %d seconds.\n', OutputSizeX, OutputSizeY, OutputSizeZ, norientations, norders,  round(elapsed) );
end