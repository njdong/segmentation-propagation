function [] = read_export_dicom(fndcm,framenum,fnout)

%------------------------------------------------------------
% INPUT:
% 
% fndcm is the name of the 4D dicom file (something.dcm)
% framenum is the time point we want to export
% fnout is the name of the output image (something.nii)
%
%------------------------------------------------------------

% Read image
I = readDicom3D(fndcm);

% Get image dimensions
x_dim = I.width;
y_dim = I.height;
z_dim = I.depth;

delta_x = I.widthspan*10/x_dim;    % resolution in x (mm)
delta_y = I.heightspan*10/y_dim;   % resolution in y (mm)
delta_z = I.depthspan*10/z_dim;    % resolution in z (mm)

% Create and save the nifti image
nii = make_nii(I.data(:,:,:,framenum),[delta_x delta_y delta_z],[1 1 1],2);
if strcmp(fnout(end-2:end),'.gz')
    fnout = fnout(1:end-3);
end
save_nii(nii,fnout);
gzip_img(fnout,'gzip');