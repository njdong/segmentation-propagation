function [ ] = dcm_2_mcnifti( fndicom, fnout )

% ------------------------------------------------------------------------
% INPUT:
%
%   fndicom is the filename of the Cartesian dicom ('something.dcm')
%   
%   fnout is the filename of the output multi-component nifti image
%       ('something.nii.gz')
%
% ------------------------------------------------------------------------

% BIN = '/Applications/usr/bin/'; 

% Read the Cartesian dicom
I = readDicom3D(fndicom);

% Image resolution
x_dim = I.width;
y_dim = I.height;
z_dim = I.depth;

delta_x = I.widthspan*10/x_dim;    % resolution in x (mm)
delta_y = I.heightspan*10/y_dim;   % resolution in y (mm)
delta_z = I.depthspan*10/z_dim;    % resolution in z (mm)

% Create and save the 4D nifti
nii = make_nii(I.data,[delta_x delta_y delta_z],[1 1 1],2);
if ~strcmp(fnout(end-2:end),'.gz')
    save_nii(nii,fnout);
    gzip_img(fnout,'gzip');
else 
    save_nii(nii,fnout(1:end-3));
    gzip_img(fnout(1:end-3),'gzip');
end

end

function fnimg_new = gzip_img(fnimg,zip_flag)

% -------------------------------------------------------------------------
% Input:    fnimg is the filename of the image to be written
%           zip_flag is either 'gunzip' or 'gzip'
%
% Output:   fnimg_new is the new filename of the image to be written with
%                   the correct extension
% -------------------------------------------------------------------------

    if strcmp(zip_flag,'gunzip')

        % only unzip if image is zipped
        if strcmp(fnimg(end-2:end),'.gz') 
            system(['gunzip ' fnimg]);
            fnimg_new = fnimg(1:end-3);
        else fnimg_new = fnimg;
        end

    elseif strcmp(zip_flag,'gzip')

        % only zip if image is unzipped
        if ~strcmp(fnimg(end-2:end),'.gz') 
            system(['gzip ' fnimg]);
            fnimg_new = [fnimg '.gz'];
        else fnimg_new = fnimg;
        end

    else error('myApp:flagChk','Invalid flag');
    end
end

