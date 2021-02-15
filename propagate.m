function [] = propagate(fndcm,outdir,tag,seg_ref,seg_ref_vtk,framenums,fref)

% -------------------------------------------------------------------------
% Propagates a reference segmentation to other frames in a series. 
%
% INPUT:
%
%   fndcm: Cartesian DICOM filename
%   outdir: name of existing directory for output files
%   tag: study identifier (e.g., 'bav08_root')
%   seg_ref: filename of reference segmentation (nii)
%   seg_ref_vtk: filename of reference segmentation (vtk mesh)
%   framenums: list of frame numbers to segment
%   fref: reference segmentation frame
%
%   Use full paths for all filenames. The script contains system calls to
%   both c3d and vtklevelset.
% -------------------------------------------------------------------------

% series of frame numbers
f = framenums;

% index of reference frame in f
fref_ind_f = find(f == fref);

% convert reference segmentation (.nii) to a vtk mesh
% seg_ref_vtk = [strtok(seg_ref,'.') '.vtk'];
% system(['vtklevelset ' seg_ref ' ' seg_ref_vtk ' 1']);

% dilate the reference segmentation and downsample to create a mask image
% mask_ref_srs = [strtok(seg_ref,'.') '_srs.nii.gz'];
mask_ref_srs = [outdir '/mask' sprintf('%02d',fref) '_' tag '_srs.nii.gz'];
mask_ref_srs_vtk = [strtok(mask_ref_srs,'.') '.vtk'];
system(['/usr/local/bin/c3d -int 0 ' seg_ref ' -threshold 1 inf 1 0 -dilate 1 10x10x10vox -resample 50% -o ' mask_ref_srs])
system(['vtklevelset ' mask_ref_srs ' ' mask_ref_srs_vtk ' 1']);

% cells of input mesh
mref = vtk_polydata_read(seg_ref_vtk);
mcells = mref.cells;

% write series of 3D images from Cartesian DICOM 
for j = 1 : length(f)
    i = f(j);
    
    fnimg = [outdir '/img' sprintf('%02d',i) '_' tag '.nii.gz'];
    read_export_dicom(fndcm,i,fnimg);
    
    fnimg_rs = [outdir '/img' sprintf('%02d',i) '_' tag '_srs.nii.gz'];
    system(['/usr/local/bin/c3d ' fnimg ' -smooth 1mm -resample 50% -o ' fnimg_rs]);
end

% The frame fref is fixed and the frame i is moving. Apply the inverse 
% transform to propagate the reference segmentation to each moving image 
% in series. The following two loops operate on downsampled images and can 
% be run in parallel.

% propagate forward
warp_str_forward = [];
for j =  fref_ind_f : length(f) - 1
    i = f(j);
    i_next = f(j+1);

        if i == fref
            mask_init = mask_ref_srs;
        else 
            mask_init = [outdir '/mask' sprintf('%02d',f(j-1)) '_to_' sprintf('%02d',i) '_' tag '_srs_reslice_init.nii.gz'];
        end

        img_fix = [outdir '/img' sprintf('%02d',i) '_' tag '_srs.nii.gz'];
        img_mov = [outdir '/img' sprintf('%02d',i_next) '_' tag '_srs.nii.gz'];

        regout_affine = [outdir '/affine' sprintf('%02d',i_next) '_to_' sprintf('%02d',i) '_srs_init.mat'];
        regout_deform = [outdir '/warp' sprintf('%02d',i_next) '_to_' sprintf('%02d',i) '_srs_init.nii.gz'];
        regout_deform_inv = [outdir '/warp' sprintf('%02d',i_next) '_to_' sprintf('%02d',i) '_srs_init_inv.nii.gz'];

        img_reslice_init = [outdir '/img' sprintf('%02d',i) '_to_' sprintf('%02d',i_next) '_' tag '_srs_reslice_init.nii.gz'];
        mask_init_reslice = [outdir '/mask' sprintf('%02d',i) '_to_' sprintf('%02d',i_next) '_' tag '_srs_reslice_init.nii.gz'];
        mask_init_reslice_vtk = [outdir '/mask' sprintf('%02d',i) '_to_' sprintf('%02d',i_next) '_' tag '_srs_reslice_init.vtk'];

        greedy_call(img_fix,img_mov,[],regout_affine,regout_deform,regout_deform_inv,mask_init);

        warp_str_forward = [warp_str_forward ' ' regout_affine ',-1 ' regout_deform_inv ' '];
        apply_warp('grayscale',img_mov,img_fix,img_reslice_init,warp_str_forward,[],[]);
        apply_warp('label',img_mov,mask_ref_srs,mask_init_reslice,warp_str_forward,[],[]);
        apply_warp('mesh',img_mov,mask_ref_srs_vtk,mask_init_reslice_vtk,warp_str_forward,[],[]);
end

        
% propagate backward
warp_str_back = [];
for j = fref_ind_f : -1 : 2
    i = f(j);
    i_next = f(j-1);

        if i == fref
            mask_init = mask_ref_srs;
        else 
            mask_init = [outdir '/mask' sprintf('%02d',f(j+1)) '_to_' sprintf('%02d',i) '_' tag '_srs_reslice_init.nii.gz'];
        end

        img_fix = [outdir '/img' sprintf('%02d',i) '_' tag '_srs.nii.gz'];
        img_mov = [outdir '/img' sprintf('%02d',i_next) '_' tag '_srs.nii.gz'];

        regout_affine = [outdir '/affine' sprintf('%02d',i_next) '_to_' sprintf('%02d',i) '_srs_init.mat'];
        regout_deform = [outdir '/warp' sprintf('%02d',i_next) '_to_' sprintf('%02d',i) '_srs_init.nii.gz'];
        regout_deform_inv = [outdir '/warp' sprintf('%02d',i_next) '_to_' sprintf('%02d',i) '_srs_init_inv.nii.gz'];

        img_reslice_init = [outdir '/img' sprintf('%02d',i) '_to_' sprintf('%02d',i_next) '_' tag '_srs_reslice_init.nii.gz'];
        mask_init_reslice = [outdir '/mask' sprintf('%02d',i) '_to_' sprintf('%02d',i_next) '_' tag '_srs_reslice_init.nii.gz'];
        mask_init_reslice_vtk = [outdir '/mask' sprintf('%02d',i) '_to_' sprintf('%02d',i_next) '_' tag '_srs_reslice_init.vtk'];

        greedy_call(img_fix,img_mov,[],regout_affine,regout_deform,regout_deform_inv,mask_init);

        warp_str_back = [warp_str_back ' ' regout_affine ',-1 ' regout_deform_inv ' '];
        apply_warp('grayscale',img_mov,img_fix,img_reslice_init,warp_str_back,[],[]);
        apply_warp('label',img_mov,mask_ref_srs,mask_init_reslice,warp_str_back,[],[]);
        apply_warp('mesh',img_mov,mask_ref_srs_vtk,mask_init_reslice_vtk,warp_str_back,[],[]);
end

% Propagate the reference segmentation to each of the other images in
% series using the masks generated above. This loop is performed at full
% resolution.
for j = 1 : length(f)
    i = f(j);
    if i ~= fref        
        img_fix = [outdir '/img' sprintf('%02d',i) '_' tag '.nii.gz'];
        img_mov = [outdir '/img' sprintf('%02d',fref) '_' tag '.nii.gz'];
        
        affine_warps = [];
        affine_warps_pts = [];
        if j > fref_ind_f
            
            % string of affine transforms from fref to i
            for k = fref_ind_f+1 : j
                regout_affine_init = [outdir '/affine' sprintf('%02d',f(k)) '_to_' sprintf('%02d',f(k-1)) '_srs_init.mat'];
                affine_warps = [affine_warps ' ' regout_affine_init ',-1 '];
                affine_warps_pts = [regout_affine_init ' ' affine_warps_pts];
            end
            
            % mask for this frame
            mask_fix_srs = [outdir '/mask' sprintf('%02d',f(j-1)) '_to_' sprintf('%02d',i) '_' tag '_srs_reslice_init.nii.gz'];
            mask_fix = [outdir '/mask' sprintf('%02d',f(j-1)) '_to_' sprintf('%02d',i) '_' tag '_reslice_init.nii.gz'];
            system(['/usr/local/bin/c3d -interpolation NearestNeighbor ' img_fix ' ' mask_fix_srs ' -reslice-identity -o ' mask_fix]);
            
        else
            
            % string of affine transforms
            for k = fref_ind_f-1 : -1: j
                regout_affine_init = [outdir '/affine' sprintf('%02d',f(k)) '_to_' sprintf('%02d',f(k+1)) '_srs_init.mat'];
                affine_warps = [affine_warps ' ' regout_affine_init ',-1 '];
                affine_warps_pts = [regout_affine_init ' ' affine_warps_pts];
            end
                
             % mask for this frame
            mask_fix_srs = [outdir '/mask' sprintf('%02d',f(j+1)) '_to_' sprintf('%02d',i) '_' tag '_srs_reslice_init.nii.gz'];
            mask_fix = [outdir '/mask' sprintf('%02d',f(j+1)) '_to_' sprintf('%02d',i) '_' tag '_reslice_init.nii.gz'];
            system(['/usr/local/bin/c3d -interpolation NearestNeighbor ' img_fix ' ' mask_fix_srs ' -reslice-identity -o ' mask_fix]);
            
        end
       
        img_reslice = [outdir '/img' sprintf('%02d',fref) '_to_' sprintf('%02d',i) '_' tag '_reslice.nii.gz'];
        seg_reslice = [outdir '/seg' sprintf('%02d',fref) '_to_' sprintf('%02d',i) '_' tag '_reslice.nii.gz'];
        seg_reslice_vtk = [outdir '/seg' sprintf('%02d',fref) '_to_' sprintf('%02d',i) '_' tag '_reslice_labeled.vtk'];
        
        regout_deform = [outdir '/warp' sprintf('%02d',fref) '_to_' sprintf('%02d',i) '.nii.gz'];
        regout_deform_inv = [outdir '/warp' sprintf('%02d',fref) '_to_' sprintf('%02d',i) '_inv.nii.gz'];
        
        greedy_call(img_fix,img_mov,affine_warps,[],regout_deform,regout_deform_inv,mask_fix);
        apply_warp('grayscale',img_fix,img_mov,img_reslice,affine_warps,regout_deform,[]);
        apply_warp('label',img_fix,seg_ref,seg_reslice,affine_warps,regout_deform,[]);
        
        mesh_warps = [affine_warps_pts regout_deform_inv];
        apply_warp('mesh',img_fix,seg_ref_vtk,seg_reslice_vtk,mesh_warps,[],[]);
        restore_vtk_data(seg_reslice_vtk,mcells);
   
    end
end

end

function [] = restore_vtk_data(fnmesh,mcells)

    m = vtk_polydata_read(fnmesh);
    m.cells = mcells;
    vtk_polydata_write(fnmesh,m);

end
