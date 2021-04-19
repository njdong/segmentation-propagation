fndcm = '/users/jileihao/playground/sandbox/bav07.dcm';
outdir = '/users/jileihao/playground/sandbox/matlab';
tag = 'matlab';
seg_ref = '/users/jileihao/playground/sandbox/seg05_bav07_root_labeled.nii.gz';
seg_ref_vtk = '/users/jileihao/playground/sandbox/seg05_bav07_root_labeled.vtk';
framenums = [2 3 5 7 9 11];
fref = 5;

addpath('vtk_scripts')
addpath('cartesian_dicom_scripts')

propagate(fndcm, outdir, tag, seg_ref, seg_ref_vtk, framenums, fref);


