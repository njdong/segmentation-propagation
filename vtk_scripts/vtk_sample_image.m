function f = vtk_sample_image(fn_image, mesh)
% Sample values from an image
% Usage:
%   f = vtk_sample_image(fn_image, mesh)
A = spm_get_space(fn_image);
X = hippo_affine_inverse(A, mesh.points);
f = spm_sample_vol(spm_vol(fn_image), X(:,1), X(:,2), X(:,3), 1);

