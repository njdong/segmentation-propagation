function f = vtk_sample_image_util(fn_image, fn_mesh_in, fn_mesh_out, attr_name)
% Sample values from an image
% Usage:
%   f = vtk_sample_image_util(fn_image, fn_mesh_in, fn_mesh_out, attr_name)
mesh = vtk_polydata_read(fn_mesh_in);
A = spm_get_space(fn_image);
X = hippo_affine_inverse(A, mesh.points);
f = spm_sample_vol(spm_vol(fn_image), X(:,1), X(:,2), X(:,3), 1);
meshout = vtk_add_point_data(mesh, attr_name, f);
vtk_polydata_write(fn_mesh_out, meshout);
