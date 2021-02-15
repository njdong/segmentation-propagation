function arr = vtk_get_cell_data(p, name)
% Get the cell_data array by a certain name from a structure
% returned by vtk_mesh_read
% Usage:
%       arr = vtk_get_cell_data(p, name)
index = find(strcmpi(name, {p.cell_data(:).name}));
if length(index) ~= 1
    error('cell data not found');
end
arr = p.cell_data(index).data;

