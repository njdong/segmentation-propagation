function p1 = vtk_set_point_data(p, name, data)
% Update an attribute array in a vtk mesh
% Usage:
%   p1 = vtk_set_field_data(p, name, data)
% Parameters
%   p        VTK mesh struct (from vtk_polydata_read)
%   name     Name of the new array (string)
%   data     An Nxk matrix of values to add

arr.name = name;
arr.type = 'field';

if size(data, 1) == size(p.points,1)
    arr.data = data;
elseif size(data, 2) == size(p.points,1)
    arr.data = data';
else
    error('Data size does not match point array size');
end


p1 = p;
if ~isfield(p1, 'point_data')
    p1.point_data(1) = arr;
else
    id = find(strcmpi(name, {p.point_data(:).name}));
    if isempty(id)
        p1.point_data(length(p1.point_data)+1) = arr;
    else
        p1.point_data(id) = arr;
    end
end
