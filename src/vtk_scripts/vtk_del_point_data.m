function p = vtk_del_point_data(p, name)
% Delete a point attribute array from a vtk mesh
% Usage:
%   p = vtk_del_field_data(p, name)
% Parameters
%   p         VTK mesh struct (from vtk_polydata_read)
%   name      Name of the new array (string)



if ~isfield(p, 'point_data')
    error('No point data in mesh');
else
    pos = strmatch(name, {p.point_data.name}, 'exact');
    if ~isempty(pos)
       p.point_data(pos(1))=[];
    else         
       error('No array by that name');
    end
end
