function A = vtk_adjacency_matrix(mesh)
% Compute adjacency matrix from ITK mesh polygon spec
% Usage:
%   A = vtk_adjacency_matrix(mesh)

T=mesh.cells.polygons;

k=size(mesh.points,1);
n=0;
for i = 1:length(T)
  n = n + length(T{i});
end

A = spalloc(k,k,2*n);
for i = 1:length(T)
  t = T{i};
  for j = 1:length(t)
    A(t(1),t(2)) = 1;
    A(t(2),t(1)) = 1;
    t = circshift(t,1);
  end
end



