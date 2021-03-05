function p = vtk_fiber_read(file, varargin)
% Read fiber data from file
% Usage:
%   p = vtk_fiber_read(file)

% Open the file
fid = fopen(file, 'r');

ifiber=1;
p.points = [];
p.hdr.name = 'vtk fibers';


% Scan file until we find 'ElementSpacing array
spacing = [1 1 1];
while ~feof(fid)
    try
        val = textscan(...
            fid, 'ElementSpacing = %f %f %f\n',1, 'ReturnOnError', 0);
        spacing = double(cell2mat(val));
        break;
    catch
    end
end
        

% Read until we encounter a tube object
while ~feof(fid)
    
    % Read next line
    try
        val = textscan(fid, 'NPoints = %d\n', 1, 'ReturnOnError', 0);
    catch
        continue;
    end

    junk = fgets(fid);
    n = val{1};        

    if n > 3            
        X = [];
        k = 1 + size(p.points,1);            
        for i = 1:val{1}           
            pt = double(cell2mat(textscan(fid,'%f', 18)));
            x = pt(1:3)';
            if isempty(X) || any(x ~= X(end,:))   
                X = [X; spacing .* x];
            end
        end

        p.points = [p.points; X];
        p.cells.lines{ifiber} = [k:(k + size(X,1) - 1)]';
        ifiber = ifiber + 1;
        fprintf('.');
        if mod(ifiber,80) == 0
            fprintf(' %d \n',ifiber);
        end
    end
    
    % Check if it's a tube
end

fprintf('\n');
fclose(fid)