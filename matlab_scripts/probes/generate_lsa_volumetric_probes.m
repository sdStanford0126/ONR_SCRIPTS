function grid = generate_lsa_volumetric_probes(output_file, write_grid)
%GENERATE_LSA_VOLUMETRIC_PROBES Define a structured LSA sampling volume.
%   grid = generate_lsa_volumetric_probes() previews the definition only.
%   grid = generate_lsa_volumetric_probes(output_file, true) writes the file.
%
%   h = 1, L = 2*h, dx = 0.05*h, dy = dz = 0.02*h.
%   x spans [0, 35*h]. Both y and z cover [-8.8*De, 8.8*De],
%   rounded OUTWARD to multiples of their spacing about zero. Thus the
%   grid includes zero and is exactly reflection-symmetric in y and z.
%
%   Text format: one descriptive comment header, then x y z columns.
%   Row order: x fastest, then y, then z, all increasing. In particular,
%       row = i + Nx*(j-1) + Nx*Ny*(k-1)
%   corresponds to (x(i), y(j), z(k)). For sequential sampled values q:
%       Q = reshape(q, [grid.Nx, grid.Ny, grid.Nz]);  % Q(i,j,k)
%   In NumPy:
%       Q = q.reshape((Nx, Ny, Nz), order='F')
%   These assume the sampling solver preserves the probe-file row order.
%   For coordinates read back into MATLAB:
%       xyz = readmatrix(output_file, 'NumHeaderLines', 1);
%       X = reshape(xyz(:,1), [grid.Nx, grid.Ny, grid.Nz]);
%       Y = reshape(xyz(:,2), [grid.Nx, grid.Ny, grid.Nz]);
%       Z = reshape(xyz(:,3), [grid.Nx, grid.Ny, grid.Nz]);
%
%   Only coordinate vectors are created in preview mode. Writing streams
%   small blocks to disk; the full point cloud is never allocated.

if nargin < 1 || isempty(output_file)
    output_file = fullfile(fileparts(mfilename('fullpath')), ...
        'lsa_volumetric_probes.txt');
end
if nargin < 2
    write_grid = false;
end
validateattributes(write_grid, {'logical'}, {'scalar'}, mfilename, 'write_grid');
if ~(ischar(output_file) && isrow(output_file)) && ...
        ~(isstring(output_file) && isscalar(output_file))
    error('LSAProbes:OutputPath', 'output_file must be a character row or string scalar.');
end
output_file = char(output_file);

h = 1;
L = 2*h;
dx = 0.05*h;
dy = 0.02*h;
dz = 0.02*h;
De = 1.3*(L*h)^0.625/(L+h)^0.25;
requested_half_extent = 8.8*De;

nx_intervals = round(35*h/dx);
ny_half = ceil(requested_half_extent/dy);
nz_half = ceil(requested_half_extent/dz);
x = (0:nx_intervals)*dx;
y = (-ny_half:ny_half)*dy;
z = (-nz_half:nz_half)*dz;
Nx = numel(x);
Ny = numel(y);
Nz = numel(z);
point_count = double(Nx)*double(Ny)*double(Nz);

header = sprintf(['# LSA volumetric probes; columns=x y z; h=%.8f; L=%.8f; ' ...
    'De=%.12g; dx=%.8f; dy=%.8f; dz=%.8f; ' ...
    'x=[%.8f,%.8f]; y=[%.8f,%.8f]; z=[%.8f,%.8f]; ' ...
    'requested_yz_extent=+/-8.8De (rounded outward); ' ...
    'Nx=%d; Ny=%d; Nz=%d; N=%.0f; ' ...
    'order=x-fastest,y-next,z-slowest; MATLAB_reshape=[Nx,Ny,Nz]'], ...
    h, L, De, dx, dy, dz, x(1), x(end), y(1), y(end), z(1), z(end), ...
    Nx, Ny, Nz, point_count);

% Exact text-size estimate for the fixed format below: each axis value is
% repeated for every combination of the other two axes, plus two spaces
% and one LF newline per row. Include the header and its LF newline.
x_characters = numel(sprintf('%.8f', x));
y_characters = numel(sprintf('%.8f', y));
z_characters = numel(sprintf('%.8f', z));
estimated_bytes = double(x_characters)*Ny*Nz + ...
    double(y_characters)*Nx*Nz + double(z_characters)*Nx*Ny + ...
    3*point_count + numel(header) + 1;

grid = struct('h', h, 'L', L, 'De', De, 'dx', dx, 'dy', dy, 'dz', dz, ...
    'x', x, 'y', y, 'z', z, 'Nx', Nx, 'Ny', Ny, 'Nz', Nz, ...
    'point_count', point_count, 'requested_half_extent', requested_half_extent, ...
    'estimated_bytes', estimated_bytes, 'output_file', output_file, ...
    'header', header, 'order', 'x-fastest,y-next,z-slowest', 'written', false);

fprintf('De/h = %.12g; requested transverse half-extent/h = %.12g\n', De/h, requested_half_extent/h);
fprintf('x/h: [%g, %g]; y/h: [%g, %g]; z/h: [%g, %g]\n', ...
    x(1)/h, x(end)/h, y(1)/h, y(end)/h, z(1)/h, z(end)/h);
fprintf('[Nx, Ny, Nz] = [%d, %d, %d]; total probes = %.0f\n', Nx, Ny, Nz, point_count);
fprintf('Text size: %.3f GB (%.3f GiB)\n', estimated_bytes/1e9, estimated_bytes/2^30);
if ~write_grid
    fprintf('Preview only: no grid file written.\n');
    return
end

if isfile(output_file)
    error('LSAProbes:FileExists', 'Output already exists: %s', output_file);
end
[fid, message] = fopen(output_file, 'w');
if fid == -1
    error('LSAProbes:OpenFailed', 'Cannot open %s: %s', output_file, message);
end
close_file = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', header);

y_per_block = 64;
progress_interval = max(1, ceil(Nz/20));
for k = 1:Nz
    for j_start = 1:y_per_block:Ny
        j_end = min(j_start + y_per_block - 1, Ny);
        % ndgrid + column-major flattening keeps x fastest across blocks.
        [X, Y] = ndgrid(x, y(j_start:j_end));
        block = [X(:).'; Y(:).'; repmat(z(k), 1, numel(X))];
        fprintf(fid, '%.8f %.8f %.8f\n', block);
    end
    [message, error_number] = ferror(fid);
    if error_number ~= 0
        error('LSAProbes:WriteFailed', 'Writing failed: %s', message);
    end
    if mod(k, progress_interval) == 0 || k == Nz
        fprintf('Written z planes: %d/%d (%.1f%%)\n', k, Nz, 100*k/Nz);
    end
end
% Close now, also triggering the cleanup handler on normal completion.
clear close_file
grid.written = true;
fprintf('Saved %s\n', output_file);
end
