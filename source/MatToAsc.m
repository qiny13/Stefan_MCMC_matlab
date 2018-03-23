% Author: Qin | Date: 2015-09-09
% MatToAsc: Read and Transform *.mat files in Matlab and output *.asc grid-files
function MatToAsc(matrix,filename,cols,rows,xllc,yllc,size,nodata)
% e.g. MatToAsc(INPUT, 'OUTPUT', 20, 10, 2000, -1500, 1000, -9999)
% Create and Open the asc file
filename = [filename '.asc'];
fid = fopen(filename,'w');
% Write the table header
fprintf(fid,'%s         %d\n', 'ncols', cols); % num of cols
fprintf(fid,'%s         %d\n', 'nrows', rows); % num of rows
fprintf(fid,'%s     %d\n', 'xllcorner', xllc); % coordinates of x left-low
fprintf(fid,'%s     %d\n', 'yllcorner', yllc); % coordinates of y left-low
fprintf(fid,'%s      %d\n', 'cellsize', size); % size of grid cell
fprintf(fid,'%s  %d\n', 'NODATA_value', nodata); % nodata value
% Input the matrix from *.mat file into the *.asc file
dlmwrite(filename, matrix, '-append', 'delimiter', ' ');
% Close the *.asc file
fclose(fid);
end