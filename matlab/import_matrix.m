function A = mm_reader(filename);

  fid = fopen(filename, 'r');
  
  dummy_string = fscanf(fid, '%s %s %s %s %s', 5);
  
  size = fscanf(fid, '%d', 1);
  size = fscanf(fid, '%d', 1);
  nz = fscanf(fid, '%d', 1);
  
  disp(sprintf('\n Size: %i', size));
  disp(sprintf(' NZel: %i\n', nz));
  
  [val, count] = fscanf(fid, '%d %d %f', [3, nz]);
  
  A = sparse(val(1,:), val(2,:), val(3,:), size, size, nz);
  
  fclose(fid);
