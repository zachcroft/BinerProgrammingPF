% This function writes the grid point values
% to a vtk file to generate the contour plots
% to be viewed by using Paraview/VisIt

function [] = write_vtk_grid_values(nx,ny,dx,dy,istep,data1)
  
  format long;

  % Open the output file
  fname = sprintf('time_%d.vtk', istep);
  out   = fopen(fname, 'w');

  nz = 1;

  npoin = nx*ny*nz;

  % Start writing the vtk ASCII file...

  % Header:
  fprintf(out, '# vtk DataFile Version 2.0\n');
  fprintf(out, 'time_10.vtk\n');
  fprintf(out, 'ASCII\n');
  fprintf(out, 'DATASET STRUCTURED_GRID\n');

  % Coordinates of grid points:
  fprintf(out,'DIMENSIONS %5d %5d %5d\n', nx, ny, nz);
  fprintf(out, 'POINTS %7d float\n', npoin);

  for i = 1:nx
    for j = 1:ny
      
      x = (i-1) * dx;
      y = (j-1) * dy;
      z = 0.0;

      fprintf(out, '%14.6e %14.6e %14.6e\n', x,y,z);

    end
  end

  % Write grid point values:
  fprintf(out,'POINT_DATA %5d\n', npoin);
  fprintf(out,'SCALARS CON float 1\n');
  fprintf(out,'LOOKUP_TABLE default\n');

  for i = 1:nx
    for j = 1:ny
      
      ii = (i-1) * nx + j;
      
      fprintf(out,'%14.6e\n', data1(i,j));

    end
  end

  fclose(out);

end
