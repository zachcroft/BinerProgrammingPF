% This function prepares the Fourier transform 
% coefficients kx, ky, k2, and k4

function [kx,ky,k2,k4] = prepare_fft(Nx,Ny,dx,dy)
  
  format long;

  Nx21 = Nx/2 + 1;
  Ny21 = Ny/2 + 1;

  Nx2 = Nx + 2;
  Ny2 = Ny + 2;

  delkx = (2.0*pi) / (Nx*dx);
  delky = (2.0*pi) / (Ny*dy);

  % Calculate kx
  for i = 1:Nx21

    fk1         = (i-1) * delkx;
    kx(i)       = fk1;
    kx(Nx2 - i) = -fk1;
  
  end

  % Calculate ky
  for i = 1:Ny21

    fk2         = (i-1) * delky;
    ky(i)       = fk2;
    ky(Ny2 - i) = -fk2;
  
  end

  % Calcualte k2
  for i = 1:Nx
    for j = 1:Ny

      k2(i,j) = kx(i)^2 + ky(j)^2;

    end
  end

  % Calculate k4
  k4 = k2 .^ 2;

end