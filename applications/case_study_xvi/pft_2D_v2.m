%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2D Semi-implicit Spectral    %
%    Phase-field Crystal Code     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include path to functions
addpath('../../functions/');

% Get initial wall time
time0 = clock();
format long;


out1 = fopen('final_conf.out','w');
out2 = fopen('energy.out','w');

% Simulation cell parameters
Nx = 64;
Ny = 64;

pix = 4.0*atan(1.0);

dx = pix / 4.0;
dy = pix / 4.0;

% Time integration parameters
nsteps = 2e5;
dtime  = 0.01;
nprint = 2000;
nstart = 1;

% Material parameters
den0   = -0.285;
tempr  = -0.25;
tempr0 = tempr;
noise  = 0.01 * den0;


% Set initial conditions
nfile = 0;
if nfile==1
  % Read from input file
  in1 = fopen('g3_2r.inp','r');
  % ... put code for reading data here
else
    % Define initial conditions
  for i = 1:Nx
    for j = 1:Ny
      den(i,j) = den0 + noise*(0.5 - rand);
    end
  end
end


% Prepare FFT
[kx,ky,k2,k4] = prepare_fft(Nx,Ny,dx,dy);


% Evolve
for istep = nstart:nsteps
  %tempr = tempr + tempr0 / nsteps;

  % Take the FFT of the density
  f_den = fft2(den);

  Linx = -k2 .* (tempr + 1.0 - 2.0*k2 + k4);

  denom = 1.0 - dtime .* Linx;

  den3 = den.^3;

  f_den3 = fft2(den3);

  Nonx = -k2 .* f_den3;

  f_den = (f_den + dtime * Nonx) ./ denom;

  den = real(ifft2(f_den));


  % Print results
  if ((mod(istep,nprint) == 0) || (istep == 1))
    fprintf('done step: %5d\n', istep);

    % Energy calculation
    ss2 = den.^2;
    ss4 = den.^4;

    f_ff = 0.5 * f_den .*(1.0 - 2.0*k2 + k4);

    ff = real(ifft2(f_ff));

    ff = ff .* den + 0.5 * tempr * ss2 + 0.25 * ss4;

    energy = sum(ff(:)) / (Nx*Ny);

    fprintf(out2, '%d %14.6e\n',istep,energy);

    % Export to vtk
    write_vtk_grid_values(Nx, Ny, dx, dy, istep, den);
    %write_vtk_grid_values(Nx, Ny, dx, dy, istep, den, ff);

    % Visualize
    imagesc(den)
    title(strcat('t = ', num2str(istep)));
    axis square
    colorbar
    drawnow
  end

  % Write Configuration
  if (istep == nsteps)
    for i = 1:Nx
      for j = 1:Ny
        fprintf(out1, '%5d %5d %14.6e\n', i, j, den(i,j));
      end
    end
  end

end % End of time step

compute_time = etime(clock(),time0);
fprintf('Compute time: %10d\n',compute_time);





