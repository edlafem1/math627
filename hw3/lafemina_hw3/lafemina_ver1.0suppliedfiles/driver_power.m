
  n = 4;
  [I, J] = ndgrid([1:n]);
  A = 1 ./ (I + J - 1);
  x = ones(n,1) / sqrt(n); % set standard initial guess
  tol = 1.0e-12;
  itmax = 50;

  tic
  [lambda, x, it] = my_power (A, x, tol, itmax);
  tsec = toc;

  resnormabs = norm (A * x - lambda * x);
  resnormrel = resnormabs / lambda; % norm(x) = 1 by construction!

  % output:
  format compact
  format long e
  n
  itmax
  it
  tol
  lambda
  tsec
  if (n < 25) % output x only if not too long
    x
  end;
  disp(' ');

