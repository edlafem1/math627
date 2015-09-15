function [lambda, x, it] = my_power (A, x, tol, itmax);
% Power method for approximation of largest eigenvalue of matrix A.
% On input, A is an n-by-n matrix, the n-vector x is the initial guess,
% tol is the tolerance, and itmax the maximum number of iterations.
% On output, lambda and x are the eigenvalue and eigenvector approximations,
% respectively, that is, A*x=lamda*x, and it is the number of iterations taken.

  % ... check sizes of x and A:
  [n, m] = size(x);
  if (m ~= 1)
    error('vector x must be column vector');
  end;
  if ( (n ~= size(A,1)) || (n ~= size(A,2)) )
    error('matrix A must be square and consistent with x');
  end;

  y = A * x;
  err = tol + 1.0; % to ensure one pass through the while loop
  it = 0;
  while ( (err > tol) & (it < itmax) )
    it = it + 1; % increment iteration counter

    if (it == 1) % special case of first iteration
      lambdaold = 0.0;
    else
      lambdaold = lambda;
    end;

    normy = norm(y); % Euclidean vector norm of vector y
    x = y / normy; % scale y by its norm, such that x has norm 1
    y =  A * x;
    lambda = x' * y; % = x'*A*x = (y'*A*y)/(y'*y) = eigenvalue approximation

    err = abs( (lambda-lambdaold) / lambda );
  end;

  if (it == itmax)
    warning ('maximum number of iterations reached');
  end;

