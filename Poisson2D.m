function [Poisson] = Poisson2D(A, B, KX, KY, Dealias)

% Takes Fourier transformed scalar functions A & B and calculates the
% Poisson bracket   {A,B} = [(del_x(A))(del_y(B)) - (del_y(A))(del_x(B)) ].
% This is returned as a dealiased matrix [Poisson] in Fourier Space.
% NOTE: A, B, KX, KY, Dealias must all be matrices of the same size

% A, B are the scalar functions we wish to calculate the bracket for
% KX, KY are the matrix of wavevectors defining the Fourier Space
% Dealias is a dealiasing matrix defined by the 2/3 rule

A_x = real(ifft2(KX.*A));
A_y = real(ifft2(KY.*A));
B_x = real(ifft2(KX.*B));
B_y = real(ifft2(KY.*B));

% Important to use 'real' in case numerical errors introduce a complex
% component through the Fourier Transform

Poisson_Real = A_x.*B_y - A_y.*B_x;
Poisson_Aliased = fft2(Poisson_Real);

Poisson = Poisson_Aliased.*Dealias;