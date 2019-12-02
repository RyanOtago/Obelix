function [Poisson] = Poisson3D(A, B, KX, KY, Dealias)

% Takes Fourier transformed scalar functions A & B and calculates the
% Poisson bracket   {A,B} = [(del_x(A))(del_y(B)) - (del_y(A))(del_x(B)) ].
% This is returned as a dealiased matrix [Poisson] in Fourier Space.
% NOTE: A, B, KX, KY, Dealias must all be matrices of the same size

% A, B are the scalar functions we wish to calculate the bracket for
% KX, KY are the matrix of wavevectors defining the Fourier Space
% Dealias is a dealiasing matrix defined by the 2/3 rule

A_x = real(ifftn(KX.*A));
A_y = real(ifftn(KY.*A));
B_x = real(ifftn(KX.*B));
B_y = real(ifftn(KY.*B));

% Important to use 'real' in case numerical errors introduce a complex
% component through the Fourier Transform

Poisson_Real = A_x.*B_y - A_y.*B_x;
Poisson_Aliased = fftn(Poisson_Real);

Poisson = Poisson_Aliased.*Dealias;