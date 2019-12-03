LX = 2*pi;
LY = 2*pi;
NX = 128;
NY = 128;
I = sqrt(-1);

a = linspace(0,LX,NX);
b = linspace(0,LY,NY);

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
[KX, KY] = ndgrid(kx, ky);
dealias = abs(KX) < (1/3)*NX & abs(KY) < (1/3)*NY;

[A,B] = ndgrid(a,b);

Asin = sin(A + B);
Bsin = sin(A - B);

Ahat = fft2(Asin);
Bhat = fft2(Bsin);

PBhat = Poisson(Ahat, Bhat, KX, KY);

analytic = -cos(2*B)-cos(2*A);

PB = real(ifft2(PBhat.*dealias));

subplot(1,2,1)
surf(PB)
pbaspect([1 1 1])
subplot(1,2,2)
surf(analytic-PB)
pbaspect([1 1 1])