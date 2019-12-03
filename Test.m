LX = 1;
LY = 1;
NX = 1024;
NY = 1024;
I = sqrt(-1);

a = linspace(0,LX,NX);
b = linspace(0,LY,NY);

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
[KX, KY] = ndgrid(kx, ky);
dealias = abs(KX) < (1/3)*NX & abs(KY) < (1/3)*NY;

[A,B] = ndgrid(a,b);

Aexp = exp(A.^2 + B.^2);
Bsin = sin(A + B);

Ahat = fft2(Aexp);
Bhat = fft2(Bsin);
PBhat = Poisson(Ahat, Bhat, KX, KY);

%analytic = 2*

PB = real(ifft2(PBhat).*dealias);
plot3(a,b,Aexp)