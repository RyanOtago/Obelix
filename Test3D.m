LX = 2*pi;
LY = 2*pi;
LZ = 2*pi;

NX = 64;
NY = 64;
NZ = 64;

I = sqrt(-1);

a = linspace(0,LX-LX/NX,NX);
b = linspace(0,LY-LY/NY,NY);
c = linspace(0,LZ-LZ/NZ,NZ);

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
kz = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];
[KX, KY, KZ] = ndgrid(kx, ky, kz);
dealias = abs(KX)<(1/3)*NX & abs(KY)<(1/3)*NY & abs(KZ)<(1/3)*NZ;

[A,B,C] = ndgrid(a,b,c);

Asin = sin(A + B + C);
Bsin = sin(A - B - C);

Ahat = fftn(Asin);
Bhat = fftn(Bsin);

PBhat = Poisson(Ahat, Bhat, KX, KY);

analytic = -2.*cos(A+B+C).*cos(A-B-C);

PB = real(ifftn(PBhat.*dealias));

subplot(1,2,1)
contourslice(PB,[],[],[1 2 3 4])

subplot(1,2,2)
contourslice(analytic-PB,[],[],[1 2 3 4])
