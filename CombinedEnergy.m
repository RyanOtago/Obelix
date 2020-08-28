Ezptot = [Ezp Ezp1];
Ezmtot = [Ezm Ezm1];
Esptot = [Esp Esp1];
Esmtot = [Esm Esm1];
ttot   = [t t(end)+t1];


subplot(2,1,1)
plot(ttot, Ezptot, ttot, Ezmtot)
title('\beta = 0.03')
legend('\zeta^+', '\zeta^-', 'Location', 'best')
ylabel('Alfven Wave Energy')

subplot(2,1,2)
plot(ttot, Esptot, ttot, Esmtot)
legend('z^+', 'z^-', 'Location', 'best')
xlabel('Time (t/t_A)')
ylabel('Slow Wave Energy')