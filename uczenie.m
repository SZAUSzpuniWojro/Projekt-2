%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.577570e+003; foe(n+1)=4.663067e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=2.535681e+003; foe(n+1)=5.097471e+003; krok(n+1)=8.062699e-004; ng(n+1)=4.904469e+003;
n=2; farx(n+1)=1.634584e+003; foe(n+1)=1.318306e+004; krok(n+1)=1.802468e-004; ng(n+1)=7.522967e+003;
n=3; farx(n+1)=4.397440e+002; foe(n+1)=6.506080e+003; krok(n+1)=3.942116e-003; ng(n+1)=5.408620e+003;
n=4; farx(n+1)=1.652027e+002; foe(n+1)=7.408698e+003; krok(n+1)=1.457334e-003; ng(n+1)=2.685830e+003;
n=5; farx(n+1)=1.370136e+002; foe(n+1)=7.425340e+003; krok(n+1)=1.146884e-003; ng(n+1)=1.866528e+003;
n=6; farx(n+1)=1.156461e+002; foe(n+1)=7.214075e+003; krok(n+1)=3.941672e-003; ng(n+1)=1.441808e+003;
n=7; farx(n+1)=8.908153e+001; foe(n+1)=3.671224e+003; krok(n+1)=1.258748e-002; ng(n+1)=1.237082e+003;
n=8; farx(n+1)=6.413435e+001; foe(n+1)=1.424763e+003; krok(n+1)=3.571281e-002; ng(n+1)=5.811825e+002;
n=9; farx(n+1)=5.205310e+001; foe(n+1)=9.630666e+002; krok(n+1)=2.243330e-002; ng(n+1)=2.699820e+002;
n=10; farx(n+1)=4.754915e+001; foe(n+1)=1.041879e+003; krok(n+1)=1.825682e-002; ng(n+1)=2.738673e+002;
n=11; farx(n+1)=4.383070e+001; foe(n+1)=9.283238e+002; krok(n+1)=3.320390e-002; ng(n+1)=2.355348e+002;
n=12; farx(n+1)=3.530537e+001; foe(n+1)=8.119266e+002; krok(n+1)=1.341318e-001; ng(n+1)=3.220794e+002;
n=13; farx(n+1)=3.049613e+001; foe(n+1)=2.060334e+003; krok(n+1)=5.753341e-002; ng(n+1)=2.191234e+002;
n=14; farx(n+1)=2.839498e+001; foe(n+1)=4.409966e+002; krok(n+1)=5.483245e-002; ng(n+1)=2.201435e+002;
n=15; farx(n+1)=2.608732e+001; foe(n+1)=8.576411e+002; krok(n+1)=5.625683e-002; ng(n+1)=2.118625e+002;
n=16; farx(n+1)=2.425618e+001; foe(n+1)=4.263422e+002; krok(n+1)=1.238533e-001; ng(n+1)=3.030630e+002;
n=17; farx(n+1)=1.820068e+001; foe(n+1)=8.479722e+002; krok(n+1)=3.730776e-001; ng(n+1)=2.267026e+002;
n=18; farx(n+1)=1.605383e+001; foe(n+1)=1.401103e+003; krok(n+1)=5.403454e-001; ng(n+1)=3.209323e+002;
n=19; farx(n+1)=1.272809e+001; foe(n+1)=1.244027e+003; krok(n+1)=1.033731e+000; ng(n+1)=1.730714e+002;
n=20; farx(n+1)=9.793506e+000; foe(n+1)=2.623621e+003; krok(n+1)=4.989677e-001; ng(n+1)=1.336306e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
