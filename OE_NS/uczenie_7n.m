%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.546551e+003; foe(n+1)=4.681583e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=3.879108e+003; foe(n+1)=3.947978e+003; krok(n+1)=5.276810e-004; ng(n+1)=3.300383e+003;
n=2; farx(n+1)=1.627745e+003; foe(n+1)=1.568873e+003; krok(n+1)=4.150488e-003; ng(n+1)=1.631599e+003;
n=3; farx(n+1)=7.021693e+002; foe(n+1)=9.399296e+002; krok(n+1)=1.733806e-004; ng(n+1)=7.922266e+003;
n=4; farx(n+1)=7.486912e+002; foe(n+1)=8.861531e+002; krok(n+1)=3.569936e-004; ng(n+1)=1.788878e+003;
n=5; farx(n+1)=6.131749e+002; foe(n+1)=8.455855e+002; krok(n+1)=1.228120e-004; ng(n+1)=2.889011e+003;
n=6; farx(n+1)=6.377637e+002; foe(n+1)=8.175256e+002; krok(n+1)=2.563847e-004; ng(n+1)=1.266400e+003;
n=7; farx(n+1)=5.473917e+002; foe(n+1)=7.913192e+002; krok(n+1)=1.123699e-004; ng(n+1)=2.149044e+003;
n=8; farx(n+1)=5.675692e+002; foe(n+1)=7.681127e+002; krok(n+1)=2.463545e-004; ng(n+1)=1.304757e+003;
n=9; farx(n+1)=4.892224e+002; foe(n+1)=7.445417e+002; krok(n+1)=1.032442e-004; ng(n+1)=2.145065e+003;
n=10; farx(n+1)=5.017988e+002; foe(n+1)=7.223373e+002; krok(n+1)=2.439194e-004; ng(n+1)=1.228184e+003;
n=11; farx(n+1)=4.315199e+002; foe(n+1)=6.989206e+002; krok(n+1)=9.476985e-005; ng(n+1)=2.209511e+003;
n=12; farx(n+1)=4.395505e+002; foe(n+1)=6.771323e+002; krok(n+1)=2.266194e-004; ng(n+1)=1.253899e+003;
n=13; farx(n+1)=3.783184e+002; foe(n+1)=6.541159e+002; krok(n+1)=8.683747e-005; ng(n+1)=2.269676e+003;
n=14; farx(n+1)=3.827701e+002; foe(n+1)=6.318051e+002; krok(n+1)=2.215284e-004; ng(n+1)=1.277080e+003;
n=15; farx(n+1)=3.279985e+002; foe(n+1)=6.080584e+002; krok(n+1)=7.868993e-005; ng(n+1)=2.422105e+003;
n=16; farx(n+1)=3.302294e+002; foe(n+1)=5.858117e+002; krok(n+1)=2.044712e-004; ng(n+1)=1.318274e+003;
n=17; farx(n+1)=2.830648e+002; foe(n+1)=5.625689e+002; krok(n+1)=7.166629e-005; ng(n+1)=2.508928e+003;
n=18; farx(n+1)=2.843371e+002; foe(n+1)=5.409880e+002; krok(n+1)=1.858107e-004; ng(n+1)=1.371581e+003;
n=19; farx(n+1)=2.442732e+002; foe(n+1)=5.185485e+002; krok(n+1)=6.446963e-005; ng(n+1)=2.589367e+003;
n=20; farx(n+1)=2.442757e+002; foe(n+1)=4.971453e+002; krok(n+1)=1.764835e-004; ng(n+1)=1.409904e+003;
n=21; farx(n+1)=2.098132e+002; foe(n+1)=4.751663e+002; krok(n+1)=5.755042e-005; ng(n+1)=2.705910e+003;
n=22; farx(n+1)=2.091877e+002; foe(n+1)=4.544684e+002; krok(n+1)=1.635278e-004; ng(n+1)=1.451956e+003;
n=23; farx(n+1)=1.798796e+002; foe(n+1)=4.335227e+002; krok(n+1)=5.161420e-005; ng(n+1)=2.790383e+003;
n=24; farx(n+1)=1.793397e+002; foe(n+1)=4.142708e+002; krok(n+1)=1.460087e-004; ng(n+1)=1.505798e+003;
n=25; farx(n+1)=1.547612e+002; foe(n+1)=3.949165e+002; krok(n+1)=4.645269e-005; ng(n+1)=2.842542e+003;
n=26; farx(n+1)=1.543315e+002; foe(n+1)=3.771323e+002; krok(n+1)=1.304810e-004; ng(n+1)=1.551750e+003;
n=27; farx(n+1)=1.337511e+002; foe(n+1)=3.594230e+002; krok(n+1)=4.176829e-005; ng(n+1)=2.879559e+003;
n=28; farx(n+1)=1.332837e+002; foe(n+1)=3.428460e+002; krok(n+1)=1.198705e-004; ng(n+1)=1.590950e+003;
n=29; farx(n+1)=1.157999e+002; foe(n+1)=3.264621e+002; krok(n+1)=3.754715e-005; ng(n+1)=2.935103e+003;
n=30; farx(n+1)=1.153924e+002; foe(n+1)=3.112276e+002; krok(n+1)=1.089132e-004; ng(n+1)=1.624550e+003;
n=31; farx(n+1)=1.005999e+002; foe(n+1)=2.963617e+002; krok(n+1)=3.400479e-005; ng(n+1)=2.955854e+003;
n=32; farx(n+1)=1.003630e+002; foe(n+1)=2.825990e+002; krok(n+1)=9.866949e-005; ng(n+1)=1.648737e+003;
n=33; farx(n+1)=8.784253e+001; foe(n+1)=2.691871e+002; krok(n+1)=3.079431e-005; ng(n+1)=2.964090e+003;
n=34; farx(n+1)=8.760273e+001; foe(n+1)=2.565706e+002; krok(n+1)=9.229511e-005; ng(n+1)=1.653402e+003;
n=35; farx(n+1)=7.684952e+001; foe(n+1)=2.443522e+002; krok(n+1)=2.797908e-005; ng(n+1)=2.985330e+003;
n=36; farx(n+1)=7.671330e+001; foe(n+1)=2.330467e+002; krok(n+1)=8.489916e-005; ng(n+1)=1.660178e+003;
n=37; farx(n+1)=6.753734e+001; foe(n+1)=2.221428e+002; krok(n+1)=2.560827e-005; ng(n+1)=2.974535e+003;
n=38; farx(n+1)=6.748347e+001; foe(n+1)=2.120816e+002; krok(n+1)=7.773395e-005; ng(n+1)=1.654912e+003;
n=39; farx(n+1)=5.968824e+001; foe(n+1)=2.025250e+002; krok(n+1)=2.362476e-005; ng(n+1)=2.931488e+003;
n=40; farx(n+1)=5.971914e+001; foe(n+1)=1.936169e+002; krok(n+1)=7.225059e-005; ng(n+1)=1.638292e+003;
n=41; farx(n+1)=5.304847e+001; foe(n+1)=1.851549e+002; krok(n+1)=2.170937e-005; ng(n+1)=2.899152e+003;
n=42; farx(n+1)=5.303757e+001; foe(n+1)=1.770615e+002; krok(n+1)=7.022416e-005; ng(n+1)=1.600925e+003;
n=43; farx(n+1)=4.722336e+001; foe(n+1)=1.694090e+002; krok(n+1)=1.995189e-005; ng(n+1)=2.898034e+003;
n=44; farx(n+1)=4.721939e+001; foe(n+1)=1.621796e+002; krok(n+1)=6.738910e-005; ng(n+1)=1.558438e+003;
n=45; farx(n+1)=4.217478e+001; foe(n+1)=1.553774e+002; krok(n+1)=1.850362e-005; ng(n+1)=2.858890e+003;
n=46; farx(n+1)=4.222627e+001; foe(n+1)=1.490574e+002; krok(n+1)=6.334606e-005; ng(n+1)=1.516256e+003;
n=47; farx(n+1)=3.789393e+001; foe(n+1)=1.431375e+002; krok(n+1)=1.726474e-005; ng(n+1)=2.775765e+003;
n=48; farx(n+1)=3.796529e+001; foe(n+1)=1.375613e+002; krok(n+1)=6.078783e-005; ng(n+1)=1.464737e+003;
n=49; farx(n+1)=3.419094e+001; foe(n+1)=1.323498e+002; krok(n+1)=1.621284e-005; ng(n+1)=2.704741e+003;
n=50; farx(n+1)=3.432223e+001; foe(n+1)=1.275597e+002; krok(n+1)=5.618497e-005; ng(n+1)=1.416882e+003;
n=51; farx(n+1)=3.107270e+001; foe(n+1)=1.231309e+002; krok(n+1)=1.549259e-005; ng(n+1)=2.572417e+003;
n=52; farx(n+1)=3.129876e+001; foe(n+1)=1.191168e+002; krok(n+1)=5.097701e-005; ng(n+1)=1.377284e+003;
n=53; farx(n+1)=2.850263e+001; foe(n+1)=1.153526e+002; krok(n+1)=1.481560e-005; ng(n+1)=2.435044e+003;
n=54; farx(n+1)=2.873724e+001; foe(n+1)=1.118472e+002; krok(n+1)=4.789661e-005; ng(n+1)=1.328110e+003;
n=55; farx(n+1)=2.631353e+001; foe(n+1)=1.086106e+002; krok(n+1)=1.416371e-005; ng(n+1)=2.312538e+003;
n=56; farx(n+1)=2.654257e+001; foe(n+1)=1.054836e+002; krok(n+1)=4.738492e-005; ng(n+1)=1.274864e+003;
n=57; farx(n+1)=2.436669e+001; foe(n+1)=1.025503e+002; krok(n+1)=1.346311e-005; ng(n+1)=2.260770e+003;
n=58; farx(n+1)=2.459088e+001; foe(n+1)=9.978003e+001; krok(n+1)=4.599898e-005; ng(n+1)=1.222285e+003;
n=59; farx(n+1)=2.264456e+001; foe(n+1)=9.718304e+001; krok(n+1)=1.302585e-005; ng(n+1)=2.175135e+003;
n=60; farx(n+1)=2.290629e+001; foe(n+1)=9.483924e+001; krok(n+1)=4.204485e-005; ng(n+1)=1.180526e+003;
n=61; farx(n+1)=2.121643e+001; foe(n+1)=9.262106e+001; krok(n+1)=1.259797e-005; ng(n+1)=2.041744e+003;
n=62; farx(n+1)=2.145292e+001; foe(n+1)=9.048811e+001; krok(n+1)=4.194015e-005; ng(n+1)=1.127328e+003;
n=63; farx(n+1)=1.992307e+001; foe(n+1)=8.848211e+001; krok(n+1)=1.216245e-005; ng(n+1)=1.978966e+003;
n=64; farx(n+1)=2.015889e+001; foe(n+1)=8.659945e+001; krok(n+1)=4.014885e-005; ng(n+1)=1.081518e+003;
n=65; farx(n+1)=1.879171e+001; foe(n+1)=8.483242e+001; krok(n+1)=1.186082e-005; ng(n+1)=1.883942e+003;
n=66; farx(n+1)=1.903068e+001; foe(n+1)=8.317081e+001; krok(n+1)=3.858161e-005; ng(n+1)=1.039382e+003;
n=67; farx(n+1)=1.779864e+001; foe(n+1)=8.159570e+001; krok(n+1)=1.157184e-005; ng(n+1)=1.802761e+003;
n=68; farx(n+1)=1.803437e+001; foe(n+1)=8.011479e+001; krok(n+1)=3.723990e-005; ng(n+1)=9.992815e+002;
n=69; farx(n+1)=1.692144e+001; foe(n+1)=7.870734e+001; krok(n+1)=1.131276e-005; ng(n+1)=1.725932e+003;
n=70; farx(n+1)=1.715106e+001; foe(n+1)=7.737575e+001; krok(n+1)=3.621888e-005; ng(n+1)=9.612937e+002;
n=71; farx(n+1)=1.613778e+001; foe(n+1)=7.610841e+001; krok(n+1)=1.109569e-005; ng(n+1)=1.657235e+003;
n=72; farx(n+1)=1.636414e+001; foe(n+1)=7.491589e+001; krok(n+1)=3.488004e-005; ng(n+1)=9.260578e+002;
n=73; farx(n+1)=1.544755e+001; foe(n+1)=7.377879e+001; krok(n+1)=1.085468e-005; ng(n+1)=1.584507e+003;
n=74; farx(n+1)=1.566104e+001; foe(n+1)=7.267390e+001; krok(n+1)=3.514435e-005; ng(n+1)=8.900918e+002;
n=75; farx(n+1)=1.481266e+001; foe(n+1)=7.161815e+001; krok(n+1)=1.057446e-005; ng(n+1)=1.546434e+003;
n=76; farx(n+1)=1.501581e+001; foe(n+1)=7.059727e+001; krok(n+1)=3.517612e-005; ng(n+1)=8.556095e+002;
n=77; farx(n+1)=1.422790e+001; foe(n+1)=6.962082e+001; krok(n+1)=1.038706e-005; ng(n+1)=1.504452e+003;
n=78; farx(n+1)=1.442876e+001; foe(n+1)=6.869982e+001; krok(n+1)=3.400479e-005; ng(n+1)=8.244965e+002;
n=79; farx(n+1)=1.370405e+001; foe(n+1)=6.781532e+001; krok(n+1)=1.029142e-005; ng(n+1)=1.442322e+003;
n=80; farx(n+1)=1.390337e+001; foe(n+1)=6.698762e+001; krok(n+1)=3.245086e-005; ng(n+1)=7.960114e+002;
n=81; farx(n+1)=1.324176e+001; foe(n+1)=6.619162e+001; krok(n+1)=1.018440e-005; ng(n+1)=1.374432e+003;
n=82; farx(n+1)=1.343453e+001; foe(n+1)=6.542549e+001; krok(n+1)=3.217628e-005; ng(n+1)=7.676470e+002;
n=83; farx(n+1)=1.281497e+001; foe(n+1)=6.468717e+001; krok(n+1)=1.009346e-005; ng(n+1)=1.334100e+003;
n=84; farx(n+1)=1.300538e+001; foe(n+1)=6.399501e+001; krok(n+1)=3.061523e-005; ng(n+1)=7.427351e+002;
n=85; farx(n+1)=1.243400e+001; foe(n+1)=6.332853e+001; krok(n+1)=1.009346e-005; ng(n+1)=1.270149e+003;
n=86; farx(n+1)=1.262322e+001; foe(n+1)=6.270067e+001; krok(n+1)=2.932016e-005; ng(n+1)=7.198587e+002;
n=87; farx(n+1)=1.209388e+001; foe(n+1)=6.209110e+001; krok(n+1)=1.003408e-005; ng(n+1)=1.217783e+003;
n=88; farx(n+1)=1.227675e+001; foe(n+1)=6.150694e+001; krok(n+1)=2.877521e-005; ng(n+1)=6.969021e+002;
n=89; farx(n+1)=1.178210e+001; foe(n+1)=6.094222e+001; krok(n+1)=9.955536e-006; ng(n+1)=1.176507e+003;
n=90; farx(n+1)=1.195877e+001; foe(n+1)=6.039654e+001; krok(n+1)=2.832743e-005; ng(n+1)=6.750447e+002;
n=91; farx(n+1)=1.149193e+001; foe(n+1)=5.987135e+001; krok(n+1)=9.952471e-006; ng(n+1)=1.138637e+003;
n=92; farx(n+1)=1.166795e+001; foe(n+1)=5.937368e+001; krok(n+1)=2.713913e-005; ng(n+1)=6.564379e+002;
n=93; farx(n+1)=1.123019e+001; foe(n+1)=5.888847e+001; krok(n+1)=9.952471e-006; ng(n+1)=1.095585e+003;
n=94; farx(n+1)=1.140135e+001; foe(n+1)=5.842728e+001; krok(n+1)=2.617217e-005; ng(n+1)=6.381740e+002;
n=95; farx(n+1)=1.099442e+001; foe(n+1)=5.797985e+001; krok(n+1)=9.888364e-006; ng(n+1)=1.051387e+003;
n=96; farx(n+1)=1.116005e+001; foe(n+1)=5.753842e+001; krok(n+1)=2.668806e-005; ng(n+1)=6.186210e+002;
n=97; farx(n+1)=1.077083e+001; foe(n+1)=5.710847e+001; krok(n+1)=9.716744e-006; ng(n+1)=1.039019e+003;
n=98; farx(n+1)=1.093050e+001; foe(n+1)=5.668561e+001; krok(n+1)=2.711505e-005; ng(n+1)=5.996425e+002;
n=99; farx(n+1)=1.055796e+001; foe(n+1)=5.627430e+001; krok(n+1)=9.602384e-006; ng(n+1)=1.022905e+003;
n=100; farx(n+1)=1.071380e+001; foe(n+1)=5.587430e+001; krok(n+1)=2.713913e-005; ng(n+1)=5.822177e+002;
n=101; farx(n+1)=1.035822e+001; foe(n+1)=5.548302e+001; krok(n+1)=9.515478e-006; ng(n+1)=1.002564e+003;
n=102; farx(n+1)=1.050920e+001; foe(n+1)=5.510501e+001; krok(n+1)=2.689247e-005; ng(n+1)=5.655985e+002;
n=103; farx(n+1)=1.017169e+001; foe(n+1)=5.473671e+001; krok(n+1)=9.471564e-006; ng(n+1)=9.750089e+002;
n=104; farx(n+1)=1.031860e+001; foe(n+1)=5.438031e+001; krok(n+1)=2.657240e-005; ng(n+1)=5.499619e+002;
n=105; farx(n+1)=9.998579e+000; foe(n+1)=5.403337e+001; krok(n+1)=9.410994e-006; ng(n+1)=9.482785e+002;
n=106; farx(n+1)=1.014180e+001; foe(n+1)=5.369401e+001; krok(n+1)=2.671200e-005; ng(n+1)=5.344428e+002;
n=107; farx(n+1)=9.834932e+000; foe(n+1)=5.336250e+001; krok(n+1)=9.331364e-006; ng(n+1)=9.311349e+002;
n=108; farx(n+1)=9.974071e+000; foe(n+1)=5.304052e+001; krok(n+1)=2.657240e-005; ng(n+1)=5.197227e+002;
n=109; farx(n+1)=9.681369e+000; foe(n+1)=5.272663e+001; krok(n+1)=9.282697e-006; ng(n+1)=9.082794e+002;
n=110; farx(n+1)=9.817241e+000; foe(n+1)=5.242128e+001; krok(n+1)=2.648762e-005; ng(n+1)=5.056682e+002;
n=111; farx(n+1)=9.537082e+000; foe(n+1)=5.212269e+001; krok(n+1)=9.225986e-006; ng(n+1)=8.886788e+002;
n=112; farx(n+1)=9.669379e+000; foe(n+1)=5.183265e+001; krok(n+1)=2.637524e-005; ng(n+1)=4.920844e+002;
n=113; farx(n+1)=9.401437e+000; foe(n+1)=5.154912e+001; krok(n+1)=9.177209e-006; ng(n+1)=8.683281e+002;
n=114; farx(n+1)=9.530239e+000; foe(n+1)=5.127329e+001; krok(n+1)=2.627836e-005; ng(n+1)=4.789795e+002;
n=115; farx(n+1)=9.273543e+000; foe(n+1)=5.100386e+001; krok(n+1)=9.137433e-006; ng(n+1)=8.487843e+002;
n=116; farx(n+1)=9.399521e+000; foe(n+1)=5.074207e+001; krok(n+1)=2.616958e-005; ng(n+1)=4.665072e+002;
n=117; farx(n+1)=9.153205e+000; foe(n+1)=5.048535e+001; krok(n+1)=9.089066e-006; ng(n+1)=8.309755e+002;
n=118; farx(n+1)=9.274467e+000; foe(n+1)=5.023607e+001; krok(n+1)=2.581106e-005; ng(n+1)=4.545717e+002;
n=119; farx(n+1)=9.039687e+000; foe(n+1)=4.999507e+001; krok(n+1)=9.109987e-006; ng(n+1)=8.047393e+002;
n=120; farx(n+1)=9.159263e+000; foe(n+1)=4.976029e+001; krok(n+1)=2.549454e-005; ng(n+1)=4.441399e+002;
n=121; farx(n+1)=8.933545e+000; foe(n+1)=4.953064e+001; krok(n+1)=9.078452e-006; ng(n+1)=7.874192e+002;
n=122; farx(n+1)=9.050151e+000; foe(n+1)=4.930739e+001; krok(n+1)=2.529606e-005; ng(n+1)=4.336639e+002;
n=123; farx(n+1)=8.833536e+000; foe(n+1)=4.908906e+001; krok(n+1)=9.050777e-006; ng(n+1)=7.689531e+002;
n=124; farx(n+1)=8.947878e+000; foe(n+1)=4.887616e+001; krok(n+1)=2.531441e-005; ng(n+1)=4.234314e+002;
n=125; farx(n+1)=8.738712e+000; foe(n+1)=4.866667e+001; krok(n+1)=8.997287e-006; ng(n+1)=7.557697e+002;
n=126; farx(n+1)=8.850232e+000; foe(n+1)=4.846302e+001; krok(n+1)=2.531441e-005; ng(n+1)=4.133624e+002;
n=127; farx(n+1)=8.648559e+000; foe(n+1)=4.826266e+001; krok(n+1)=8.958287e-006; ng(n+1)=7.410928e+002;
n=128; farx(n+1)=8.756609e+000; foe(n+1)=4.806809e+001; krok(n+1)=2.509886e-005; ng(n+1)=4.036539e+002;
n=129; farx(n+1)=8.563340e+000; foe(n+1)=4.787819e+001; krok(n+1)=8.943311e-006; ng(n+1)=7.217944e+002;
n=130; farx(n+1)=8.669385e+000; foe(n+1)=4.769210e+001; krok(n+1)=2.518354e-005; ng(n+1)=3.942519e+002;
n=131; farx(n+1)=8.481683e+000; foe(n+1)=4.750929e+001; krok(n+1)=8.913402e-006; ng(n+1)=7.110606e+002;
n=132; farx(n+1)=8.584533e+000; foe(n+1)=4.733269e+001; krok(n+1)=2.468559e-005; ng(n+1)=3.858486e+002;
n=133; farx(n+1)=8.404966e+000; foe(n+1)=4.716043e+001; krok(n+1)=8.935769e-006; ng(n+1)=6.892289e+002;
n=134; farx(n+1)=8.505910e+000; foe(n+1)=4.699245e+001; krok(n+1)=2.450546e-005; ng(n+1)=3.779409e+002;
n=135; farx(n+1)=8.332805e+000; foe(n+1)=4.682756e+001; krok(n+1)=8.887598e-006; ng(n+1)=6.757347e+002;
n=136; farx(n+1)=8.431609e+000; foe(n+1)=4.666538e+001; krok(n+1)=2.480406e-005; ng(n+1)=3.697216e+002;
n=137; farx(n+1)=8.263326e+000; foe(n+1)=4.650581e+001; krok(n+1)=8.818566e-006; ng(n+1)=6.678365e+002;
n=138; farx(n+1)=8.359669e+000; foe(n+1)=4.634943e+001; krok(n+1)=2.496301e-005; ng(n+1)=3.617089e+002;
n=139; farx(n+1)=8.196349e+000; foe(n+1)=4.619589e+001; krok(n+1)=8.786087e-006; ng(n+1)=6.571017e+002;
n=140; farx(n+1)=8.290732e+000; foe(n+1)=4.604637e+001; krok(n+1)=2.490114e-005; ng(n+1)=3.541807e+002;
n=141; farx(n+1)=8.132598e+000; foe(n+1)=4.589899e+001; krok(n+1)=8.753864e-006; ng(n+1)=6.452512e+002;
n=142; farx(n+1)=8.223949e+000; foe(n+1)=4.575548e+001; krok(n+1)=2.468559e-005; ng(n+1)=3.467561e+002;
n=143; farx(n+1)=8.071863e+000; foe(n+1)=4.561555e+001; krok(n+1)=8.756691e-006; ng(n+1)=6.284414e+002;
n=144; farx(n+1)=8.162067e+000; foe(n+1)=4.547814e+001; krok(n+1)=2.480406e-005; ng(n+1)=3.399277e+002;
n=145; farx(n+1)=8.013457e+000; foe(n+1)=4.534244e+001; krok(n+1)=8.718949e-006; ng(n+1)=6.215662e+002;
n=146; farx(n+1)=8.101011e+000; foe(n+1)=4.521118e+001; krok(n+1)=2.443486e-005; ng(n+1)=3.333862e+002;
n=147; farx(n+1)=7.958308e+000; foe(n+1)=4.508251e+001; krok(n+1)=8.720009e-006; ng(n+1)=6.045359e+002;
n=148; farx(n+1)=8.044061e+000; foe(n+1)=4.495618e+001; krok(n+1)=2.450546e-005; ng(n+1)=3.268111e+002;
n=149; farx(n+1)=7.905297e+000; foe(n+1)=4.483208e+001; krok(n+1)=8.686387e-006; ng(n+1)=5.954068e+002;
n=150; farx(n+1)=7.988937e+000; foe(n+1)=4.471060e+001; krok(n+1)=2.443486e-005; ng(n+1)=3.208682e+002;
n=151; farx(n+1)=7.854543e+000; foe(n+1)=4.459158e+001; krok(n+1)=8.668026e-006; ng(n+1)=5.837713e+002;
n=152; farx(n+1)=7.936382e+000; foe(n+1)=4.447469e+001; krok(n+1)=2.443486e-005; ng(n+1)=3.152119e+002;
n=153; farx(n+1)=7.805960e+000; foe(n+1)=4.436009e+001; krok(n+1)=8.632369e-006; ng(n+1)=5.740319e+002;
n=154; farx(n+1)=7.886063e+000; foe(n+1)=4.424709e+001; krok(n+1)=2.458493e-005; ng(n+1)=3.094374e+002;
n=155; farx(n+1)=7.759039e+000; foe(n+1)=4.413628e+001; krok(n+1)=8.593476e-006; ng(n+1)=5.661352e+002;
n=156; farx(n+1)=7.837934e+000; foe(n+1)=4.402723e+001; krok(n+1)=2.480406e-005; ng(n+1)=3.040833e+002;
n=157; farx(n+1)=7.714269e+000; foe(n+1)=4.391942e+001; krok(n+1)=8.500298e-006; ng(n+1)=5.607271e+002;
n=158; farx(n+1)=7.791473e+000; foe(n+1)=4.381206e+001; krok(n+1)=2.556514e-005; ng(n+1)=2.978665e+002;
n=159; farx(n+1)=7.669607e+000; foe(n+1)=4.370627e+001; krok(n+1)=8.445281e-006; ng(n+1)=5.587901e+002;
n=160; farx(n+1)=7.745063e+000; foe(n+1)=4.360322e+001; krok(n+1)=2.541507e-005; ng(n+1)=2.930108e+002;
n=161; farx(n+1)=7.626959e+000; foe(n+1)=4.350172e+001; krok(n+1)=8.429495e-006; ng(n+1)=5.478287e+002;
n=162; farx(n+1)=7.700598e+000; foe(n+1)=4.340231e+001; krok(n+1)=2.537849e-005; ng(n+1)=2.880770e+002;
n=163; farx(n+1)=7.586370e+000; foe(n+1)=4.330471e+001; krok(n+1)=8.379583e-006; ng(n+1)=5.378056e+002;
n=164; farx(n+1)=7.658048e+000; foe(n+1)=4.320725e+001; krok(n+1)=2.581106e-005; ng(n+1)=2.825348e+002;
n=165; farx(n+1)=7.545998e+000; foe(n+1)=4.311253e+001; krok(n+1)=8.378583e-006; ng(n+1)=5.314696e+002;
n=166; farx(n+1)=7.616745e+000; foe(n+1)=4.302008e+001; krok(n+1)=2.543436e-005; ng(n+1)=2.786220e+002;
n=167; farx(n+1)=7.507907e+000; foe(n+1)=4.292895e+001; krok(n+1)=8.364627e-006; ng(n+1)=5.217959e+002;
n=168; farx(n+1)=7.577155e+000; foe(n+1)=4.283972e+001; krok(n+1)=2.541507e-005; ng(n+1)=2.740971e+002;
n=169; farx(n+1)=7.471364e+000; foe(n+1)=4.275178e+001; krok(n+1)=8.333245e-006; ng(n+1)=5.133252e+002;
n=170; farx(n+1)=7.538999e+000; foe(n+1)=4.266506e+001; krok(n+1)=2.555890e-005; ng(n+1)=2.693938e+002;
n=171; farx(n+1)=7.435616e+000; foe(n+1)=4.257999e+001; krok(n+1)=8.329312e-006; ng(n+1)=5.058599e+002;
n=172; farx(n+1)=7.501735e+000; foe(n+1)=4.249696e+001; krok(n+1)=2.519594e-005; ng(n+1)=2.654829e+002;
n=173; farx(n+1)=7.401350e+000; foe(n+1)=4.241553e+001; krok(n+1)=8.357584e-006; ng(n+1)=4.946579e+002;
n=174; farx(n+1)=7.465959e+000; foe(n+1)=4.233630e+001; krok(n+1)=2.468559e-005; ng(n+1)=2.618936e+002;
n=175; farx(n+1)=7.368396e+000; foe(n+1)=4.225860e+001; krok(n+1)=8.423638e-006; ng(n+1)=4.824637e+002;
n=176; farx(n+1)=7.431363e+000; foe(n+1)=4.218382e+001; krok(n+1)=2.371798e-005; ng(n+1)=2.590269e+002;
n=177; farx(n+1)=7.337771e+000; foe(n+1)=4.211054e+001; krok(n+1)=8.481911e-006; ng(n+1)=4.662380e+002;
n=178; farx(n+1)=7.399848e+000; foe(n+1)=4.203855e+001; krok(n+1)=2.360381e-005; ng(n+1)=2.554893e+002;
n=179; farx(n+1)=7.308553e+000; foe(n+1)=4.196745e+001; krok(n+1)=8.445859e-006; ng(n+1)=4.599968e+002;
n=180; farx(n+1)=7.369641e+000; foe(n+1)=4.189713e+001; krok(n+1)=2.390184e-005; ng(n+1)=2.515097e+002;
n=181; farx(n+1)=7.279689e+000; foe(n+1)=4.182754e+001; krok(n+1)=8.423638e-006; ng(n+1)=4.567660e+002;
n=182; farx(n+1)=7.339442e+000; foe(n+1)=4.175979e+001; krok(n+1)=2.360381e-005; ng(n+1)=2.481718e+002;
n=183; farx(n+1)=7.252581e+000; foe(n+1)=4.169288e+001; krok(n+1)=8.379583e-006; ng(n+1)=4.475271e+002;
n=184; farx(n+1)=7.311626e+000; foe(n+1)=4.162570e+001; krok(n+1)=2.444110e-005; ng(n+1)=2.439853e+002;
n=185; farx(n+1)=7.225150e+000; foe(n+1)=4.155916e+001; krok(n+1)=8.316918e-006; ng(n+1)=4.494643e+002;
n=186; farx(n+1)=7.282743e+000; foe(n+1)=4.149416e+001; krok(n+1)=2.419290e-005; ng(n+1)=2.406677e+002;
n=187; farx(n+1)=7.198808e+000; foe(n+1)=4.143021e+001; krok(n+1)=8.322079e-006; ng(n+1)=4.402414e+002;
n=188; farx(n+1)=7.255349e+000; foe(n+1)=4.136708e+001; krok(n+1)=2.419290e-005; ng(n+1)=2.373617e+002;
n=189; farx(n+1)=7.173138e+000; foe(n+1)=4.130496e+001; krok(n+1)=8.318815e-006; ng(n+1)=4.344977e+002;
n=190; farx(n+1)=7.229114e+000; foe(n+1)=4.124397e+001; krok(n+1)=2.418894e-005; ng(n+1)=2.342748e+002;
n=191; farx(n+1)=7.148409e+000; foe(n+1)=4.118324e+001; krok(n+1)=8.284610e-006; ng(n+1)=4.307185e+002;
n=192; farx(n+1)=7.202580e+000; foe(n+1)=4.112378e+001; krok(n+1)=2.394830e-005; ng(n+1)=2.310394e+002;
n=193; farx(n+1)=7.124230e+000; foe(n+1)=4.106567e+001; krok(n+1)=8.334046e-006; ng(n+1)=4.203595e+002;
n=194; farx(n+1)=7.177715e+000; foe(n+1)=4.100861e+001; krok(n+1)=2.362476e-005; ng(n+1)=2.283431e+002;
n=195; farx(n+1)=7.101269e+000; foe(n+1)=4.095229e+001; krok(n+1)=8.322079e-006; ng(n+1)=4.139908e+002;
n=196; farx(n+1)=7.154120e+000; foe(n+1)=4.089659e+001; krok(n+1)=2.382783e-005; ng(n+1)=2.256232e+002;
n=197; farx(n+1)=7.078895e+000; foe(n+1)=4.084130e+001; krok(n+1)=8.274713e-006; ng(n+1)=4.115900e+002;
n=198; farx(n+1)=7.131222e+000; foe(n+1)=4.078672e+001; krok(n+1)=2.418894e-005; ng(n+1)=2.228256e+002;
n=199; farx(n+1)=7.056911e+000; foe(n+1)=4.073214e+001; krok(n+1)=8.219763e-006; ng(n+1)=4.106982e+002;
n=200; farx(n+1)=7.107418e+000; foe(n+1)=4.067870e+001; krok(n+1)=2.394830e-005; ng(n+1)=2.200916e+002;
n=201; farx(n+1)=7.035247e+000; foe(n+1)=4.062652e+001; krok(n+1)=8.284610e-006; ng(n+1)=4.004277e+002;
n=202; farx(n+1)=7.085032e+000; foe(n+1)=4.057537e+001; krok(n+1)=2.343658e-005; ng(n+1)=2.182843e+002;
n=203; farx(n+1)=7.014595e+000; foe(n+1)=4.052501e+001; krok(n+1)=8.316918e-006; ng(n+1)=3.929771e+002;
n=204; farx(n+1)=7.063832e+000; foe(n+1)=4.047565e+001; krok(n+1)=2.322634e-005; ng(n+1)=2.162810e+002;
n=205; farx(n+1)=6.994708e+000; foe(n+1)=4.042657e+001; krok(n+1)=8.316918e-006; ng(n+1)=3.883295e+002;
n=206; farx(n+1)=7.042680e+000; foe(n+1)=4.037864e+001; krok(n+1)=2.286087e-005; ng(n+1)=2.143808e+002;
n=207; farx(n+1)=6.975667e+000; foe(n+1)=4.033148e+001; krok(n+1)=8.333245e-006; ng(n+1)=3.797398e+002;
n=208; farx(n+1)=7.023183e+000; foe(n+1)=4.028456e+001; krok(n+1)=2.307378e-005; ng(n+1)=2.122012e+002;
n=209; farx(n+1)=6.956813e+000; foe(n+1)=4.023807e+001; krok(n+1)=8.316918e-006; ng(n+1)=3.783625e+002;
n=210; farx(n+1)=7.003381e+000; foe(n+1)=4.019262e+001; krok(n+1)=2.273062e-005; ng(n+1)=2.105279e+002;
n=211; farx(n+1)=6.938659e+000; foe(n+1)=4.014771e+001; krok(n+1)=8.342260e-006; ng(n+1)=3.713466e+002;
n=212; farx(n+1)=6.984451e+000; foe(n+1)=4.010352e+001; krok(n+1)=2.254192e-005; ng(n+1)=2.087748e+002;
n=213; farx(n+1)=6.921204e+000; foe(n+1)=4.005986e+001; krok(n+1)=8.333245e-006; ng(n+1)=3.661061e+002;
n=214; farx(n+1)=6.966557e+000; foe(n+1)=4.001650e+001; krok(n+1)=2.277170e-005; ng(n+1)=2.066574e+002;
n=215; farx(n+1)=6.904078e+000; foe(n+1)=3.997342e+001; krok(n+1)=8.287718e-006; ng(n+1)=3.649910e+002;
n=216; farx(n+1)=6.948631e+000; foe(n+1)=3.993083e+001; krok(n+1)=2.285025e-005; ng(n+1)=2.045270e+002;
n=217; farx(n+1)=6.887260e+000; foe(n+1)=3.988871e+001; krok(n+1)=8.274713e-006; ng(n+1)=3.613613e+002;
n=218; farx(n+1)=6.931138e+000; foe(n+1)=3.984703e+001; krok(n+1)=2.285025e-005; ng(n+1)=2.026152e+002;
n=219; farx(n+1)=6.870773e+000; foe(n+1)=3.980583e+001; krok(n+1)=8.272636e-006; ng(n+1)=3.577950e+002;
n=220; farx(n+1)=6.913868e+000; foe(n+1)=3.976522e+001; krok(n+1)=2.265271e-005; ng(n+1)=2.009871e+002;
n=221; farx(n+1)=6.854758e+000; foe(n+1)=3.972518e+001; krok(n+1)=8.287718e-006; ng(n+1)=3.525869e+002;
n=222; farx(n+1)=6.896832e+000; foe(n+1)=3.968562e+001; krok(n+1)=2.230894e-005; ng(n+1)=1.994451e+002;
n=223; farx(n+1)=6.839193e+000; foe(n+1)=3.964698e+001; krok(n+1)=8.333245e-006; ng(n+1)=3.454994e+002;
n=224; farx(n+1)=6.880993e+000; foe(n+1)=3.960857e+001; krok(n+1)=2.219549e-005; ng(n+1)=1.981343e+002;
n=225; farx(n+1)=6.824153e+000; foe(n+1)=3.957062e+001; krok(n+1)=8.318815e-006; ng(n+1)=3.431074e+002;
n=226; farx(n+1)=6.865271e+000; foe(n+1)=3.953305e+001; krok(n+1)=2.209045e-005; ng(n+1)=1.965377e+002;
n=227; farx(n+1)=6.809493e+000; foe(n+1)=3.949601e+001; krok(n+1)=8.322079e-006; ng(n+1)=3.389989e+002;
n=228; farx(n+1)=6.850359e+000; foe(n+1)=3.945927e+001; krok(n+1)=2.219138e-005; ng(n+1)=1.950421e+002;
n=229; farx(n+1)=6.795256e+000; foe(n+1)=3.942268e+001; krok(n+1)=8.274713e-006; ng(n+1)=3.380370e+002;
n=230; farx(n+1)=6.835831e+000; foe(n+1)=3.938640e+001; krok(n+1)=2.253084e-005; ng(n+1)=1.932477e+002;
n=231; farx(n+1)=6.781267e+000; foe(n+1)=3.935008e+001; krok(n+1)=8.206521e-006; ng(n+1)=3.382793e+002;
n=232; farx(n+1)=6.821028e+000; foe(n+1)=3.931410e+001; krok(n+1)=2.273062e-005; ng(n+1)=1.913773e+002;
n=233; farx(n+1)=6.767395e+000; foe(n+1)=3.927849e+001; krok(n+1)=8.186754e-006; ng(n+1)=3.353177e+002;
n=234; farx(n+1)=6.806840e+000; foe(n+1)=3.924306e+001; krok(n+1)=2.296532e-005; ng(n+1)=1.899009e+002;
n=235; farx(n+1)=6.753688e+000; foe(n+1)=3.920782e+001; krok(n+1)=8.153596e-006; ng(n+1)=3.347877e+002;
n=236; farx(n+1)=6.792403e+000; foe(n+1)=3.917305e+001; krok(n+1)=2.286087e-005; ng(n+1)=1.885489e+002;
n=237; farx(n+1)=6.740323e+000; foe(n+1)=3.913869e+001; krok(n+1)=8.155065e-006; ng(n+1)=3.304955e+002;
n=238; farx(n+1)=6.778504e+000; foe(n+1)=3.910455e+001; krok(n+1)=2.286087e-005; ng(n+1)=1.872169e+002;
n=239; farx(n+1)=6.727122e+000; foe(n+1)=3.907085e+001; krok(n+1)=8.164371e-006; ng(n+1)=3.277411e+002;
n=240; farx(n+1)=6.764902e+000; foe(n+1)=3.903763e+001; krok(n+1)=2.266263e-005; ng(n+1)=1.861811e+002;
n=241; farx(n+1)=6.714361e+000; foe(n+1)=3.900462e+001; krok(n+1)=8.164371e-006; ng(n+1)=3.245472e+002;
n=242; farx(n+1)=6.751157e+000; foe(n+1)=3.897208e+001; krok(n+1)=2.230894e-005; ng(n+1)=1.849766e+002;
n=243; farx(n+1)=6.701839e+000; foe(n+1)=3.894024e+001; krok(n+1)=8.226685e-006; ng(n+1)=3.177633e+002;
n=244; farx(n+1)=6.738470e+000; foe(n+1)=3.890869e+001; krok(n+1)=2.208984e-005; ng(n+1)=1.841499e+002;
n=245; farx(n+1)=6.689873e+000; foe(n+1)=3.887738e+001; krok(n+1)=8.206521e-006; ng(n+1)=3.156244e+002;
n=246; farx(n+1)=6.726201e+000; foe(n+1)=3.884627e+001; krok(n+1)=2.226248e-005; ng(n+1)=1.828390e+002;
n=247; farx(n+1)=6.678056e+000; foe(n+1)=3.881529e+001; krok(n+1)=8.178765e-006; ng(n+1)=3.149687e+002;
n=248; farx(n+1)=6.713926e+000; foe(n+1)=3.878471e+001; krok(n+1)=2.226248e-005; ng(n+1)=1.816346e+002;
n=249; farx(n+1)=6.666550e+000; foe(n+1)=3.875427e+001; krok(n+1)=8.155065e-006; ng(n+1)=3.126397e+002;
n=250; farx(n+1)=6.701783e+000; foe(n+1)=3.872402e+001; krok(n+1)=2.230894e-005; ng(n+1)=1.802692e+002;
n=251; farx(n+1)=6.655101e+000; foe(n+1)=3.869418e+001; krok(n+1)=8.166804e-006; ng(n+1)=3.096814e+002;
n=252; farx(n+1)=6.690141e+000; foe(n+1)=3.866463e+001; krok(n+1)=2.226248e-005; ng(n+1)=1.793147e+002;
n=253; farx(n+1)=6.644074e+000; foe(n+1)=3.863517e+001; krok(n+1)=8.129247e-006; ng(n+1)=3.084040e+002;
n=254; farx(n+1)=6.678408e+000; foe(n+1)=3.860587e+001; krok(n+1)=2.230894e-005; ng(n+1)=1.779231e+002;
n=255; farx(n+1)=6.632971e+000; foe(n+1)=3.857704e+001; krok(n+1)=8.164371e-006; ng(n+1)=3.050894e+002;
n=256; farx(n+1)=6.667040e+000; foe(n+1)=3.854861e+001; krok(n+1)=2.198334e-005; ng(n+1)=1.772509e+002;
n=257; farx(n+1)=6.622409e+000; foe(n+1)=3.852037e+001; krok(n+1)=8.155065e-006; ng(n+1)=3.020472e+002;
n=258; farx(n+1)=6.656242e+000; foe(n+1)=3.849228e+001; krok(n+1)=2.219138e-005; ng(n+1)=1.760255e+002;
n=259; farx(n+1)=6.611923e+000; foe(n+1)=3.846426e+001; krok(n+1)=8.128882e-006; ng(n+1)=3.018904e+002;
n=260; farx(n+1)=6.645266e+000; foe(n+1)=3.843663e+001; krok(n+1)=2.208984e-005; ng(n+1)=1.749809e+002;
n=261; farx(n+1)=6.601568e+000; foe(n+1)=3.840918e+001; krok(n+1)=8.151193e-006; ng(n+1)=2.989038e+002;
n=262; farx(n+1)=6.634663e+000; foe(n+1)=3.838222e+001; krok(n+1)=2.190752e-005; ng(n+1)=1.742234e+002;
n=263; farx(n+1)=6.591720e+000; foe(n+1)=3.835521e+001; krok(n+1)=8.114575e-006; ng(n+1)=2.967287e+002;
n=264; farx(n+1)=6.624450e+000; foe(n+1)=3.832832e+001; krok(n+1)=2.226248e-005; ng(n+1)=1.728285e+002;
n=265; farx(n+1)=6.581703e+000; foe(n+1)=3.830152e+001; krok(n+1)=8.106422e-006; ng(n+1)=2.966332e+002;
n=266; farx(n+1)=6.613600e+000; foe(n+1)=3.827525e+001; krok(n+1)=2.170937e-005; ng(n+1)=1.720745e+002;
n=267; farx(n+1)=6.572020e+000; foe(n+1)=3.824945e+001; krok(n+1)=8.172419e-006; ng(n+1)=2.896936e+002;
n=268; farx(n+1)=6.603472e+000; foe(n+1)=3.822378e+001; krok(n+1)=2.139878e-005; ng(n+1)=1.713533e+002;
n=269; farx(n+1)=6.562599e+000; foe(n+1)=3.819858e+001; krok(n+1)=8.206521e-006; ng(n+1)=2.859085e+002;
n=270; farx(n+1)=6.593956e+000; foe(n+1)=3.817348e+001; krok(n+1)=2.132008e-005; ng(n+1)=1.706442e+002;
n=271; farx(n+1)=6.553379e+000; foe(n+1)=3.814854e+001; krok(n+1)=8.205260e-006; ng(n+1)=2.850596e+002;
n=272; farx(n+1)=6.584567e+000; foe(n+1)=3.812396e+001; krok(n+1)=2.122479e-005; ng(n+1)=1.699135e+002;
n=273; farx(n+1)=6.544488e+000; foe(n+1)=3.809934e+001; krok(n+1)=8.178765e-006; ng(n+1)=2.837687e+002;
n=274; farx(n+1)=6.575027e+000; foe(n+1)=3.807497e+001; krok(n+1)=2.114892e-005; ng(n+1)=1.688239e+002;
n=275; farx(n+1)=6.535647e+000; foe(n+1)=3.805090e+001; krok(n+1)=8.206521e-006; ng(n+1)=2.801130e+002;
n=276; farx(n+1)=6.566214e+000; foe(n+1)=3.802698e+001; krok(n+1)=2.122479e-005; ng(n+1)=1.680680e+002;
n=277; farx(n+1)=6.527084e+000; foe(n+1)=3.800301e+001; krok(n+1)=8.153596e-006; ng(n+1)=2.806829e+002;
n=278; farx(n+1)=6.557183e+000; foe(n+1)=3.797923e+001; krok(n+1)=2.132008e-005; ng(n+1)=1.669390e+002;
n=279; farx(n+1)=6.518533e+000; foe(n+1)=3.795562e+001; krok(n+1)=8.153596e-006; ng(n+1)=2.787773e+002;
n=280; farx(n+1)=6.548169e+000; foe(n+1)=3.793220e+001; krok(n+1)=2.114892e-005; ng(n+1)=1.661067e+002;
n=281; farx(n+1)=6.510107e+000; foe(n+1)=3.790908e+001; krok(n+1)=8.182517e-006; ng(n+1)=2.756093e+002;
n=282; farx(n+1)=6.539489e+000; foe(n+1)=3.788612e+001; krok(n+1)=2.097007e-005; ng(n+1)=1.654485e+002;
n=283; farx(n+1)=6.501839e+000; foe(n+1)=3.786337e+001; krok(n+1)=8.205260e-006; ng(n+1)=2.735028e+002;
n=284; farx(n+1)=6.531028e+000; foe(n+1)=3.784090e+001; krok(n+1)=2.077412e-005; ng(n+1)=1.648683e+002;
n=285; farx(n+1)=6.493841e+000; foe(n+1)=3.781849e+001; krok(n+1)=8.205260e-006; ng(n+1)=2.716350e+002;
n=286; farx(n+1)=6.522610e+000; foe(n+1)=3.779633e+001; krok(n+1)=2.063670e-005; ng(n+1)=1.640633e+002;
n=287; farx(n+1)=6.485974e+000; foe(n+1)=3.777433e+001; krok(n+1)=8.223627e-006; ng(n+1)=2.688586e+002;
n=288; farx(n+1)=6.514398e+000; foe(n+1)=3.775252e+001; krok(n+1)=2.047978e-005; ng(n+1)=1.633464e+002;
n=289; farx(n+1)=6.478309e+000; foe(n+1)=3.773091e+001; krok(n+1)=8.229150e-006; ng(n+1)=2.663828e+002;
n=290; farx(n+1)=6.506474e+000; foe(n+1)=3.770935e+001; krok(n+1)=2.047978e-005; ng(n+1)=1.625300e+002;
n=291; farx(n+1)=6.470657e+000; foe(n+1)=3.768800e+001; krok(n+1)=8.245016e-006; ng(n+1)=2.651166e+002;
n=292; farx(n+1)=6.498545e+000; foe(n+1)=3.766689e+001; krok(n+1)=2.020969e-005; ng(n+1)=1.620108e+002;
n=293; farx(n+1)=6.463322e+000; foe(n+1)=3.764592e+001; krok(n+1)=8.245594e-006; ng(n+1)=2.624601e+002;
n=294; farx(n+1)=6.491068e+000; foe(n+1)=3.762498e+001; krok(n+1)=2.036925e-005; ng(n+1)=1.611237e+002;
n=295; farx(n+1)=6.456079e+000; foe(n+1)=3.760413e+001; krok(n+1)=8.206521e-006; ng(n+1)=2.625763e+002;
n=296; farx(n+1)=6.483708e+000; foe(n+1)=3.758333e+001; krok(n+1)=2.058283e-005; ng(n+1)=1.602052e+002;
n=297; farx(n+1)=6.448932e+000; foe(n+1)=3.756255e+001; krok(n+1)=8.155065e-006; ng(n+1)=2.630691e+002;
n=298; farx(n+1)=6.476163e+000; foe(n+1)=3.754183e+001; krok(n+1)=2.067984e-005; ng(n+1)=1.592028e+002;
n=299; farx(n+1)=6.441725e+000; foe(n+1)=3.752132e+001; krok(n+1)=8.166804e-006; ng(n+1)=2.614868e+002;
n=300; farx(n+1)=6.468673e+000; foe(n+1)=3.750095e+001; krok(n+1)=2.047978e-005; ng(n+1)=1.586165e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)