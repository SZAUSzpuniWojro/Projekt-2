%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.314400e+003; foe(n+1)=4.422941e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=3.755043e+003; foe(n+1)=3.883368e+003; krok(n+1)=5.733303e-004; ng(n+1)=3.123445e+003;
n=2; farx(n+1)=1.384838e+003; foe(n+1)=1.541883e+003; krok(n+1)=2.678603e-003; ng(n+1)=2.362319e+003;
n=3; farx(n+1)=5.380479e+002; foe(n+1)=8.612495e+002; krok(n+1)=1.535829e-004; ng(n+1)=9.489461e+003;
n=4; farx(n+1)=5.258916e+002; foe(n+1)=7.740722e+002; krok(n+1)=6.459813e-004; ng(n+1)=1.472164e+003;
n=5; farx(n+1)=4.130942e+002; foe(n+1)=7.232027e+002; krok(n+1)=8.459568e-005; ng(n+1)=4.141327e+003;
n=6; farx(n+1)=4.078790e+002; foe(n+1)=6.932906e+002; krok(n+1)=2.860226e-004; ng(n+1)=1.143212e+003;
n=7; farx(n+1)=3.486736e+002; foe(n+1)=6.662725e+002; krok(n+1)=7.818173e-005; ng(n+1)=2.742160e+003;
n=8; farx(n+1)=3.564053e+002; foe(n+1)=6.433042e+002; krok(n+1)=2.093774e-004; ng(n+1)=1.566962e+003;
n=9; farx(n+1)=3.049196e+002; foe(n+1)=6.193429e+002; krok(n+1)=7.389767e-005; ng(n+1)=2.771109e+003;
n=10; farx(n+1)=3.096411e+002; foe(n+1)=5.975224e+002; krok(n+1)=1.874926e-004; ng(n+1)=1.467547e+003;
n=11; farx(n+1)=2.662970e+002; foe(n+1)=5.747239e+002; krok(n+1)=6.783250e-005; ng(n+1)=2.827433e+003;
n=12; farx(n+1)=2.702108e+002; foe(n+1)=5.538037e+002; krok(n+1)=1.733806e-004; ng(n+1)=1.555733e+003;
n=13; farx(n+1)=2.329402e+002; foe(n+1)=5.317313e+002; krok(n+1)=6.158863e-005; ng(n+1)=2.942487e+003;
n=14; farx(n+1)=2.351171e+002; foe(n+1)=5.110298e+002; krok(n+1)=1.662532e-004; ng(n+1)=1.569774e+003;
n=15; farx(n+1)=2.030458e+002; foe(n+1)=4.897701e+002; krok(n+1)=5.603333e-005; ng(n+1)=3.065103e+003;
n=16; farx(n+1)=2.046352e+002; foe(n+1)=4.703013e+002; krok(n+1)=1.526551e-004; ng(n+1)=1.608591e+003;
n=17; farx(n+1)=1.773825e+002; foe(n+1)=4.503767e+002; krok(n+1)=5.113027e-005; ng(n+1)=3.141858e+003;
n=18; farx(n+1)=1.785757e+002; foe(n+1)=4.322974e+002; krok(n+1)=1.387513e-004; ng(n+1)=1.639517e+003;
n=19; farx(n+1)=1.555356e+002; foe(n+1)=4.139723e+002; krok(n+1)=4.685521e-005; ng(n+1)=3.182499e+003;
n=20; farx(n+1)=1.565312e+002; foe(n+1)=3.974154e+002; krok(n+1)=1.251495e-004; ng(n+1)=1.673714e+003;
n=21; farx(n+1)=1.369690e+002; foe(n+1)=3.807596e+002; krok(n+1)=4.334515e-005; ng(n+1)=3.195958e+003;
n=22; farx(n+1)=1.380446e+002; foe(n+1)=3.660197e+002; krok(n+1)=1.094979e-004; ng(n+1)=1.709283e+003;
n=23; farx(n+1)=1.215787e+002; foe(n+1)=3.512401e+002; krok(n+1)=4.037383e-005; ng(n+1)=3.152009e+003;
n=24; farx(n+1)=1.226168e+002; foe(n+1)=3.379769e+002; krok(n+1)=9.814688e-005; ng(n+1)=1.737746e+003;
n=25; farx(n+1)=1.086188e+002; foe(n+1)=3.247178e+002; krok(n+1)=3.742939e-005; ng(n+1)=3.130898e+003;
n=26; farx(n+1)=1.094244e+002; foe(n+1)=3.124947e+002; krok(n+1)=9.144347e-005; ng(n+1)=1.743277e+003;
n=27; farx(n+1)=9.731624e+001; foe(n+1)=3.004313e+002; krok(n+1)=3.474555e-005; ng(n+1)=3.133113e+003;
n=28; farx(n+1)=9.800091e+001; foe(n+1)=2.893063e+002; krok(n+1)=8.459568e-005; ng(n+1)=1.743570e+003;
n=29; farx(n+1)=8.756048e+001; foe(n+1)=2.784101e+002; krok(n+1)=3.223481e-005; ng(n+1)=3.114869e+003;
n=30; farx(n+1)=8.807310e+001; foe(n+1)=2.681289e+002; krok(n+1)=8.155422e-005; ng(n+1)=1.727805e+003;
n=31; farx(n+1)=7.882822e+001; foe(n+1)=2.580222e+002; krok(n+1)=2.992557e-005; ng(n+1)=3.146650e+003;
n=32; farx(n+1)=7.929102e+001; foe(n+1)=2.487719e+002; krok(n+1)=7.528795e-005; ng(n+1)=1.722752e+003;
n=33; farx(n+1)=7.127509e+001; foe(n+1)=2.397474e+002; krok(n+1)=2.809223e-005; ng(n+1)=3.101358e+003;
n=34; farx(n+1)=7.169184e+001; foe(n+1)=2.314053e+002; krok(n+1)=7.046195e-005; ng(n+1)=1.711602e+003;
n=35; farx(n+1)=6.469011e+001; foe(n+1)=2.232769e+002; krok(n+1)=2.630754e-005; ng(n+1)=3.067664e+003;
n=36; farx(n+1)=6.501774e+001; foe(n+1)=2.156714e+002; krok(n+1)=6.703667e-005; ng(n+1)=1.688893e+003;
n=37; farx(n+1)=5.883816e+001; foe(n+1)=2.083228e+002; krok(n+1)=2.480406e-005; ng(n+1)=3.035096e+003;
n=38; farx(n+1)=5.916858e+001; foe(n+1)=2.015356e+002; krok(n+1)=6.234382e-005; ng(n+1)=1.671395e+003;
n=39; farx(n+1)=5.373837e+001; foe(n+1)=1.949713e+002; krok(n+1)=2.352777e-005; ng(n+1)=2.973573e+003;
n=40; farx(n+1)=5.406666e+001; foe(n+1)=1.888994e+002; krok(n+1)=5.843669e-005; ng(n+1)=1.648921e+003;
n=41; farx(n+1)=4.928094e+001; foe(n+1)=1.830057e+002; krok(n+1)=2.226608e-005; ng(n+1)=2.914941e+003;
n=42; farx(n+1)=4.956037e+001; foe(n+1)=1.774906e+002; krok(n+1)=5.575933e-005; ng(n+1)=1.615664e+003;
n=43; farx(n+1)=4.529918e+001; foe(n+1)=1.721758e+002; krok(n+1)=2.122479e-005; ng(n+1)=2.859404e+003;
n=44; farx(n+1)=4.559269e+001; foe(n+1)=1.672667e+002; krok(n+1)=5.210339e-005; ng(n+1)=1.589785e+003;
n=45; farx(n+1)=4.181941e+001; foe(n+1)=1.625097e+002; krok(n+1)=2.026093e-005; ng(n+1)=2.785445e+003;
n=46; farx(n+1)=4.208916e+001; foe(n+1)=1.580727e+002; krok(n+1)=4.960813e-005; ng(n+1)=1.555072e+003;
n=47; farx(n+1)=3.873206e+001; foe(n+1)=1.537850e+002; krok(n+1)=1.931405e-005; ng(n+1)=2.719590e+003;
n=48; farx(n+1)=3.896693e+001; foe(n+1)=1.497445e+002; krok(n+1)=4.780311e-005; ng(n+1)=1.515372e+003;
n=49; farx(n+1)=3.595876e+001; foe(n+1)=1.458584e+002; krok(n+1)=1.846630e-005; ng(n+1)=2.659750e+003;
n=50; farx(n+1)=3.617811e+001; foe(n+1)=1.422065e+002; krok(n+1)=4.570049e-005; ng(n+1)=1.477995e+003;
n=51; farx(n+1)=3.348316e+001; foe(n+1)=1.386981e+002; krok(n+1)=1.773106e-005; ng(n+1)=2.590555e+003;
n=52; farx(n+1)=3.368866e+001; foe(n+1)=1.354011e+002; krok(n+1)=4.341873e-005; ng(n+1)=1.440675e+003;
n=53; farx(n+1)=3.128314e+001; foe(n+1)=1.322563e+002; krok(n+1)=1.709205e-005; ng(n+1)=2.507455e+003;
n=54; farx(n+1)=3.147690e+001; foe(n+1)=1.292651e+002; krok(n+1)=4.211523e-005; ng(n+1)=1.401674e+003;
n=55; farx(n+1)=2.929613e+001; foe(n+1)=1.263886e+002; krok(n+1)=1.643953e-005; ng(n+1)=2.453467e+003;
n=56; farx(n+1)=2.948063e+001; foe(n+1)=1.236874e+002; krok(n+1)=4.014885e-005; ng(n+1)=1.364281e+003;
n=57; farx(n+1)=2.752177e+001; foe(n+1)=1.211017e+002; krok(n+1)=1.593982e-005; ng(n+1)=2.372307e+003;
n=58; farx(n+1)=2.770084e+001; foe(n+1)=1.186597e+002; krok(n+1)=3.862984e-005; ng(n+1)=1.327150e+003;
n=59; farx(n+1)=2.592187e+001; foe(n+1)=1.163103e+002; krok(n+1)=1.549259e-005; ng(n+1)=2.309052e+003;
n=60; farx(n+1)=2.610310e+001; foe(n+1)=1.141201e+002; krok(n+1)=3.650218e-005; ng(n+1)=1.294862e+003;
n=61; farx(n+1)=2.450179e+001; foe(n+1)=1.120121e+002; krok(n+1)=1.511881e-005; ng(n+1)=2.222871e+003;
n=62; farx(n+1)=2.467554e+001; foe(n+1)=1.100246e+002; krok(n+1)=3.514435e-005; ng(n+1)=1.259835e+003;
n=63; farx(n+1)=2.321829e+001; foe(n+1)=1.081093e+002; krok(n+1)=1.476011e-005; ng(n+1)=2.155604e+003;
n=64; farx(n+1)=2.339009e+001; foe(n+1)=1.063124e+002; krok(n+1)=3.369455e-005; ng(n+1)=1.227291e+003;
n=65; farx(n+1)=2.206361e+001; foe(n+1)=1.045718e+002; krok(n+1)=1.442554e-005; ng(n+1)=2.087231e+003;
n=66; farx(n+1)=2.222841e+001; foe(n+1)=1.029374e+002; krok(n+1)=3.242569e-005; ng(n+1)=1.193719e+003;
n=67; farx(n+1)=2.101836e+001; foe(n+1)=1.013551e+002; krok(n+1)=1.414532e-005; ng(n+1)=2.018051e+003;
n=68; farx(n+1)=2.117735e+001; foe(n+1)=9.986912e+001; krok(n+1)=3.104372e-005; ng(n+1)=1.161688e+003;
n=69; farx(n+1)=2.007589e+001; foe(n+1)=9.843530e+001; krok(n+1)=1.390226e-005; ng(n+1)=1.943172e+003;
n=70; farx(n+1)=2.022751e+001; foe(n+1)=9.707587e+001; krok(n+1)=3.008156e-005; ng(n+1)=1.128956e+003;
n=71; farx(n+1)=1.921940e+001; foe(n+1)=9.576371e+001; krok(n+1)=1.361563e-005; ng(n+1)=1.883231e+003;
n=72; farx(n+1)=1.936082e+001; foe(n+1)=9.451158e+001; krok(n+1)=2.941717e-005; ng(n+1)=1.095050e+003;
n=73; farx(n+1)=1.843123e+001; foe(n+1)=9.330475e+001; krok(n+1)=1.338477e-005; ng(n+1)=1.830730e+003;
n=74; farx(n+1)=1.856850e+001; foe(n+1)=9.216061e+001; krok(n+1)=2.832743e-005; ng(n+1)=1.065373e+003;
n=75; farx(n+1)=1.771499e+001; foe(n+1)=9.105825e+001; krok(n+1)=1.321254e-005; ng(n+1)=1.767422e+003;
n=76; farx(n+1)=1.784902e+001; foe(n+1)=9.000811e+001; krok(n+1)=2.769106e-005; ng(n+1)=1.035842e+003;
n=77; farx(n+1)=1.705522e+001; foe(n+1)=8.898643e+001; krok(n+1)=1.302585e-005; ng(n+1)=1.723669e+003;
n=78; farx(n+1)=1.718612e+001; foe(n+1)=8.802690e+001; krok(n+1)=2.639619e-005; ng(n+1)=1.009456e+003;
n=79; farx(n+1)=1.645426e+001; foe(n+1)=8.709823e+001; krok(n+1)=1.302585e-005; ng(n+1)=1.652208e+003;
n=80; farx(n+1)=1.658596e+001; foe(n+1)=8.622446e+001; krok(n+1)=2.503699e-005; ng(n+1)=9.868852e+002;
n=81; farx(n+1)=1.590989e+001; foe(n+1)=8.537691e+001; krok(n+1)=1.302585e-005; ng(n+1)=1.584809e+003;
n=82; farx(n+1)=1.604184e+001; foe(n+1)=8.457534e+001; krok(n+1)=2.418894e-005; ng(n+1)=9.639527e+002;
n=83; farx(n+1)=1.541503e+001; foe(n+1)=8.378952e+001; krok(n+1)=1.284475e-005; ng(n+1)=1.538205e+003;
n=84; farx(n+1)=1.553781e+001; foe(n+1)=8.304092e+001; krok(n+1)=2.375465e-005; ng(n+1)=9.353821e+002;
n=85; farx(n+1)=1.495535e+001; foe(n+1)=8.231212e+001; krok(n+1)=1.270754e-005; ng(n+1)=1.491543e+003;
n=86; farx(n+1)=1.507207e+001; foe(n+1)=8.161427e+001; krok(n+1)=2.334786e-005; ng(n+1)=9.091161e+002;
n=87; farx(n+1)=1.452887e+001; foe(n+1)=8.093510e+001; krok(n+1)=1.255575e-005; ng(n+1)=1.450868e+003;
n=88; farx(n+1)=1.463992e+001; foe(n+1)=8.028332e+001; krok(n+1)=2.299949e-005; ng(n+1)=8.835198e+002;
n=89; farx(n+1)=1.413071e+001; foe(n+1)=7.964895e+001; krok(n+1)=1.245057e-005; ng(n+1)=1.412697e+003;
n=90; farx(n+1)=1.423786e+001; foe(n+1)=7.904293e+001; krok(n+1)=2.230894e-005; ng(n+1)=8.611586e+002;
n=91; farx(n+1)=1.376551e+001; foe(n+1)=7.845473e+001; krok(n+1)=1.234279e-005; ng(n+1)=1.364243e+003;
n=92; farx(n+1)=1.386809e+001; foe(n+1)=7.788417e+001; krok(n+1)=2.253084e-005; ng(n+1)=8.353647e+002;
n=93; farx(n+1)=1.342422e+001; foe(n+1)=7.732396e+001; krok(n+1)=1.197415e-005; ng(n+1)=1.349809e+003;
n=94; farx(n+1)=1.351598e+001; foe(n+1)=7.677909e+001; krok(n+1)=2.286087e-005; ng(n+1)=8.053694e+002;
n=95; farx(n+1)=1.309622e+001; foe(n+1)=7.625015e+001; krok(n+1)=1.184623e-005; ng(n+1)=1.323719e+003;
n=96; farx(n+1)=1.318720e+001; foe(n+1)=7.574261e+001; krok(n+1)=2.253084e-005; ng(n+1)=7.851494e+002;
n=97; farx(n+1)=1.279369e+001; foe(n+1)=7.524405e+001; krok(n+1)=1.164400e-005; ng(n+1)=1.294864e+003;
n=98; farx(n+1)=1.287676e+001; foe(n+1)=7.476283e+001; krok(n+1)=2.230894e-005; ng(n+1)=7.600918e+002;
n=99; farx(n+1)=1.250812e+001; foe(n+1)=7.429834e+001; krok(n+1)=1.161523e-005; ng(n+1)=1.252786e+003;
n=100; farx(n+1)=1.259030e+001; foe(n+1)=7.384814e+001; krok(n+1)=2.208984e-005; ng(n+1)=7.407364e+002;
n=101; farx(n+1)=1.224313e+001; foe(n+1)=7.340866e+001; krok(n+1)=1.143043e-005; ng(n+1)=1.227679e+003;
n=102; farx(n+1)=1.232099e+001; foe(n+1)=7.298163e+001; krok(n+1)=2.226248e-005; ng(n+1)=7.179888e+002;
n=103; farx(n+1)=1.199487e+001; foe(n+1)=7.256382e+001; krok(n+1)=1.115447e-005; ng(n+1)=1.207281e+003;
n=104; farx(n+1)=1.206570e+001; foe(n+1)=7.215253e+001; krok(n+1)=2.296646e-005; ng(n+1)=6.924566e+002;
n=105; farx(n+1)=1.175711e+001; foe(n+1)=7.175155e+001; krok(n+1)=1.085468e-005; ng(n+1)=1.195787e+003;
n=106; farx(n+1)=1.182219e+001; foe(n+1)=7.135523e+001; krok(n+1)=2.382783e-005; ng(n+1)=6.684726e+002;
n=107; farx(n+1)=1.152803e+001; foe(n+1)=7.096756e+001; krok(n+1)=1.057446e-005; ng(n+1)=1.189826e+003;
n=108; farx(n+1)=1.158875e+001; foe(n+1)=7.058788e+001; krok(n+1)=2.444110e-005; ng(n+1)=6.472912e+002;
n=109; farx(n+1)=1.130941e+001; foe(n+1)=7.021532e+001; krok(n+1)=1.033992e-005; ng(n+1)=1.177293e+003;
n=110; farx(n+1)=1.136534e+001; foe(n+1)=6.985144e+001; krok(n+1)=2.481646e-005; ng(n+1)=6.272045e+002;
n=111; farx(n+1)=1.109892e+001; foe(n+1)=6.949760e+001; krok(n+1)=1.029142e-005; ng(n+1)=1.154336e+003;
n=112; farx(n+1)=1.115592e+001; foe(n+1)=6.916124e+001; krok(n+1)=2.395025e-005; ng(n+1)=6.124488e+002;
n=113; farx(n+1)=1.090611e+001; foe(n+1)=6.883251e+001; krok(n+1)=1.029142e-005; ng(n+1)=1.113375e+003;
n=114; farx(n+1)=1.096212e+001; foe(n+1)=6.851736e+001; krok(n+1)=2.343658e-005; ng(n+1)=5.969742e+002;
n=115; farx(n+1)=1.072668e+001; foe(n+1)=6.821025e+001; krok(n+1)=1.029142e-005; ng(n+1)=1.077373e+003;
n=116; farx(n+1)=1.078195e+001; foe(n+1)=6.791521e+001; krok(n+1)=2.286087e-005; ng(n+1)=5.824165e+002;
n=117; farx(n+1)=1.056123e+001; foe(n+1)=6.762840e+001; krok(n+1)=1.023989e-005; ng(n+1)=1.041411e+003;
n=118; farx(n+1)=1.061471e+001; foe(n+1)=6.734737e+001; krok(n+1)=2.314368e-005; ng(n+1)=5.666354e+002;
n=119; farx(n+1)=1.040386e+001; foe(n+1)=6.707252e+001; krok(n+1)=1.010485e-005; ng(n+1)=1.026067e+003;
n=120; farx(n+1)=1.045573e+001; foe(n+1)=6.680564e+001; krok(n+1)=2.322634e-005; ng(n+1)=5.515711e+002;
n=121; farx(n+1)=1.025558e+001; foe(n+1)=6.654403e+001; krok(n+1)=9.975945e-006; ng(n+1)=1.005519e+003;
n=122; farx(n+1)=1.030488e+001; foe(n+1)=6.628870e+001; krok(n+1)=2.343658e-005; ng(n+1)=5.362506e+002;
n=123; farx(n+1)=1.011320e+001; foe(n+1)=6.603993e+001; krok(n+1)=9.952471e-006; ng(n+1)=9.844521e+002;
n=124; farx(n+1)=1.016278e+001; foe(n+1)=6.580205e+001; krok(n+1)=2.276873e-005; ng(n+1)=5.238405e+002;
n=125; farx(n+1)=9.982830e+000; foe(n+1)=6.556934e+001; krok(n+1)=9.888364e-006; ng(n+1)=9.519255e+002;
n=126; farx(n+1)=1.003030e+001; foe(n+1)=6.534094e+001; krok(n+1)=2.322634e-005; ng(n+1)=5.093906e+002;
n=127; farx(n+1)=9.856958e+000; foe(n+1)=6.511731e+001; krok(n+1)=9.809075e-006; ng(n+1)=9.386680e+002;
n=128; farx(n+1)=9.903728e+000; foe(n+1)=6.490276e+001; krok(n+1)=2.273062e-005; ng(n+1)=4.970256e+002;
n=129; farx(n+1)=9.740345e+000; foe(n+1)=6.469338e+001; krok(n+1)=9.772716e-006; ng(n+1)=9.078597e+002;
n=130; farx(n+1)=9.785494e+000; foe(n+1)=6.448860e+001; krok(n+1)=2.286087e-005; ng(n+1)=4.838502e+002;
n=131; farx(n+1)=9.628444e+000; foe(n+1)=6.428932e+001; krok(n+1)=9.766516e-006; ng(n+1)=8.885392e+002;
n=132; farx(n+1)=9.673959e+000; foe(n+1)=6.409841e+001; krok(n+1)=2.226248e-005; ng(n+1)=4.729088e+002;
n=133; farx(n+1)=9.524839e+000; foe(n+1)=6.391117e+001; krok(n+1)=9.766516e-006; ng(n+1)=8.610526e+002;
n=134; farx(n+1)=9.569444e+000; foe(n+1)=6.373074e+001; krok(n+1)=2.198334e-005; ng(n+1)=4.613628e+002;
n=135; farx(n+1)=9.427311e+000; foe(n+1)=6.355418e+001; krok(n+1)=9.766516e-006; ng(n+1)=8.363682e+002;
n=136; farx(n+1)=9.470560e+000; foe(n+1)=6.338394e+001; krok(n+1)=2.139878e-005; ng(n+1)=4.503430e+002;
n=137; farx(n+1)=9.335484e+000; foe(n+1)=6.321961e+001; krok(n+1)=9.884595e-006; ng(n+1)=8.037961e+002;
n=138; farx(n+1)=9.379692e+000; foe(n+1)=6.306104e+001; krok(n+1)=2.092147e-005; ng(n+1)=4.409723e+002;
n=139; farx(n+1)=9.250761e+000; foe(n+1)=6.290501e+001; krok(n+1)=9.820196e-006; ng(n+1)=7.842580e+002;
n=140; farx(n+1)=9.293544e+000; foe(n+1)=6.275333e+001; krok(n+1)=2.103720e-005; ng(n+1)=4.296942e+002;
n=141; farx(n+1)=9.170343e+000; foe(n+1)=6.260456e+001; krok(n+1)=9.716744e-006; ng(n+1)=7.675971e+002;
n=142; farx(n+1)=9.212166e+000; foe(n+1)=6.245861e+001; krok(n+1)=2.167257e-005; ng(n+1)=4.181618e+002;
n=143; farx(n+1)=9.092651e+000; foe(n+1)=6.231385e+001; krok(n+1)=9.540944e-006; ng(n+1)=7.631496e+002;
n=144; farx(n+1)=9.132258e+000; foe(n+1)=6.217349e+001; krok(n+1)=2.170937e-005; ng(n+1)=4.068044e+002;
n=145; farx(n+1)=9.018077e+000; foe(n+1)=6.203701e+001; krok(n+1)=9.549067e-006; ng(n+1)=7.406081e+002;
n=146; farx(n+1)=9.057539e+000; foe(n+1)=6.190438e+001; krok(n+1)=2.167257e-005; ng(n+1)=3.969939e+002;
n=147; farx(n+1)=8.947726e+000; foe(n+1)=6.177388e+001; krok(n+1)=9.471564e-006; ng(n+1)=7.262395e+002;
n=148; farx(n+1)=8.986397e+000; foe(n+1)=6.164726e+001; krok(n+1)=2.190752e-005; ng(n+1)=3.866757e+002;
n=149; farx(n+1)=8.880739e+000; foe(n+1)=6.152203e+001; krok(n+1)=9.357347e-006; ng(n+1)=7.140164e+002;
n=150; farx(n+1)=8.917903e+000; foe(n+1)=6.140023e+001; krok(n+1)=2.219138e-005; ng(n+1)=3.760284e+002;
n=151; farx(n+1)=8.816330e+000; foe(n+1)=6.128073e+001; krok(n+1)=9.290308e-006; ng(n+1)=6.989938e+002;
n=152; farx(n+1)=8.852450e+000; foe(n+1)=6.116448e+001; krok(n+1)=2.226608e-005; ng(n+1)=3.660713e+002;
n=153; farx(n+1)=8.754818e+000; foe(n+1)=6.105087e+001; krok(n+1)=9.254712e-006; ng(n+1)=6.824607e+002;
n=154; farx(n+1)=8.790335e+000; foe(n+1)=6.094053e+001; krok(n+1)=2.226248e-005; ng(n+1)=3.566764e+002;
n=155; farx(n+1)=8.696432e+000; foe(n+1)=6.083236e+001; krok(n+1)=9.207477e-006; ng(n+1)=6.671004e+002;
n=156; farx(n+1)=8.730981e+000; foe(n+1)=6.072735e+001; krok(n+1)=2.219549e-005; ng(n+1)=3.473359e+002;
n=157; farx(n+1)=8.640868e+000; foe(n+1)=6.062498e+001; krok(n+1)=9.195979e-006; ng(n+1)=6.490992e+002;
n=158; farx(n+1)=8.675011e+000; foe(n+1)=6.052548e+001; krok(n+1)=2.219138e-005; ng(n+1)=3.384987e+002;
n=159; farx(n+1)=8.588195e+000; foe(n+1)=6.042781e+001; krok(n+1)=9.139546e-006; ng(n+1)=6.353446e+002;
n=160; farx(n+1)=8.621356e+000; foe(n+1)=6.033289e+001; krok(n+1)=2.219549e-005; ng(n+1)=3.294848e+002;
n=161; farx(n+1)=8.537866e+000; foe(n+1)=6.024029e+001; krok(n+1)=9.125546e-006; ng(n+1)=6.192662e+002;
n=162; farx(n+1)=8.570513e+000; foe(n+1)=6.015034e+001; krok(n+1)=2.212034e-005; ng(n+1)=3.211466e+002;
n=163; farx(n+1)=8.490166e+000; foe(n+1)=6.006231e+001; krok(n+1)=9.089066e-006; ng(n+1)=6.047465e+002;
n=164; farx(n+1)=8.522131e+000; foe(n+1)=5.997662e+001; krok(n+1)=2.219138e-005; ng(n+1)=3.129287e+002;
n=165; farx(n+1)=8.444603e+000; foe(n+1)=5.989270e+001; krok(n+1)=9.049168e-006; ng(n+1)=5.917352e+002;
n=166; farx(n+1)=8.475714e+000; foe(n+1)=5.981114e+001; krok(n+1)=2.209045e-005; ng(n+1)=3.049787e+002;
n=167; farx(n+1)=8.401171e+000; foe(n+1)=5.973177e+001; krok(n+1)=9.050777e-006; ng(n+1)=5.754913e+002;
n=168; farx(n+1)=8.431841e+000; foe(n+1)=5.965451e+001; krok(n+1)=2.199768e-005; ng(n+1)=2.975714e+002;
n=169; farx(n+1)=8.359898e+000; foe(n+1)=5.957902e+001; krok(n+1)=9.031324e-006; ng(n+1)=5.622990e+002;
n=170; farx(n+1)=8.389687e+000; foe(n+1)=5.950567e+001; krok(n+1)=2.170937e-005; ng(n+1)=2.902472e+002;
n=171; farx(n+1)=8.320608e+000; foe(n+1)=5.943472e+001; krok(n+1)=9.081248e-006; ng(n+1)=5.441654e+002;
n=172; farx(n+1)=8.350572e+000; foe(n+1)=5.936566e+001; krok(n+1)=2.167257e-005; ng(n+1)=2.837390e+002;
n=173; farx(n+1)=8.283577e+000; foe(n+1)=5.929741e+001; krok(n+1)=8.996611e-006; ng(n+1)=5.361018e+002;
n=174; farx(n+1)=8.312403e+000; foe(n+1)=5.923105e+001; krok(n+1)=2.170937e-005; ng(n+1)=2.761369e+002;
n=175; farx(n+1)=8.247871e+000; foe(n+1)=5.916645e+001; krok(n+1)=9.006227e-006; ng(n+1)=5.214557e+002;
n=176; farx(n+1)=8.275982e+000; foe(n+1)=5.910351e+001; krok(n+1)=2.139878e-005; ng(n+1)=2.695383e+002;
n=177; farx(n+1)=8.213849e+000; foe(n+1)=5.904267e+001; krok(n+1)=9.062615e-006; ng(n+1)=5.050843e+002;
n=178; farx(n+1)=8.241997e+000; foe(n+1)=5.898341e+001; krok(n+1)=2.122479e-005; ng(n+1)=2.636700e+002;
n=179; farx(n+1)=8.181804e+000; foe(n+1)=5.892516e+001; krok(n+1)=9.014671e-006; ng(n+1)=4.954175e+002;
n=180; farx(n+1)=8.209046e+000; foe(n+1)=5.886840e+001; krok(n+1)=2.114892e-005; ng(n+1)=2.569716e+002;
n=181; farx(n+1)=8.150922e+000; foe(n+1)=5.881324e+001; krok(n+1)=9.049168e-006; ng(n+1)=4.815583e+002;
n=182; farx(n+1)=8.177900e+000; foe(n+1)=5.875955e+001; krok(n+1)=2.092658e-005; ng(n+1)=2.512473e+002;
n=183; farx(n+1)=8.121677e+000; foe(n+1)=5.870701e+001; krok(n+1)=9.041389e-006; ng(n+1)=4.701481e+002;
n=184; farx(n+1)=8.148023e+000; foe(n+1)=5.865579e+001; krok(n+1)=2.078166e-005; ng(n+1)=2.452462e+002;
n=185; farx(n+1)=8.093616e+000; foe(n+1)=5.860590e+001; krok(n+1)=9.071830e-006; ng(n+1)=4.577228e+002;
n=186; farx(n+1)=8.119550e+000; foe(n+1)=5.855735e+001; krok(n+1)=2.047978e-005; ng(n+1)=2.398456e+002;
n=187; farx(n+1)=8.066981e+000; foe(n+1)=5.851001e+001; krok(n+1)=9.098641e-006; ng(n+1)=4.451737e+002;
n=188; farx(n+1)=8.092619e+000; foe(n+1)=5.846379e+001; krok(n+1)=2.038856e-005; ng(n+1)=2.344544e+002;
n=189; farx(n+1)=8.041579e+000; foe(n+1)=5.841846e+001; krok(n+1)=9.078452e-006; ng(n+1)=4.360649e+002;
n=190; farx(n+1)=8.066694e+000; foe(n+1)=5.837422e+001; krok(n+1)=2.034152e-005; ng(n+1)=2.289240e+002;
n=191; farx(n+1)=8.017163e+000; foe(n+1)=5.833090e+001; krok(n+1)=9.078452e-006; ng(n+1)=4.262628e+002;
n=192; farx(n+1)=8.041811e+000; foe(n+1)=5.828864e+001; krok(n+1)=2.020969e-005; ng(n+1)=2.237240e+002;
n=193; farx(n+1)=7.993810e+000; foe(n+1)=5.824730e+001; krok(n+1)=9.081828e-006; ng(n+1)=4.162828e+002;
n=194; farx(n+1)=8.018053e+000; foe(n+1)=5.820685e+001; krok(n+1)=2.016710e-005; ng(n+1)=2.188169e+002;
n=195; farx(n+1)=7.971356e+000; foe(n+1)=5.816722e+001; krok(n+1)=9.078452e-006; ng(n+1)=4.080623e+002;
n=196; farx(n+1)=7.995371e+000; foe(n+1)=5.812853e+001; krok(n+1)=2.018691e-005; ng(n+1)=2.142196e+002;
n=197; farx(n+1)=7.949906e+000; foe(n+1)=5.809033e+001; krok(n+1)=9.024321e-006; ng(n+1)=4.017101e+002;
n=198; farx(n+1)=7.973349e+000; foe(n+1)=5.805297e+001; krok(n+1)=2.030432e-005; ng(n+1)=2.092771e+002;
n=199; farx(n+1)=7.929095e+000; foe(n+1)=5.801627e+001; krok(n+1)=9.006227e-006; ng(n+1)=3.942957e+002;
n=200; farx(n+1)=7.952122e+000; foe(n+1)=5.798037e+001; krok(n+1)=2.028847e-005; ng(n+1)=2.048644e+002;
n=201; farx(n+1)=7.909064e+000; foe(n+1)=5.794513e+001; krok(n+1)=8.995004e-006; ng(n+1)=3.866830e+002;
n=202; farx(n+1)=7.931746e+000; foe(n+1)=5.791064e+001; krok(n+1)=2.030432e-005; ng(n+1)=2.006077e+002;
n=203; farx(n+1)=7.889839e+000; foe(n+1)=5.787671e+001; krok(n+1)=8.960031e-006; ng(n+1)=3.799560e+002;
n=204; farx(n+1)=7.912107e+000; foe(n+1)=5.784338e+001; krok(n+1)=2.044242e-005; ng(n+1)=1.962648e+002;
n=205; farx(n+1)=7.871130e+000; foe(n+1)=5.781063e+001; krok(n+1)=8.945225e-006; ng(n+1)=3.739893e+002;
n=206; farx(n+1)=7.893069e+000; foe(n+1)=5.777864e+001; krok(n+1)=2.034511e-005; ng(n+1)=1.925453e+002;
n=207; farx(n+1)=7.853196e+000; foe(n+1)=5.774717e+001; krok(n+1)=8.938206e-006; ng(n+1)=3.666290e+002;
n=208; farx(n+1)=7.874806e+000; foe(n+1)=5.771634e+001; krok(n+1)=2.036925e-005; ng(n+1)=1.887278e+002;
n=209; farx(n+1)=7.835840e+000; foe(n+1)=5.768598e+001; krok(n+1)=8.932063e-006; ng(n+1)=3.604616e+002;
n=210; farx(n+1)=7.857075e+000; foe(n+1)=5.765634e+001; krok(n+1)=2.020969e-005; ng(n+1)=1.852434e+002;
n=211; farx(n+1)=7.819135e+000; foe(n+1)=5.762722e+001; krok(n+1)=8.958287e-006; ng(n+1)=3.526187e+002;
n=212; farx(n+1)=7.840127e+000; foe(n+1)=5.759876e+001; krok(n+1)=2.006817e-005; ng(n+1)=1.820312e+002;
n=213; farx(n+1)=7.803192e+000; foe(n+1)=5.757072e+001; krok(n+1)=8.943311e-006; ng(n+1)=3.459722e+002;
n=214; farx(n+1)=7.823897e+000; foe(n+1)=5.754313e+001; krok(n+1)=2.020511e-005; ng(n+1)=1.787357e+002;
n=215; farx(n+1)=7.787800e+000; foe(n+1)=5.751591e+001; krok(n+1)=8.887598e-006; ng(n+1)=3.414081e+002;
n=216; farx(n+1)=7.808283e+000; foe(n+1)=5.748902e+001; krok(n+1)=2.058283e-005; ng(n+1)=1.754530e+002;
n=217; farx(n+1)=7.772720e+000; foe(n+1)=5.746236e+001; krok(n+1)=8.818436e-006; ng(n+1)=3.391499e+002;
n=218; farx(n+1)=7.792847e+000; foe(n+1)=5.743623e+001; krok(n+1)=2.068400e-005; ng(n+1)=1.725098e+002;
n=219; farx(n+1)=7.758089e+000; foe(n+1)=5.741041e+001; krok(n+1)=8.794029e-006; ng(n+1)=3.338990e+002;
n=220; farx(n+1)=7.777823e+000; foe(n+1)=5.738502e+001; krok(n+1)=2.067984e-005; ng(n+1)=1.698161e+002;
n=221; farx(n+1)=7.743860e+000; foe(n+1)=5.736011e+001; krok(n+1)=8.807743e-006; ng(n+1)=3.277161e+002;
n=222; farx(n+1)=7.763446e+000; foe(n+1)=5.733561e+001; krok(n+1)=2.063495e-005; ng(n+1)=1.677632e+002;
n=223; farx(n+1)=7.730141e+000; foe(n+1)=5.731142e+001; krok(n+1)=8.799185e-006; ng(n+1)=3.230066e+002;
n=224; farx(n+1)=7.749401e+000; foe(n+1)=5.728771e+001; krok(n+1)=2.047978e-005; ng(n+1)=1.656512e+002;
n=225; farx(n+1)=7.716881e+000; foe(n+1)=5.726440e+001; krok(n+1)=8.829868e-006; ng(n+1)=3.164640e+002;
n=226; farx(n+1)=7.736004e+000; foe(n+1)=5.724150e+001; krok(n+1)=2.036925e-005; ng(n+1)=1.638555e+002;
n=227; farx(n+1)=7.704124e+000; foe(n+1)=5.721887e+001; krok(n+1)=8.829868e-006; ng(n+1)=3.117585e+002;
n=228; farx(n+1)=7.723010e+000; foe(n+1)=5.719667e+001; krok(n+1)=2.027520e-005; ng(n+1)=1.619129e+002;
n=229; farx(n+1)=7.691800e+000; foe(n+1)=5.717474e+001; krok(n+1)=8.838430e-006; ng(n+1)=3.066965e+002;
n=230; farx(n+1)=7.710445e+000; foe(n+1)=5.715319e+001; krok(n+1)=2.016710e-005; ng(n+1)=1.600395e+002;
n=231; farx(n+1)=7.679951e+000; foe(n+1)=5.713194e+001; krok(n+1)=8.833238e-006; ng(n+1)=3.016167e+002;
n=232; farx(n+1)=7.698438e+000; foe(n+1)=5.711090e+001; krok(n+1)=2.031518e-005; ng(n+1)=1.579176e+002;
n=233; farx(n+1)=7.668379e+000; foe(n+1)=5.709009e+001; krok(n+1)=8.807743e-006; ng(n+1)=2.990324e+002;
n=234; farx(n+1)=7.686659e+000; foe(n+1)=5.706959e+001; krok(n+1)=2.028847e-005; ng(n+1)=1.560907e+002;
n=235; farx(n+1)=7.657155e+000; foe(n+1)=5.704934e+001; krok(n+1)=8.807743e-006; ng(n+1)=2.949974e+002;
n=236; farx(n+1)=7.675218e+000; foe(n+1)=5.702937e+001; krok(n+1)=2.020969e-005; ng(n+1)=1.545737e+002;
n=237; farx(n+1)=7.646243e+000; foe(n+1)=5.700968e+001; krok(n+1)=8.829868e-006; ng(n+1)=2.906354e+002;
n=238; farx(n+1)=7.664184e+000; foe(n+1)=5.699030e+001; krok(n+1)=2.006817e-005; ng(n+1)=1.534151e+002;
n=239; farx(n+1)=7.635788e+000; foe(n+1)=5.697111e+001; krok(n+1)=8.818566e-006; ng(n+1)=2.866493e+002;
n=240; farx(n+1)=7.653579e+000; foe(n+1)=5.695212e+001; krok(n+1)=2.020511e-005; ng(n+1)=1.518260e+002;
n=241; farx(n+1)=7.625554e+000; foe(n+1)=5.693330e+001; krok(n+1)=8.799185e-006; ng(n+1)=2.843420e+002;
n=242; farx(n+1)=7.643086e+000; foe(n+1)=5.691474e+001; krok(n+1)=2.007443e-005; ng(n+1)=1.504918e+002;
n=243; farx(n+1)=7.615585e+000; foe(n+1)=5.689645e+001; krok(n+1)=8.838430e-006; ng(n+1)=2.796052e+002;
n=244; farx(n+1)=7.633033e+000; foe(n+1)=5.687842e+001; krok(n+1)=1.991107e-005; ng(n+1)=1.495686e+002;
n=245; farx(n+1)=7.605961e+000; foe(n+1)=5.686057e+001; krok(n+1)=8.856543e-006; ng(n+1)=2.760322e+002;
n=246; farx(n+1)=7.623204e+000; foe(n+1)=5.684300e+001; krok(n+1)=1.968224e-005; ng(n+1)=1.485880e+002;
n=247; farx(n+1)=7.596612e+000; foe(n+1)=5.682564e+001; krok(n+1)=8.913402e-006; ng(n+1)=2.712787e+002;
n=248; farx(n+1)=7.613767e+000; foe(n+1)=5.680858e+001; krok(n+1)=1.941361e-005; ng(n+1)=1.479212e+002;
n=249; farx(n+1)=7.587720e+000; foe(n+1)=5.679167e+001; krok(n+1)=8.914462e-006; ng(n+1)=2.672593e+002;
n=250; farx(n+1)=7.604819e+000; foe(n+1)=5.677487e+001; krok(n+1)=1.961815e-005; ng(n+1)=1.465531e+002;
n=251; farx(n+1)=7.579018e+000; foe(n+1)=5.675814e+001; krok(n+1)=8.865531e-006; ng(n+1)=2.664626e+002;
n=252; farx(n+1)=7.595945e+000; foe(n+1)=5.674160e+001; krok(n+1)=1.969916e-005; ng(n+1)=1.452089e+002;
n=253; farx(n+1)=7.570554e+000; foe(n+1)=5.672517e+001; krok(n+1)=8.833238e-006; ng(n+1)=2.640882e+002;
n=254; farx(n+1)=7.587348e+000; foe(n+1)=5.670881e+001; krok(n+1)=1.991107e-005; ng(n+1)=1.437493e+002;
n=255; farx(n+1)=7.562202e+000; foe(n+1)=5.669259e+001; krok(n+1)=8.807743e-006; ng(n+1)=2.627910e+002;
n=256; farx(n+1)=7.578946e+000; foe(n+1)=5.667651e+001; krok(n+1)=1.999909e-005; ng(n+1)=1.426565e+002;
n=257; farx(n+1)=7.554126e+000; foe(n+1)=5.666049e+001; krok(n+1)=8.756691e-006; ng(n+1)=2.613993e+002;
n=258; farx(n+1)=7.570685e+000; foe(n+1)=5.664454e+001; krok(n+1)=2.020969e-005; ng(n+1)=1.411625e+002;
n=259; farx(n+1)=7.546115e+000; foe(n+1)=5.662874e+001; krok(n+1)=8.746968e-006; ng(n+1)=2.597217e+002;
n=260; farx(n+1)=7.562580e+000; foe(n+1)=5.661308e+001; krok(n+1)=2.016710e-005; ng(n+1)=1.402339e+002;
n=261; farx(n+1)=7.538326e+000; foe(n+1)=5.659755e+001; krok(n+1)=8.746968e-006; ng(n+1)=2.573264e+002;
n=262; farx(n+1)=7.554746e+000; foe(n+1)=5.658215e+001; krok(n+1)=2.020511e-005; ng(n+1)=1.394491e+002;
n=263; farx(n+1)=7.530766e+000; foe(n+1)=5.656681e+001; krok(n+1)=8.720009e-006; ng(n+1)=2.558027e+002;
n=264; farx(n+1)=7.547080e+000; foe(n+1)=5.655158e+001; krok(n+1)=2.031518e-005; ng(n+1)=1.386039e+002;
n=265; farx(n+1)=7.523309e+000; foe(n+1)=5.653644e+001; krok(n+1)=8.718949e-006; ng(n+1)=2.543505e+002;
n=266; farx(n+1)=7.539537e+000; foe(n+1)=5.652149e+001; krok(n+1)=2.020511e-005; ng(n+1)=1.380555e+002;
n=267; farx(n+1)=7.516172e+000; foe(n+1)=5.650661e+001; krok(n+1)=8.690837e-006; ng(n+1)=2.519020e+002;
n=268; farx(n+1)=7.532364e+000; foe(n+1)=5.649173e+001; krok(n+1)=2.058283e-005; ng(n+1)=1.370322e+002;
n=269; farx(n+1)=7.509032e+000; foe(n+1)=5.647688e+001; krok(n+1)=8.657728e-006; ng(n+1)=2.524628e+002;
n=270; farx(n+1)=7.525090e+000; foe(n+1)=5.646224e+001; krok(n+1)=2.044242e-005; ng(n+1)=1.364751e+002;
n=271; farx(n+1)=7.502082e+000; foe(n+1)=5.644767e+001; krok(n+1)=8.682944e-006; ng(n+1)=2.495504e+002;
n=272; farx(n+1)=7.518059e+000; foe(n+1)=5.643328e+001; krok(n+1)=2.031518e-005; ng(n+1)=1.359914e+002;
n=273; farx(n+1)=7.495361e+000; foe(n+1)=5.641895e+001; krok(n+1)=8.686387e-006; ng(n+1)=2.471716e+002;
n=274; farx(n+1)=7.511243e+000; foe(n+1)=5.640472e+001; krok(n+1)=2.031935e-005; ng(n+1)=1.353559e+002;
n=275; farx(n+1)=7.488786e+000; foe(n+1)=5.639058e+001; krok(n+1)=8.686387e-006; ng(n+1)=2.454108e+002;
n=276; farx(n+1)=7.504609e+000; foe(n+1)=5.637654e+001; krok(n+1)=2.031518e-005; ng(n+1)=1.347934e+002;
n=277; farx(n+1)=7.482392e+000; foe(n+1)=5.636256e+001; krok(n+1)=8.671956e-006; ng(n+1)=2.439287e+002;
n=278; farx(n+1)=7.498118e+000; foe(n+1)=5.634865e+001; krok(n+1)=2.036880e-005; ng(n+1)=1.341247e+002;
n=279; farx(n+1)=7.476075e+000; foe(n+1)=5.633483e+001; krok(n+1)=8.682944e-006; ng(n+1)=2.424596e+002;
n=280; farx(n+1)=7.491730e+000; foe(n+1)=5.632114e+001; krok(n+1)=2.020969e-005; ng(n+1)=1.337442e+002;
n=281; farx(n+1)=7.469904e+000; foe(n+1)=5.630753e+001; krok(n+1)=8.718949e-006; ng(n+1)=2.401694e+002;
n=282; farx(n+1)=7.485528e+000; foe(n+1)=5.629410e+001; krok(n+1)=1.999909e-005; ng(n+1)=1.335052e+002;
n=283; farx(n+1)=7.464059e+000; foe(n+1)=5.628069e+001; krok(n+1)=8.690837e-006; ng(n+1)=2.379920e+002;
n=284; farx(n+1)=7.479616e+000; foe(n+1)=5.626725e+001; krok(n+1)=2.034511e-005; ng(n+1)=1.325953e+002;
n=285; farx(n+1)=7.458265e+000; foe(n+1)=5.625386e+001; krok(n+1)=8.632369e-006; ng(n+1)=2.383049e+002;
n=286; farx(n+1)=7.473813e+000; foe(n+1)=5.624044e+001; krok(n+1)=2.068400e-005; ng(n+1)=1.318012e+002;
n=287; farx(n+1)=7.452599e+000; foe(n+1)=5.622705e+001; krok(n+1)=8.546025e-006; ng(n+1)=2.390063e+002;
n=288; farx(n+1)=7.468042e+000; foe(n+1)=5.621354e+001; krok(n+1)=2.114892e-005; ng(n+1)=1.308245e+002;
n=289; farx(n+1)=7.446900e+000; foe(n+1)=5.620016e+001; krok(n+1)=8.500298e-006; ng(n+1)=2.394583e+002;
n=290; farx(n+1)=7.462558e+000; foe(n+1)=5.618670e+001; krok(n+1)=2.167257e-005; ng(n+1)=1.302118e+002;
n=291; farx(n+1)=7.441224e+000; foe(n+1)=5.617313e+001; krok(n+1)=8.423638e-006; ng(n+1)=2.427276e+002;
n=292; farx(n+1)=7.456764e+000; foe(n+1)=5.615976e+001; krok(n+1)=2.167257e-005; ng(n+1)=1.296877e+002;
n=293; farx(n+1)=7.435772e+000; foe(n+1)=5.614635e+001; krok(n+1)=8.379583e-006; ng(n+1)=2.409239e+002;
n=294; farx(n+1)=7.451283e+000; foe(n+1)=5.613286e+001; krok(n+1)=2.226248e-005; ng(n+1)=1.287811e+002;
n=295; farx(n+1)=7.430273e+000; foe(n+1)=5.611938e+001; krok(n+1)=8.322079e-006; ng(n+1)=2.424073e+002;
n=296; farx(n+1)=7.445620e+000; foe(n+1)=5.610592e+001; krok(n+1)=2.230894e-005; ng(n+1)=1.282199e+002;
n=297; farx(n+1)=7.424838e+000; foe(n+1)=5.609261e+001; krok(n+1)=8.333245e-006; ng(n+1)=2.404060e+002;
n=298; farx(n+1)=7.440170e+000; foe(n+1)=5.607929e+001; krok(n+1)=2.230894e-005; ng(n+1)=1.278822e+002;
n=299; farx(n+1)=7.419510e+000; foe(n+1)=5.606607e+001; krok(n+1)=8.334046e-006; ng(n+1)=2.395058e+002;
n=300; farx(n+1)=7.435022e+000; foe(n+1)=5.605290e+001; krok(n+1)=2.253084e-005; ng(n+1)=1.275923e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
