%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.826699e+003; foe(n+1)=4.563749e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=4.344247e+003; foe(n+1)=3.767270e+003; krok(n+1)=6.106204e-004; ng(n+1)=4.041246e+003;
n=2; farx(n+1)=1.625138e+003; foe(n+1)=1.068633e+003; krok(n+1)=2.053340e-003; ng(n+1)=4.060366e+003;
n=3; farx(n+1)=1.989018e+003; foe(n+1)=9.256317e+002; krok(n+1)=2.171980e-004; ng(n+1)=5.013401e+003;
n=4; farx(n+1)=1.635456e+003; foe(n+1)=9.035006e+002; krok(n+1)=7.102557e-004; ng(n+1)=9.534498e+002;
n=5; farx(n+1)=1.110325e+003; foe(n+1)=7.310721e+002; krok(n+1)=3.559503e-003; ng(n+1)=2.213682e+003;
n=6; farx(n+1)=9.157956e+002; foe(n+1)=6.989480e+002; krok(n+1)=9.232346e-004; ng(n+1)=1.504446e+003;
n=7; farx(n+1)=8.516994e+001; foe(n+1)=2.672040e+002; krok(n+1)=4.181157e-003; ng(n+1)=2.434021e+003;
n=8; farx(n+1)=8.379476e+001; foe(n+1)=2.650178e+002; krok(n+1)=6.070296e-006; ng(n+1)=5.084251e+003;
n=9; farx(n+1)=8.416114e+001; foe(n+1)=2.642222e+002; krok(n+1)=7.066590e-005; ng(n+1)=3.872540e+003;
n=10; farx(n+1)=8.387140e+001; foe(n+1)=2.483513e+002; krok(n+1)=1.631264e-003; ng(n+1)=3.810467e+003;
n=11; farx(n+1)=7.861450e+001; foe(n+1)=2.354117e+002; krok(n+1)=8.259539e-004; ng(n+1)=1.339502e+003;
n=12; farx(n+1)=6.652824e+001; foe(n+1)=2.202186e+002; krok(n+1)=4.800146e-004; ng(n+1)=2.696627e+003;
n=13; farx(n+1)=4.824925e+001; foe(n+1)=1.990324e+002; krok(n+1)=1.511984e-003; ng(n+1)=2.006301e+003;
n=14; farx(n+1)=4.956420e+001; foe(n+1)=1.756706e+002; krok(n+1)=1.875750e-003; ng(n+1)=3.628714e+003;
n=15; farx(n+1)=4.659882e+001; foe(n+1)=1.563630e+002; krok(n+1)=1.861026e-003; ng(n+1)=2.590637e+003;
n=16; farx(n+1)=2.007365e+001; foe(n+1)=9.720486e+001; krok(n+1)=1.005869e-003; ng(n+1)=3.212940e+003;
n=17; farx(n+1)=1.850161e+001; foe(n+1)=8.881833e+001; krok(n+1)=3.487580e-005; ng(n+1)=7.982435e+003;
n=18; farx(n+1)=1.539302e+001; foe(n+1)=7.755187e+001; krok(n+1)=1.679702e-003; ng(n+1)=4.322746e+003;
n=19; farx(n+1)=1.174127e+001; foe(n+1)=6.789512e+001; krok(n+1)=3.949694e-004; ng(n+1)=4.567790e+003;
n=20; farx(n+1)=9.134401e+000; foe(n+1)=5.984463e+001; krok(n+1)=1.524981e-003; ng(n+1)=2.658065e+003;
n=21; farx(n+1)=8.756529e+000; foe(n+1)=5.715296e+001; krok(n+1)=4.017841e-004; ng(n+1)=2.300120e+003;
n=22; farx(n+1)=9.566085e+000; foe(n+1)=5.379414e+001; krok(n+1)=1.587460e-003; ng(n+1)=2.023607e+003;
n=23; farx(n+1)=9.593683e+000; foe(n+1)=5.031221e+001; krok(n+1)=1.353531e-003; ng(n+1)=1.726244e+003;
n=24; farx(n+1)=8.695716e+000; foe(n+1)=4.641777e+001; krok(n+1)=3.671029e-003; ng(n+1)=1.251968e+003;
n=25; farx(n+1)=8.574706e+000; foe(n+1)=4.496445e+001; krok(n+1)=1.548092e-003; ng(n+1)=1.297116e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=8.572403e+000; foe(n+1)=4.450521e+001; krok(n+1)=2.991239e-006; ng(n+1)=1.752537e+003;
n=27; farx(n+1)=8.654878e+000; foe(n+1)=4.410668e+001; krok(n+1)=4.006473e-005; ng(n+1)=5.074359e+002;
n=28; farx(n+1)=8.594453e+000; foe(n+1)=4.377640e+001; krok(n+1)=2.885108e-005; ng(n+1)=5.219764e+002;
n=29; farx(n+1)=8.565450e+000; foe(n+1)=4.235657e+001; krok(n+1)=3.127269e-004; ng(n+1)=4.057461e+002;
n=30; farx(n+1)=7.619427e+000; foe(n+1)=4.091355e+001; krok(n+1)=8.336543e-004; ng(n+1)=3.360129e+002;
n=31; farx(n+1)=6.892694e+000; foe(n+1)=3.582161e+001; krok(n+1)=3.325065e-004; ng(n+1)=6.069204e+002;
n=32; farx(n+1)=6.375370e+000; foe(n+1)=2.938680e+001; krok(n+1)=1.563831e-003; ng(n+1)=2.158110e+003;
n=33; farx(n+1)=6.103927e+000; foe(n+1)=2.576664e+001; krok(n+1)=6.200466e-004; ng(n+1)=2.485376e+003;
n=34; farx(n+1)=6.139186e+000; foe(n+1)=2.301580e+001; krok(n+1)=3.303815e-003; ng(n+1)=8.941276e+002;
n=35; farx(n+1)=5.800973e+000; foe(n+1)=1.916004e+001; krok(n+1)=1.207660e-003; ng(n+1)=2.003438e+003;
n=36; farx(n+1)=5.717526e+000; foe(n+1)=1.783304e+001; krok(n+1)=1.772228e-003; ng(n+1)=2.447236e+003;
n=37; farx(n+1)=5.578349e+000; foe(n+1)=1.707575e+001; krok(n+1)=9.446471e-004; ng(n+1)=8.235394e+002;
n=38; farx(n+1)=5.153160e+000; foe(n+1)=1.459872e+001; krok(n+1)=2.586841e-003; ng(n+1)=1.917354e+003;
n=39; farx(n+1)=4.298296e+000; foe(n+1)=1.193434e+001; krok(n+1)=2.914778e-003; ng(n+1)=1.504967e+003;
n=40; farx(n+1)=4.149643e+000; foe(n+1)=1.144916e+001; krok(n+1)=1.369522e-003; ng(n+1)=8.261374e+002;
n=41; farx(n+1)=3.748515e+000; foe(n+1)=1.033264e+001; krok(n+1)=2.104259e-003; ng(n+1)=9.461208e+002;
n=42; farx(n+1)=3.738055e+000; foe(n+1)=9.704337e+000; krok(n+1)=6.752063e-003; ng(n+1)=2.808320e+002;
n=43; farx(n+1)=4.035915e+000; foe(n+1)=8.628703e+000; krok(n+1)=1.785641e-002; ng(n+1)=4.787741e+002;
n=44; farx(n+1)=4.123784e+000; foe(n+1)=8.265650e+000; krok(n+1)=1.351485e-002; ng(n+1)=9.889873e+002;
n=45; farx(n+1)=4.134537e+000; foe(n+1)=8.000815e+000; krok(n+1)=3.478135e-003; ng(n+1)=7.188187e+002;
n=46; farx(n+1)=4.242004e+000; foe(n+1)=7.546279e+000; krok(n+1)=4.648167e-003; ng(n+1)=6.911431e+002;
n=47; farx(n+1)=4.263074e+000; foe(n+1)=7.081420e+000; krok(n+1)=4.150488e-003; ng(n+1)=1.098553e+003;
n=48; farx(n+1)=4.243423e+000; foe(n+1)=6.821014e+000; krok(n+1)=1.142218e-002; ng(n+1)=3.139831e+002;
n=49; farx(n+1)=4.053340e+000; foe(n+1)=6.515420e+000; krok(n+1)=1.333847e-002; ng(n+1)=7.164474e+002;
n=50; farx(n+1)=3.501249e+000; foe(n+1)=5.742107e+000; krok(n+1)=1.895365e-002; ng(n+1)=1.376490e+003;
%odnowa zmiennej metryki
n=51; farx(n+1)=3.496456e+000; foe(n+1)=5.673616e+000; krok(n+1)=9.808732e-006; ng(n+1)=3.972651e+002;
n=52; farx(n+1)=3.479916e+000; foe(n+1)=5.589073e+000; krok(n+1)=9.504428e-006; ng(n+1)=4.404217e+002;
n=53; farx(n+1)=3.469636e+000; foe(n+1)=5.514592e+000; krok(n+1)=2.230894e-005; ng(n+1)=2.727895e+002;
n=54; farx(n+1)=3.485506e+000; foe(n+1)=5.372840e+000; krok(n+1)=3.125285e-004; ng(n+1)=9.965526e+001;
n=55; farx(n+1)=3.498967e+000; foe(n+1)=5.309469e+000; krok(n+1)=3.870864e-004; ng(n+1)=7.667637e+001;
n=56; farx(n+1)=3.510592e+000; foe(n+1)=5.300404e+000; krok(n+1)=7.616673e-005; ng(n+1)=7.451666e+001;
n=57; farx(n+1)=3.474742e+000; foe(n+1)=5.067045e+000; krok(n+1)=5.915258e-003; ng(n+1)=5.361235e+001;
n=58; farx(n+1)=3.350171e+000; foe(n+1)=4.807255e+000; krok(n+1)=5.852382e-003; ng(n+1)=2.875555e+002;
n=59; farx(n+1)=3.102177e+000; foe(n+1)=4.658979e+000; krok(n+1)=8.402186e-003; ng(n+1)=6.708573e+002;
n=60; farx(n+1)=2.898620e+000; foe(n+1)=4.513877e+000; krok(n+1)=7.242402e-003; ng(n+1)=3.777820e+002;
n=61; farx(n+1)=2.604750e+000; foe(n+1)=4.287855e+000; krok(n+1)=7.007869e-003; ng(n+1)=3.792175e+002;
n=62; farx(n+1)=2.293324e+000; foe(n+1)=4.100021e+000; krok(n+1)=5.945944e-003; ng(n+1)=3.517440e+002;
n=63; farx(n+1)=1.905094e+000; foe(n+1)=3.769722e+000; krok(n+1)=3.400118e-003; ng(n+1)=6.838991e+002;
n=64; farx(n+1)=1.788536e+000; foe(n+1)=3.659513e+000; krok(n+1)=2.162140e-003; ng(n+1)=5.511129e+002;
n=65; farx(n+1)=1.643299e+000; foe(n+1)=3.550728e+000; krok(n+1)=1.089156e-002; ng(n+1)=1.869549e+002;
n=66; farx(n+1)=1.436240e+000; foe(n+1)=3.392431e+000; krok(n+1)=1.363820e-002; ng(n+1)=4.120404e+002;
n=67; farx(n+1)=1.352094e+000; foe(n+1)=3.310749e+000; krok(n+1)=2.216378e-003; ng(n+1)=7.338612e+002;
n=68; farx(n+1)=1.244663e+000; foe(n+1)=3.169821e+000; krok(n+1)=1.881905e-002; ng(n+1)=3.230431e+002;
n=69; farx(n+1)=1.126201e+000; foe(n+1)=3.030172e+000; krok(n+1)=1.548285e-002; ng(n+1)=4.946370e+002;
n=70; farx(n+1)=1.046552e+000; foe(n+1)=2.868033e+000; krok(n+1)=2.008828e-002; ng(n+1)=4.683794e+002;
n=71; farx(n+1)=9.933023e-001; foe(n+1)=2.726707e+000; krok(n+1)=6.273397e-003; ng(n+1)=6.228364e+002;
n=72; farx(n+1)=9.777827e-001; foe(n+1)=2.652521e+000; krok(n+1)=6.035177e-003; ng(n+1)=3.524410e+002;
n=73; farx(n+1)=9.621321e-001; foe(n+1)=2.543132e+000; krok(n+1)=9.010365e-003; ng(n+1)=4.577659e+002;
n=74; farx(n+1)=9.609666e-001; foe(n+1)=2.475519e+000; krok(n+1)=1.039581e-002; ng(n+1)=6.160451e+002;
n=75; farx(n+1)=8.632605e-001; foe(n+1)=2.309971e+000; krok(n+1)=5.983917e-002; ng(n+1)=4.948036e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=8.618564e-001; foe(n+1)=2.283084e+000; krok(n+1)=1.416880e-006; ng(n+1)=6.109238e+002;
n=77; farx(n+1)=8.576723e-001; foe(n+1)=2.255649e+000; krok(n+1)=1.654527e-005; ng(n+1)=1.786960e+002;
n=78; farx(n+1)=8.572871e-001; foe(n+1)=2.253531e+000; krok(n+1)=5.807097e-006; ng(n+1)=8.677391e+001;
n=79; farx(n+1)=8.595327e-001; foe(n+1)=2.248227e+000; krok(n+1)=2.045211e-004; ng(n+1)=2.386017e+001;
n=80; farx(n+1)=8.611790e-001; foe(n+1)=2.237242e+000; krok(n+1)=1.213830e-004; ng(n+1)=3.972617e+001;
n=81; farx(n+1)=8.657947e-001; foe(n+1)=2.221270e+000; krok(n+1)=2.995066e-004; ng(n+1)=3.305864e+001;
n=82; farx(n+1)=8.681716e-001; foe(n+1)=2.182552e+000; krok(n+1)=6.180497e-004; ng(n+1)=3.730814e+001;
n=83; farx(n+1)=8.748750e-001; foe(n+1)=2.135342e+000; krok(n+1)=4.352152e-003; ng(n+1)=1.024106e+002;
n=84; farx(n+1)=8.761996e-001; foe(n+1)=2.101698e+000; krok(n+1)=4.723705e-003; ng(n+1)=2.932170e+002;
n=85; farx(n+1)=8.608356e-001; foe(n+1)=2.033153e+000; krok(n+1)=8.744518e-003; ng(n+1)=2.153635e+002;
n=86; farx(n+1)=8.484630e-001; foe(n+1)=1.945347e+000; krok(n+1)=3.128460e-003; ng(n+1)=2.643984e+002;
n=87; farx(n+1)=8.392706e-001; foe(n+1)=1.911081e+000; krok(n+1)=2.841023e-003; ng(n+1)=6.129638e+002;
n=88; farx(n+1)=8.293419e-001; foe(n+1)=1.863400e+000; krok(n+1)=8.928204e-003; ng(n+1)=2.792341e+002;
n=89; farx(n+1)=8.149717e-001; foe(n+1)=1.845258e+000; krok(n+1)=3.485331e-003; ng(n+1)=2.230708e+002;
n=90; farx(n+1)=7.967868e-001; foe(n+1)=1.823074e+000; krok(n+1)=7.481826e-003; ng(n+1)=1.531403e+002;
n=91; farx(n+1)=7.672271e-001; foe(n+1)=1.768102e+000; krok(n+1)=2.489395e-002; ng(n+1)=3.713201e+002;
n=92; farx(n+1)=7.511105e-001; foe(n+1)=1.732031e+000; krok(n+1)=7.632150e-003; ng(n+1)=3.150601e+002;
n=93; farx(n+1)=7.336886e-001; foe(n+1)=1.707546e+000; krok(n+1)=1.618451e-002; ng(n+1)=3.755994e+002;
n=94; farx(n+1)=7.244988e-001; foe(n+1)=1.684153e+000; krok(n+1)=1.834657e-002; ng(n+1)=2.516731e+002;
n=95; farx(n+1)=7.197377e-001; foe(n+1)=1.651470e+000; krok(n+1)=1.857911e-002; ng(n+1)=2.131264e+002;
n=96; farx(n+1)=7.079009e-001; foe(n+1)=1.611840e+000; krok(n+1)=1.661865e-002; ng(n+1)=4.929762e+002;
n=97; farx(n+1)=6.861729e-001; foe(n+1)=1.555416e+000; krok(n+1)=2.366211e-002; ng(n+1)=7.261568e+002;
n=98; farx(n+1)=6.647227e-001; foe(n+1)=1.521375e+000; krok(n+1)=4.739825e-002; ng(n+1)=4.198900e+002;
n=99; farx(n+1)=6.375156e-001; foe(n+1)=1.453405e+000; krok(n+1)=2.702970e-002; ng(n+1)=2.878947e+002;
n=100; farx(n+1)=5.603775e-001; foe(n+1)=1.321329e+000; krok(n+1)=5.753341e-002; ng(n+1)=3.165678e+002;
%odnowa zmiennej metryki
n=101; farx(n+1)=5.603247e-001; foe(n+1)=1.312959e+000; krok(n+1)=7.602531e-007; ng(n+1)=4.690402e+002;
n=102; farx(n+1)=5.605232e-001; foe(n+1)=1.306930e+000; krok(n+1)=2.832743e-005; ng(n+1)=7.254232e+001;
n=103; farx(n+1)=5.606969e-001; foe(n+1)=1.302873e+000; krok(n+1)=3.779992e-006; ng(n+1)=1.566934e+002;
n=104; farx(n+1)=5.598251e-001; foe(n+1)=1.297865e+000; krok(n+1)=2.139878e-005; ng(n+1)=7.048592e+001;
n=105; farx(n+1)=5.592418e-001; foe(n+1)=1.293054e+000; krok(n+1)=3.200520e-005; ng(n+1)=5.136412e+001;
n=106; farx(n+1)=5.584733e-001; foe(n+1)=1.289984e+000; krok(n+1)=1.109501e-004; ng(n+1)=2.439125e+001;
n=107; farx(n+1)=5.583079e-001; foe(n+1)=1.283798e+000; krok(n+1)=6.328553e-004; ng(n+1)=1.693101e+001;
n=108; farx(n+1)=5.626746e-001; foe(n+1)=1.273585e+000; krok(n+1)=3.253129e-003; ng(n+1)=1.081312e+001;
n=109; farx(n+1)=5.712776e-001; foe(n+1)=1.258204e+000; krok(n+1)=4.375565e-003; ng(n+1)=3.687977e+001;
n=110; farx(n+1)=5.806793e-001; foe(n+1)=1.246113e+000; krok(n+1)=2.051078e-003; ng(n+1)=1.036363e+002;
n=111; farx(n+1)=5.880073e-001; foe(n+1)=1.227790e+000; krok(n+1)=7.530795e-003; ng(n+1)=1.422615e+002;
n=112; farx(n+1)=5.834121e-001; foe(n+1)=1.201343e+000; krok(n+1)=8.685192e-003; ng(n+1)=2.137974e+002;
n=113; farx(n+1)=5.783115e-001; foe(n+1)=1.161606e+000; krok(n+1)=1.564801e-002; ng(n+1)=2.532756e+002;
n=114; farx(n+1)=5.773108e-001; foe(n+1)=1.143662e+000; krok(n+1)=3.544455e-003; ng(n+1)=2.276843e+002;
n=115; farx(n+1)=5.784051e-001; foe(n+1)=1.124213e+000; krok(n+1)=5.318174e-003; ng(n+1)=4.323285e+002;
n=116; farx(n+1)=5.755671e-001; foe(n+1)=1.116351e+000; krok(n+1)=3.271540e-003; ng(n+1)=1.981408e+002;
n=117; farx(n+1)=5.648848e-001; foe(n+1)=1.090817e+000; krok(n+1)=1.263902e-002; ng(n+1)=1.403514e+002;
n=118; farx(n+1)=5.583704e-001; foe(n+1)=1.072097e+000; krok(n+1)=6.905382e-003; ng(n+1)=6.178316e+002;
n=119; farx(n+1)=5.484939e-001; foe(n+1)=1.054175e+000; krok(n+1)=1.895365e-002; ng(n+1)=3.798240e+002;
n=120; farx(n+1)=5.401831e-001; foe(n+1)=1.037626e+000; krok(n+1)=2.390821e-002; ng(n+1)=3.635678e+002;
n=121; farx(n+1)=5.354712e-001; foe(n+1)=1.020678e+000; krok(n+1)=7.242402e-003; ng(n+1)=4.192338e+002;
n=122; farx(n+1)=5.246048e-001; foe(n+1)=9.974471e-001; krok(n+1)=5.576355e-002; ng(n+1)=3.201492e+002;
n=123; farx(n+1)=5.203267e-001; foe(n+1)=9.717560e-001; krok(n+1)=4.076532e-002; ng(n+1)=6.243116e+002;
n=124; farx(n+1)=5.159381e-001; foe(n+1)=9.401874e-001; krok(n+1)=2.243330e-002; ng(n+1)=3.576676e+002;
n=125; farx(n+1)=5.150264e-001; foe(n+1)=9.185095e-001; krok(n+1)=2.085765e-002; ng(n+1)=2.135022e+002;
%odnowa zmiennej metryki
n=126; farx(n+1)=5.149699e-001; foe(n+1)=9.092153e-001; krok(n+1)=8.284739e-007; ng(n+1)=4.385353e+002;
n=127; farx(n+1)=5.149467e-001; foe(n+1)=9.078001e-001; krok(n+1)=9.130878e-007; ng(n+1)=1.771315e+002;
n=128; farx(n+1)=5.145372e-001; foe(n+1)=9.041126e-001; krok(n+1)=2.069140e-005; ng(n+1)=6.271083e+001;
n=129; farx(n+1)=5.145611e-001; foe(n+1)=9.023522e-001; krok(n+1)=1.373458e-005; ng(n+1)=5.566010e+001;
n=130; farx(n+1)=5.144396e-001; foe(n+1)=8.990273e-001; krok(n+1)=4.508384e-005; ng(n+1)=3.747417e+001;
n=131; farx(n+1)=5.142475e-001; foe(n+1)=8.931741e-001; krok(n+1)=3.087231e-004; ng(n+1)=2.374485e+001;
n=132; farx(n+1)=5.146045e-001; foe(n+1)=8.838022e-001; krok(n+1)=2.375655e-004; ng(n+1)=2.845210e+001;
n=133; farx(n+1)=5.143351e-001; foe(n+1)=8.765836e-001; krok(n+1)=1.143304e-003; ng(n+1)=2.073726e+001;
n=134; farx(n+1)=5.154968e-001; foe(n+1)=8.704916e-001; krok(n+1)=1.478815e-003; ng(n+1)=3.263551e+001;
n=135; farx(n+1)=5.206928e-001; foe(n+1)=8.575802e-001; krok(n+1)=2.000182e-002; ng(n+1)=7.262277e+001;
n=136; farx(n+1)=5.231005e-001; foe(n+1)=8.422925e-001; krok(n+1)=6.390530e-003; ng(n+1)=1.303216e+002;
n=137; farx(n+1)=5.207922e-001; foe(n+1)=8.375484e-001; krok(n+1)=2.957629e-003; ng(n+1)=2.189540e+002;
n=138; farx(n+1)=5.139808e-001; foe(n+1)=8.287230e-001; krok(n+1)=3.426502e-003; ng(n+1)=1.985234e+002;
n=139; farx(n+1)=5.129228e-001; foe(n+1)=8.201163e-001; krok(n+1)=1.405854e-002; ng(n+1)=2.723022e+002;
n=140; farx(n+1)=5.138078e-001; foe(n+1)=8.093674e-001; krok(n+1)=2.038894e-002; ng(n+1)=3.257925e+002;
n=141; farx(n+1)=5.174100e-001; foe(n+1)=8.018671e-001; krok(n+1)=8.202462e-003; ng(n+1)=1.466112e+002;
n=142; farx(n+1)=5.164649e-001; foe(n+1)=7.920553e-001; krok(n+1)=1.525498e-002; ng(n+1)=3.131558e+002;
n=143; farx(n+1)=5.146374e-001; foe(n+1)=7.794559e-001; krok(n+1)=1.450364e-002; ng(n+1)=2.124485e+002;
n=144; farx(n+1)=5.156545e-001; foe(n+1)=7.739996e-001; krok(n+1)=9.778245e-003; ng(n+1)=2.010392e+002;
n=145; farx(n+1)=5.226205e-001; foe(n+1)=7.670464e-001; krok(n+1)=1.830545e-002; ng(n+1)=3.621932e+002;
n=146; farx(n+1)=5.203319e-001; foe(n+1)=7.619436e-001; krok(n+1)=3.552033e-002; ng(n+1)=1.399297e+002;
n=147; farx(n+1)=5.104624e-001; foe(n+1)=7.411130e-001; krok(n+1)=5.113139e-002; ng(n+1)=3.422119e+002;
n=148; farx(n+1)=5.078905e-001; foe(n+1)=7.258221e-001; krok(n+1)=2.490870e-002; ng(n+1)=3.553308e+002;
n=149; farx(n+1)=5.082324e-001; foe(n+1)=7.108598e-001; krok(n+1)=3.847317e-002; ng(n+1)=1.759242e+002;
n=150; farx(n+1)=5.066858e-001; foe(n+1)=6.957185e-001; krok(n+1)=1.432736e-002; ng(n+1)=4.067700e+002;
%odnowa zmiennej metryki
n=151; farx(n+1)=5.066669e-001; foe(n+1)=6.917447e-001; krok(n+1)=6.395081e-007; ng(n+1)=3.353304e+002;
n=152; farx(n+1)=5.065860e-001; foe(n+1)=6.905494e-001; krok(n+1)=1.457486e-006; ng(n+1)=1.235375e+002;
n=153; farx(n+1)=5.064663e-001; foe(n+1)=6.882219e-001; krok(n+1)=3.838059e-006; ng(n+1)=1.373298e+002;
n=154; farx(n+1)=5.065342e-001; foe(n+1)=6.843712e-001; krok(n+1)=3.034095e-005; ng(n+1)=6.024834e+001;
n=155; farx(n+1)=5.065063e-001; foe(n+1)=6.806463e-001; krok(n+1)=4.229784e-005; ng(n+1)=4.317483e+001;
n=156; farx(n+1)=5.066020e-001; foe(n+1)=6.635317e-001; krok(n+1)=4.078161e-004; ng(n+1)=3.306661e+001;
n=157; farx(n+1)=5.067984e-001; foe(n+1)=6.603232e-001; krok(n+1)=9.960454e-005; ng(n+1)=2.698571e+001;
n=158; farx(n+1)=5.065425e-001; foe(n+1)=6.552636e-001; krok(n+1)=2.722891e-003; ng(n+1)=1.133612e+001;
n=159; farx(n+1)=5.063149e-001; foe(n+1)=6.537509e-001; krok(n+1)=7.010405e-004; ng(n+1)=1.850264e+001;
n=160; farx(n+1)=5.065937e-001; foe(n+1)=6.516059e-001; krok(n+1)=2.209886e-003; ng(n+1)=2.935515e+001;
n=161; farx(n+1)=5.065077e-001; foe(n+1)=6.338235e-001; krok(n+1)=8.705227e-003; ng(n+1)=4.981089e+001;
n=162; farx(n+1)=5.029965e-001; foe(n+1)=6.251712e-001; krok(n+1)=1.501192e-002; ng(n+1)=4.604037e+002;
n=163; farx(n+1)=5.011074e-001; foe(n+1)=6.192501e-001; krok(n+1)=3.416397e-003; ng(n+1)=2.400027e+002;
n=164; farx(n+1)=5.028271e-001; foe(n+1)=6.137409e-001; krok(n+1)=1.170476e-002; ng(n+1)=1.767109e+002;
n=165; farx(n+1)=5.044898e-001; foe(n+1)=6.014453e-001; krok(n+1)=3.201356e-002; ng(n+1)=1.964416e+002;
n=166; farx(n+1)=5.035662e-001; foe(n+1)=5.928847e-001; krok(n+1)=7.191611e-003; ng(n+1)=1.564910e+002;
n=167; farx(n+1)=5.015834e-001; foe(n+1)=5.879612e-001; krok(n+1)=1.664870e-002; ng(n+1)=1.654279e+002;
n=168; farx(n+1)=5.028460e-001; foe(n+1)=5.851897e-001; krok(n+1)=2.831403e-003; ng(n+1)=2.389793e+002;
n=169; farx(n+1)=5.086598e-001; foe(n+1)=5.769634e-001; krok(n+1)=2.977641e-002; ng(n+1)=1.430036e+002;
n=170; farx(n+1)=5.115711e-001; foe(n+1)=5.696668e-001; krok(n+1)=3.104078e-002; ng(n+1)=2.890347e+002;
n=171; farx(n+1)=5.218811e-001; foe(n+1)=5.583274e-001; krok(n+1)=2.946581e-002; ng(n+1)=4.205908e+002;
n=172; farx(n+1)=5.222597e-001; foe(n+1)=5.464204e-001; krok(n+1)=2.857938e-002; ng(n+1)=1.485964e+002;
n=173; farx(n+1)=5.206898e-001; foe(n+1)=5.417965e-001; krok(n+1)=1.171709e-002; ng(n+1)=2.411459e+002;
n=174; farx(n+1)=5.195412e-001; foe(n+1)=5.374877e-001; krok(n+1)=8.099144e-003; ng(n+1)=2.585593e+002;
n=175; farx(n+1)=5.201083e-001; foe(n+1)=5.252907e-001; krok(n+1)=3.899500e-002; ng(n+1)=2.250688e+002;
%odnowa zmiennej metryki
n=176; farx(n+1)=5.201049e-001; foe(n+1)=5.250214e-001; krok(n+1)=1.122542e-006; ng(n+1)=7.777306e+001;
n=177; farx(n+1)=5.200649e-001; foe(n+1)=5.240722e-001; krok(n+1)=7.194757e-006; ng(n+1)=5.261283e+001;
n=178; farx(n+1)=5.200263e-001; foe(n+1)=5.232964e-001; krok(n+1)=7.611413e-007; ng(n+1)=1.466660e+002;
n=179; farx(n+1)=5.198927e-001; foe(n+1)=5.223361e-001; krok(n+1)=1.386876e-005; ng(n+1)=4.034545e+001;
n=180; farx(n+1)=5.198135e-001; foe(n+1)=5.207380e-001; krok(n+1)=1.602043e-005; ng(n+1)=4.821811e+001;
n=181; farx(n+1)=5.186449e-001; foe(n+1)=5.159302e-001; krok(n+1)=5.252174e-004; ng(n+1)=1.675175e+001;
n=182; farx(n+1)=5.185082e-001; foe(n+1)=5.155813e-001; krok(n+1)=5.788477e-005; ng(n+1)=1.392645e+001;
n=183; farx(n+1)=5.159494e-001; foe(n+1)=5.110508e-001; krok(n+1)=4.564205e-003; ng(n+1)=7.232403e+000;
n=184; farx(n+1)=5.153082e-001; foe(n+1)=5.091153e-001; krok(n+1)=5.747417e-004; ng(n+1)=2.873619e+001;
n=185; farx(n+1)=5.149409e-001; foe(n+1)=5.065151e-001; krok(n+1)=2.914397e-003; ng(n+1)=4.907739e+001;
n=186; farx(n+1)=5.147654e-001; foe(n+1)=5.034315e-001; krok(n+1)=2.487486e-003; ng(n+1)=1.022770e+002;
n=187; farx(n+1)=5.174482e-001; foe(n+1)=4.993417e-001; krok(n+1)=7.805422e-003; ng(n+1)=2.069533e+002;
n=188; farx(n+1)=5.167743e-001; foe(n+1)=4.958266e-001; krok(n+1)=1.355378e-002; ng(n+1)=4.573611e+002;
n=189; farx(n+1)=5.154676e-001; foe(n+1)=4.935537e-001; krok(n+1)=3.411474e-003; ng(n+1)=3.660946e+002;
n=190; farx(n+1)=5.182503e-001; foe(n+1)=4.901710e-001; krok(n+1)=2.959126e-002; ng(n+1)=3.186733e+002;
n=191; farx(n+1)=5.162667e-001; foe(n+1)=4.838931e-001; krok(n+1)=1.382426e-002; ng(n+1)=2.519428e+002;
n=192; farx(n+1)=5.137437e-001; foe(n+1)=4.799851e-001; krok(n+1)=1.095618e-002; ng(n+1)=3.146737e+002;
n=193; farx(n+1)=5.111316e-001; foe(n+1)=4.752770e-001; krok(n+1)=1.660195e-002; ng(n+1)=4.226838e+002;
n=194; farx(n+1)=5.067716e-001; foe(n+1)=4.641145e-001; krok(n+1)=3.388472e-002; ng(n+1)=1.606684e+002;
n=195; farx(n+1)=5.086171e-001; foe(n+1)=4.561087e-001; krok(n+1)=2.812841e-002; ng(n+1)=2.184440e+002;
n=196; farx(n+1)=5.075257e-001; foe(n+1)=4.465225e-001; krok(n+1)=2.643052e-002; ng(n+1)=1.678396e+002;
n=197; farx(n+1)=5.061845e-001; foe(n+1)=4.358017e-001; krok(n+1)=1.076348e-002; ng(n+1)=5.124162e+002;
n=198; farx(n+1)=5.035283e-001; foe(n+1)=4.302397e-001; krok(n+1)=1.647425e-002; ng(n+1)=1.663686e+002;
n=199; farx(n+1)=5.020874e-001; foe(n+1)=4.253495e-001; krok(n+1)=2.277959e-002; ng(n+1)=2.284511e+002;
n=200; farx(n+1)=4.989388e-001; foe(n+1)=4.147051e-001; krok(n+1)=2.941687e-002; ng(n+1)=2.811863e+002;
%odnowa zmiennej metryki
n=201; farx(n+1)=4.989242e-001; foe(n+1)=4.096102e-001; krok(n+1)=3.479426e-007; ng(n+1)=5.842290e+002;
n=202; farx(n+1)=4.988413e-001; foe(n+1)=4.076133e-001; krok(n+1)=3.744041e-006; ng(n+1)=1.095572e+002;
n=203; farx(n+1)=4.988313e-001; foe(n+1)=4.073102e-001; krok(n+1)=9.412859e-007; ng(n+1)=9.214737e+001;
n=204; farx(n+1)=4.988177e-001; foe(n+1)=4.025817e-001; krok(n+1)=4.765565e-005; ng(n+1)=5.125740e+001;
n=205; farx(n+1)=4.987803e-001; foe(n+1)=4.018840e-001; krok(n+1)=1.726474e-005; ng(n+1)=3.561741e+001;
n=206; farx(n+1)=4.988002e-001; foe(n+1)=3.959477e-001; krok(n+1)=3.296469e-004; ng(n+1)=2.415865e+001;
n=207; farx(n+1)=4.988360e-001; foe(n+1)=3.945475e-001; krok(n+1)=8.923578e-005; ng(n+1)=2.095476e+001;
n=208; farx(n+1)=4.990015e-001; foe(n+1)=3.881617e-001; krok(n+1)=4.927090e-004; ng(n+1)=2.145271e+001;
n=209; farx(n+1)=5.000211e-001; foe(n+1)=3.858392e-001; krok(n+1)=9.053002e-004; ng(n+1)=1.857936e+001;
n=210; farx(n+1)=5.026990e-001; foe(n+1)=3.826357e-001; krok(n+1)=5.548179e-003; ng(n+1)=2.242646e+001;
n=211; farx(n+1)=5.040277e-001; foe(n+1)=3.798750e-001; krok(n+1)=5.207429e-003; ng(n+1)=1.951898e+001;
n=212; farx(n+1)=5.050761e-001; foe(n+1)=3.773187e-001; krok(n+1)=1.704444e-003; ng(n+1)=4.502564e+001;
n=213; farx(n+1)=5.048563e-001; foe(n+1)=3.735901e-001; krok(n+1)=8.516524e-003; ng(n+1)=1.361278e+002;
n=214; farx(n+1)=5.053642e-001; foe(n+1)=3.700437e-001; krok(n+1)=2.643052e-002; ng(n+1)=3.523379e+002;
n=215; farx(n+1)=5.079841e-001; foe(n+1)=3.663288e-001; krok(n+1)=1.321526e-002; ng(n+1)=8.146598e+001;
n=216; farx(n+1)=5.074658e-001; foe(n+1)=3.610702e-001; krok(n+1)=2.896961e-002; ng(n+1)=2.006168e+002;
n=217; farx(n+1)=5.073046e-001; foe(n+1)=3.555727e-001; krok(n+1)=1.199493e-002; ng(n+1)=3.561440e+002;
n=218; farx(n+1)=5.086376e-001; foe(n+1)=3.532709e-001; krok(n+1)=4.653602e-003; ng(n+1)=2.710836e+002;
n=219; farx(n+1)=5.080478e-001; foe(n+1)=3.488970e-001; krok(n+1)=1.586357e-002; ng(n+1)=3.406120e+002;
n=220; farx(n+1)=5.070858e-001; foe(n+1)=3.444482e-001; krok(n+1)=2.680031e-002; ng(n+1)=1.162207e+002;
n=221; farx(n+1)=5.075430e-001; foe(n+1)=3.414617e-001; krok(n+1)=1.296098e-002; ng(n+1)=1.446103e+002;
n=222; farx(n+1)=5.090573e-001; foe(n+1)=3.395386e-001; krok(n+1)=1.153579e-002; ng(n+1)=1.590361e+002;
n=223; farx(n+1)=5.110057e-001; foe(n+1)=3.370337e-001; krok(n+1)=2.440235e-002; ng(n+1)=1.616891e+002;
n=224; farx(n+1)=5.124063e-001; foe(n+1)=3.303560e-001; krok(n+1)=3.023343e-002; ng(n+1)=1.456342e+002;
n=225; farx(n+1)=5.098952e-001; foe(n+1)=3.228864e-001; krok(n+1)=5.121364e-002; ng(n+1)=2.737476e+002;
%odnowa zmiennej metryki
n=226; farx(n+1)=5.099046e-001; foe(n+1)=3.224215e-001; krok(n+1)=3.588312e-007; ng(n+1)=1.794068e+002;
n=227; farx(n+1)=5.099082e-001; foe(n+1)=3.217600e-001; krok(n+1)=1.779653e-006; ng(n+1)=8.936572e+001;
n=228; farx(n+1)=5.099134e-001; foe(n+1)=3.214075e-001; krok(n+1)=1.001277e-006; ng(n+1)=1.001143e+002;
n=229; farx(n+1)=5.098658e-001; foe(n+1)=3.188974e-001; krok(n+1)=5.493831e-005; ng(n+1)=3.114257e+001;
n=230; farx(n+1)=5.098303e-001; foe(n+1)=3.178803e-001; krok(n+1)=1.891271e-005; ng(n+1)=3.621746e+001;
n=231; farx(n+1)=5.097409e-001; foe(n+1)=3.144521e-001; krok(n+1)=9.199795e-005; ng(n+1)=3.454599e+001;
n=232; farx(n+1)=5.096755e-001; foe(n+1)=3.136200e-001; krok(n+1)=1.022356e-004; ng(n+1)=1.647407e+001;
n=233; farx(n+1)=5.090950e-001; foe(n+1)=3.106862e-001; krok(n+1)=5.823131e-004; ng(n+1)=1.173366e+001;
n=234; farx(n+1)=5.087429e-001; foe(n+1)=3.080585e-001; krok(n+1)=1.291963e-003; ng(n+1)=1.057513e+001;
n=235; farx(n+1)=5.081584e-001; foe(n+1)=3.073641e-001; krok(n+1)=1.350179e-003; ng(n+1)=2.912644e+001;
n=236; farx(n+1)=5.074750e-001; foe(n+1)=3.064259e-001; krok(n+1)=2.804162e-003; ng(n+1)=3.852924e+001;
n=237; farx(n+1)=5.064953e-001; foe(n+1)=3.052239e-001; krok(n+1)=5.946996e-003; ng(n+1)=5.457645e+001;
n=238; farx(n+1)=5.076061e-001; foe(n+1)=3.039765e-001; krok(n+1)=8.894860e-003; ng(n+1)=1.026390e+002;
n=239; farx(n+1)=5.069782e-001; foe(n+1)=3.033090e-001; krok(n+1)=1.121665e-002; ng(n+1)=1.983597e+002;
n=240; farx(n+1)=5.062727e-001; foe(n+1)=3.011527e-001; krok(n+1)=3.277333e-002; ng(n+1)=2.582331e+002;
n=241; farx(n+1)=5.063546e-001; foe(n+1)=2.995452e-001; krok(n+1)=1.547432e-002; ng(n+1)=3.276720e+002;
n=242; farx(n+1)=5.050372e-001; foe(n+1)=2.956944e-001; krok(n+1)=5.016986e-002; ng(n+1)=2.490822e+002;
n=243; farx(n+1)=5.055119e-001; foe(n+1)=2.923223e-001; krok(n+1)=1.501192e-002; ng(n+1)=1.132026e+002;
n=244; farx(n+1)=5.048017e-001; foe(n+1)=2.909656e-001; krok(n+1)=1.064021e-002; ng(n+1)=2.665297e+002;
n=245; farx(n+1)=5.027379e-001; foe(n+1)=2.891067e-001; krok(n+1)=2.398987e-002; ng(n+1)=1.432175e+002;
n=246; farx(n+1)=5.025539e-001; foe(n+1)=2.878568e-001; krok(n+1)=1.650507e-002; ng(n+1)=1.649524e+002;
n=247; farx(n+1)=5.017493e-001; foe(n+1)=2.849793e-001; krok(n+1)=6.754317e-002; ng(n+1)=8.833947e+001;
n=248; farx(n+1)=5.017761e-001; foe(n+1)=2.829679e-001; krok(n+1)=1.540867e-002; ng(n+1)=1.644796e+002;
n=249; farx(n+1)=5.042260e-001; foe(n+1)=2.764110e-001; krok(n+1)=4.614317e-002; ng(n+1)=1.670336e+002;
n=250; farx(n+1)=5.054311e-001; foe(n+1)=2.730160e-001; krok(n+1)=1.878731e-002; ng(n+1)=3.007422e+002;
%odnowa zmiennej metryki
n=251; farx(n+1)=5.054348e-001; foe(n+1)=2.723550e-001; krok(n+1)=2.125074e-006; ng(n+1)=9.168504e+001;
n=252; farx(n+1)=5.054361e-001; foe(n+1)=2.717422e-001; krok(n+1)=6.042203e-007; ng(n+1)=1.603351e+002;
n=253; farx(n+1)=5.054418e-001; foe(n+1)=2.714037e-001; krok(n+1)=6.150701e-007; ng(n+1)=1.044063e+002;
n=254; farx(n+1)=5.054426e-001; foe(n+1)=2.711059e-001; krok(n+1)=3.873148e-006; ng(n+1)=3.648785e+001;
n=255; farx(n+1)=5.055927e-001; foe(n+1)=2.700832e-001; krok(n+1)=2.247399e-004; ng(n+1)=1.236000e+001;
n=256; farx(n+1)=5.056875e-001; foe(n+1)=2.696419e-001; krok(n+1)=7.681907e-005; ng(n+1)=1.230451e+001;
n=257; farx(n+1)=5.058808e-001; foe(n+1)=2.689591e-001; krok(n+1)=2.501060e-004; ng(n+1)=9.278166e+000;
n=258; farx(n+1)=5.062627e-001; foe(n+1)=2.676589e-001; krok(n+1)=1.711903e-004; ng(n+1)=1.340504e+001;
n=259; farx(n+1)=5.063353e-001; foe(n+1)=2.670450e-001; krok(n+1)=9.208067e-004; ng(n+1)=6.234136e+000;
n=260; farx(n+1)=5.064843e-001; foe(n+1)=2.667528e-001; krok(n+1)=1.420511e-003; ng(n+1)=7.013267e+000;
n=261; farx(n+1)=5.065926e-001; foe(n+1)=2.661507e-001; krok(n+1)=2.464960e-003; ng(n+1)=1.054909e+001;
n=262; farx(n+1)=5.064168e-001; foe(n+1)=2.649534e-001; krok(n+1)=5.548179e-003; ng(n+1)=2.196949e+001;
n=263; farx(n+1)=5.070592e-001; foe(n+1)=2.642298e-001; krok(n+1)=1.884952e-002; ng(n+1)=7.433623e+001;
n=264; farx(n+1)=5.070042e-001; foe(n+1)=2.636814e-001; krok(n+1)=1.000726e-002; ng(n+1)=1.410993e+002;
n=265; farx(n+1)=5.085863e-001; foe(n+1)=2.628418e-001; krok(n+1)=1.450364e-002; ng(n+1)=1.609029e+002;
n=266; farx(n+1)=5.098029e-001; foe(n+1)=2.622095e-001; krok(n+1)=1.251384e-002; ng(n+1)=1.160929e+002;
n=267; farx(n+1)=5.100755e-001; foe(n+1)=2.607905e-001; krok(n+1)=4.034388e-002; ng(n+1)=1.280481e+002;
n=268; farx(n+1)=5.114783e-001; foe(n+1)=2.596105e-001; krok(n+1)=1.019133e-002; ng(n+1)=2.108174e+002;
n=269; farx(n+1)=5.120662e-001; foe(n+1)=2.588465e-001; krok(n+1)=1.296098e-002; ng(n+1)=1.857875e+002;
n=270; farx(n+1)=5.121823e-001; foe(n+1)=2.580119e-001; krok(n+1)=1.586441e-002; ng(n+1)=7.486274e+001;
n=271; farx(n+1)=5.144908e-001; foe(n+1)=2.558613e-001; krok(n+1)=9.139037e-002; ng(n+1)=1.609721e+002;
n=272; farx(n+1)=5.153969e-001; foe(n+1)=2.547189e-001; krok(n+1)=1.936662e-002; ng(n+1)=1.191431e+002;
n=273; farx(n+1)=5.183586e-001; foe(n+1)=2.520627e-001; krok(n+1)=4.773765e-002; ng(n+1)=2.839688e+002;
n=274; farx(n+1)=5.190174e-001; foe(n+1)=2.506979e-001; krok(n+1)=2.757245e-002; ng(n+1)=1.566395e+002;
n=275; farx(n+1)=5.204797e-001; foe(n+1)=2.488464e-001; krok(n+1)=3.552033e-002; ng(n+1)=1.129536e+002;
%odnowa zmiennej metryki
n=276; farx(n+1)=5.204811e-001; foe(n+1)=2.485212e-001; krok(n+1)=2.293407e-007; ng(n+1)=1.767627e+002;
n=277; farx(n+1)=5.204810e-001; foe(n+1)=2.484932e-001; krok(n+1)=4.165909e-007; ng(n+1)=3.701572e+001;
n=278; farx(n+1)=5.204815e-001; foe(n+1)=2.484770e-001; krok(n+1)=8.445859e-006; ng(n+1)=8.335489e+000;
n=279; farx(n+1)=5.204738e-001; foe(n+1)=2.483967e-001; krok(n+1)=1.944453e-005; ng(n+1)=1.056244e+001;
n=280; farx(n+1)=5.204598e-001; foe(n+1)=2.483103e-001; krok(n+1)=1.641052e-005; ng(n+1)=1.223688e+001;
n=281; farx(n+1)=5.204209e-001; foe(n+1)=2.481526e-001; krok(n+1)=9.874235e-005; ng(n+1)=6.165263e+000;
n=282; farx(n+1)=5.203388e-001; foe(n+1)=2.480067e-001; krok(n+1)=1.823889e-004; ng(n+1)=4.426502e+000;
n=283; farx(n+1)=5.200911e-001; foe(n+1)=2.476485e-001; krok(n+1)=2.154098e-004; ng(n+1)=6.769031e+000;
n=284; farx(n+1)=5.196859e-001; foe(n+1)=2.472945e-001; krok(n+1)=1.234612e-003; ng(n+1)=3.912992e+000;
n=285; farx(n+1)=5.190688e-001; foe(n+1)=2.469562e-001; krok(n+1)=1.427772e-003; ng(n+1)=6.958335e+000;
n=286; farx(n+1)=5.186902e-001; foe(n+1)=2.466035e-001; krok(n+1)=1.457334e-003; ng(n+1)=1.110488e+001;
n=287; farx(n+1)=5.187290e-001; foe(n+1)=2.461259e-001; krok(n+1)=8.003391e-003; ng(n+1)=1.698152e+001;
n=288; farx(n+1)=5.193321e-001; foe(n+1)=2.449699e-001; krok(n+1)=1.321323e-002; ng(n+1)=1.736593e+001;
n=289; farx(n+1)=5.191991e-001; foe(n+1)=2.444300e-001; krok(n+1)=6.606617e-003; ng(n+1)=5.162087e+001;
n=290; farx(n+1)=5.196483e-001; foe(n+1)=2.440632e-001; krok(n+1)=6.353015e-003; ng(n+1)=8.182993e+001;
n=291; farx(n+1)=5.194562e-001; foe(n+1)=2.434971e-001; krok(n+1)=2.173419e-002; ng(n+1)=8.896534e+001;
n=292; farx(n+1)=5.180054e-001; foe(n+1)=2.426045e-001; krok(n+1)=5.234463e-002; ng(n+1)=8.567086e+001;
n=293; farx(n+1)=5.172982e-001; foe(n+1)=2.417148e-001; krok(n+1)=4.388300e-002; ng(n+1)=1.064826e+002;
n=294; farx(n+1)=5.165787e-001; foe(n+1)=2.409967e-001; krok(n+1)=4.331299e-002; ng(n+1)=3.823394e+001;
n=295; farx(n+1)=5.160559e-001; foe(n+1)=2.406855e-001; krok(n+1)=1.129495e-002; ng(n+1)=1.328543e+002;
n=296; farx(n+1)=5.157907e-001; foe(n+1)=2.402665e-001; krok(n+1)=1.501192e-002; ng(n+1)=7.082982e+001;
n=297; farx(n+1)=5.146281e-001; foe(n+1)=2.387066e-001; krok(n+1)=4.130221e-002; ng(n+1)=8.294024e+001;
n=298; farx(n+1)=5.139376e-001; foe(n+1)=2.369341e-001; krok(n+1)=5.381807e-002; ng(n+1)=7.769603e+001;
n=299; farx(n+1)=5.139877e-001; foe(n+1)=2.358154e-001; krok(n+1)=4.523999e-002; ng(n+1)=1.905643e+002;
n=300; farx(n+1)=5.135710e-001; foe(n+1)=2.343582e-001; krok(n+1)=5.716001e-002; ng(n+1)=7.540582e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
