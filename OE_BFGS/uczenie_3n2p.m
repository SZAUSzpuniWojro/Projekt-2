%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.665853e+003; foe(n+1)=4.677540e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=4.024334e+003; foe(n+1)=3.987108e+003; krok(n+1)=5.113081e-004; ng(n+1)=1.813951e+003;
n=2; farx(n+1)=9.905605e+002; foe(n+1)=8.394121e+002; krok(n+1)=7.558357e-003; ng(n+1)=7.248714e+002;
n=3; farx(n+1)=1.043583e+003; foe(n+1)=7.735954e+002; krok(n+1)=1.274265e-004; ng(n+1)=1.906503e+003;
n=4; farx(n+1)=1.197922e+003; foe(n+1)=7.586368e+002; krok(n+1)=7.432430e-004; ng(n+1)=1.506579e+003;
n=5; farx(n+1)=8.583979e+002; foe(n+1)=5.996686e+002; krok(n+1)=2.011815e-003; ng(n+1)=7.039803e+002;
n=6; farx(n+1)=5.941880e+002; foe(n+1)=5.201514e+002; krok(n+1)=1.111417e-003; ng(n+1)=1.211380e+003;
n=7; farx(n+1)=4.286016e+002; foe(n+1)=4.788254e+002; krok(n+1)=5.126539e-004; ng(n+1)=2.098539e+003;
n=8; farx(n+1)=3.451916e+002; foe(n+1)=4.598329e+002; krok(n+1)=9.671447e-004; ng(n+1)=1.790612e+003;
n=9; farx(n+1)=2.653030e+002; foe(n+1)=4.284017e+002; krok(n+1)=1.442683e-003; ng(n+1)=1.056182e+003;
n=10; farx(n+1)=1.590490e+002; foe(n+1)=3.924425e+002; krok(n+1)=3.049962e-003; ng(n+1)=1.642995e+003;
n=11; farx(n+1)=1.413460e+002; foe(n+1)=3.857971e+002; krok(n+1)=7.046696e-004; ng(n+1)=2.620436e+003;
n=12; farx(n+1)=1.338564e+002; foe(n+1)=3.825189e+002; krok(n+1)=2.859446e-004; ng(n+1)=3.506549e+003;
n=13; farx(n+1)=1.197352e+002; foe(n+1)=3.622053e+002; krok(n+1)=1.184956e-002; ng(n+1)=3.743054e+003;
n=14; farx(n+1)=1.161094e+002; foe(n+1)=3.527766e+002; krok(n+1)=5.643882e-004; ng(n+1)=3.638174e+003;
n=15; farx(n+1)=1.276446e+002; foe(n+1)=3.406179e+002; krok(n+1)=2.316999e-003; ng(n+1)=5.911106e+003;
n=16; farx(n+1)=1.310570e+002; foe(n+1)=3.309508e+002; krok(n+1)=1.634950e-002; ng(n+1)=3.517696e+003;
n=17; farx(n+1)=1.459242e+002; foe(n+1)=2.992679e+002; krok(n+1)=4.111242e-002; ng(n+1)=4.392970e+003;
n=18; farx(n+1)=1.158236e+002; foe(n+1)=2.568386e+002; krok(n+1)=1.210679e-001; ng(n+1)=3.016792e+003;
n=19; farx(n+1)=7.370049e+001; foe(n+1)=2.208149e+002; krok(n+1)=6.024636e-002; ng(n+1)=2.689463e+003;
n=20; farx(n+1)=4.768321e+001; foe(n+1)=1.698591e+002; krok(n+1)=1.486329e-001; ng(n+1)=4.949771e+003;
n=21; farx(n+1)=3.814432e+001; foe(n+1)=1.510106e+002; krok(n+1)=1.067078e-001; ng(n+1)=8.928882e+003;
n=22; farx(n+1)=2.753693e+001; foe(n+1)=1.182198e+002; krok(n+1)=3.655615e-001; ng(n+1)=5.913041e+003;
n=23; farx(n+1)=2.537488e+001; foe(n+1)=9.180074e+001; krok(n+1)=3.319551e-001; ng(n+1)=8.500054e+003;
n=24; farx(n+1)=1.862906e+001; foe(n+1)=8.116001e+001; krok(n+1)=3.550835e-001; ng(n+1)=1.536681e+003;
n=25; farx(n+1)=1.766652e+001; foe(n+1)=7.738685e+001; krok(n+1)=3.005969e-001; ng(n+1)=1.727090e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=1.765234e+001; foe(n+1)=7.733687e+001; krok(n+1)=1.357816e-006; ng(n+1)=5.542490e+002;
n=27; farx(n+1)=1.798336e+001; foe(n+1)=7.689811e+001; krok(n+1)=1.219587e-005; ng(n+1)=4.963476e+002;
n=28; farx(n+1)=1.792861e+001; foe(n+1)=7.635609e+001; krok(n+1)=1.906434e-004; ng(n+1)=1.553186e+002;
n=29; farx(n+1)=1.718072e+001; foe(n+1)=7.521123e+001; krok(n+1)=4.966995e-004; ng(n+1)=1.349958e+002;
n=30; farx(n+1)=1.776512e+001; foe(n+1)=7.368937e+001; krok(n+1)=8.566255e-004; ng(n+1)=1.441859e+002;
n=31; farx(n+1)=1.740360e+001; foe(n+1)=7.222165e+001; krok(n+1)=8.522871e-004; ng(n+1)=1.543377e+003;
n=32; farx(n+1)=1.411540e+001; foe(n+1)=6.746674e+001; krok(n+1)=6.607631e-003; ng(n+1)=3.882985e+003;
n=33; farx(n+1)=1.306376e+001; foe(n+1)=6.606065e+001; krok(n+1)=3.766475e-004; ng(n+1)=3.652945e+003;
n=34; farx(n+1)=1.274462e+001; foe(n+1)=6.429713e+001; krok(n+1)=3.584297e-003; ng(n+1)=1.589440e+003;
n=35; farx(n+1)=1.176875e+001; foe(n+1)=6.310795e+001; krok(n+1)=2.286608e-003; ng(n+1)=1.979199e+003;
n=36; farx(n+1)=1.254863e+001; foe(n+1)=6.199352e+001; krok(n+1)=2.316852e-002; ng(n+1)=4.085911e+003;
n=37; farx(n+1)=1.227289e+001; foe(n+1)=5.869999e+001; krok(n+1)=2.972664e-002; ng(n+1)=4.162191e+003;
n=38; farx(n+1)=1.116309e+001; foe(n+1)=5.687693e+001; krok(n+1)=3.280985e-002; ng(n+1)=4.713917e+002;
n=39; farx(n+1)=1.029303e+001; foe(n+1)=5.615384e+001; krok(n+1)=2.067140e-002; ng(n+1)=9.861467e+002;
n=40; farx(n+1)=9.154090e+000; foe(n+1)=5.409963e+001; krok(n+1)=5.155493e-001; ng(n+1)=1.715951e+003;
n=41; farx(n+1)=8.084888e+000; foe(n+1)=5.301195e+001; krok(n+1)=4.331299e-002; ng(n+1)=2.531740e+003;
n=42; farx(n+1)=7.729590e+000; foe(n+1)=5.002479e+001; krok(n+1)=5.714745e-001; ng(n+1)=3.174880e+003;
n=43; farx(n+1)=7.983332e+000; foe(n+1)=4.917986e+001; krok(n+1)=1.794664e-001; ng(n+1)=2.850447e+003;
n=44; farx(n+1)=8.592271e+000; foe(n+1)=4.605225e+001; krok(n+1)=3.485301e-001; ng(n+1)=2.481209e+003;
n=45; farx(n+1)=8.648841e+000; foe(n+1)=4.559015e+001; krok(n+1)=3.245113e-002; ng(n+1)=1.223298e+003;
n=46; farx(n+1)=8.239803e+000; foe(n+1)=4.219214e+001; krok(n+1)=2.935451e-001; ng(n+1)=1.490923e+003;
n=47; farx(n+1)=8.340962e+000; foe(n+1)=4.128890e+001; krok(n+1)=1.383770e-001; ng(n+1)=1.711782e+003;
n=48; farx(n+1)=8.358040e+000; foe(n+1)=4.109039e+001; krok(n+1)=1.091588e-002; ng(n+1)=1.282220e+003;
n=49; farx(n+1)=7.924788e+000; foe(n+1)=3.962596e+001; krok(n+1)=1.664220e-001; ng(n+1)=1.391421e+003;
n=50; farx(n+1)=7.819109e+000; foe(n+1)=3.911494e+001; krok(n+1)=7.689657e-002; ng(n+1)=5.945150e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=7.832558e+000; foe(n+1)=3.896602e+001; krok(n+1)=7.645508e-007; ng(n+1)=1.279309e+003;
n=52; farx(n+1)=7.834253e+000; foe(n+1)=3.859723e+001; krok(n+1)=1.052881e-005; ng(n+1)=5.940285e+002;
n=53; farx(n+1)=7.779140e+000; foe(n+1)=3.852200e+001; krok(n+1)=1.308609e-005; ng(n+1)=2.043604e+002;
n=54; farx(n+1)=7.436137e+000; foe(n+1)=3.828011e+001; krok(n+1)=4.494798e-004; ng(n+1)=9.221828e+001;
n=55; farx(n+1)=7.587974e+000; foe(n+1)=3.807180e+001; krok(n+1)=1.563635e-004; ng(n+1)=1.354843e+002;
n=56; farx(n+1)=7.723219e+000; foe(n+1)=3.795802e+001; krok(n+1)=2.247378e-004; ng(n+1)=1.370003e+002;
n=57; farx(n+1)=7.749747e+000; foe(n+1)=3.783860e+001; krok(n+1)=4.526501e-004; ng(n+1)=1.340542e+002;
n=58; farx(n+1)=7.561232e+000; foe(n+1)=3.748223e+001; krok(n+1)=7.476124e-004; ng(n+1)=2.038633e+002;
n=59; farx(n+1)=7.506869e+000; foe(n+1)=3.702838e+001; krok(n+1)=7.883344e-003; ng(n+1)=1.082074e+003;
n=60; farx(n+1)=7.308026e+000; foe(n+1)=3.648386e+001; krok(n+1)=9.809382e-003; ng(n+1)=1.669143e+003;
n=61; farx(n+1)=7.307161e+000; foe(n+1)=3.615903e+001; krok(n+1)=4.403244e-003; ng(n+1)=2.941776e+003;
n=62; farx(n+1)=7.360956e+000; foe(n+1)=3.540794e+001; krok(n+1)=1.783149e-002; ng(n+1)=2.343835e+003;
n=63; farx(n+1)=7.354374e+000; foe(n+1)=3.537579e+001; krok(n+1)=1.524981e-003; ng(n+1)=3.592153e+003;
n=64; farx(n+1)=7.387994e+000; foe(n+1)=3.502051e+001; krok(n+1)=2.075244e-003; ng(n+1)=3.710809e+003;
n=65; farx(n+1)=7.511323e+000; foe(n+1)=3.470496e+001; krok(n+1)=1.750226e-002; ng(n+1)=3.084215e+003;
n=66; farx(n+1)=7.517110e+000; foe(n+1)=3.459455e+001; krok(n+1)=1.640492e-002; ng(n+1)=3.590949e+003;
n=67; farx(n+1)=7.468111e+000; foe(n+1)=3.390925e+001; krok(n+1)=1.940886e-001; ng(n+1)=2.766745e+003;
n=68; farx(n+1)=7.244395e+000; foe(n+1)=3.263711e+001; krok(n+1)=5.954236e-002; ng(n+1)=3.506536e+003;
n=69; farx(n+1)=7.179813e+000; foe(n+1)=3.188692e+001; krok(n+1)=2.601807e-001; ng(n+1)=1.324727e+003;
n=70; farx(n+1)=7.146622e+000; foe(n+1)=3.163751e+001; krok(n+1)=1.124683e-001; ng(n+1)=3.384797e+002;
n=71; farx(n+1)=7.050548e+000; foe(n+1)=3.123359e+001; krok(n+1)=6.995615e-002; ng(n+1)=8.000550e+002;
n=72; farx(n+1)=7.102045e+000; foe(n+1)=3.096868e+001; krok(n+1)=9.137743e-002; ng(n+1)=6.053009e+002;
n=73; farx(n+1)=7.191997e+000; foe(n+1)=3.051114e+001; krok(n+1)=4.228884e-001; ng(n+1)=5.967631e+002;
n=74; farx(n+1)=7.343332e+000; foe(n+1)=3.020289e+001; krok(n+1)=1.056270e+000; ng(n+1)=9.202192e+002;
n=75; farx(n+1)=7.479831e+000; foe(n+1)=3.011726e+001; krok(n+1)=4.936410e-001; ng(n+1)=5.701399e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=7.479693e+000; foe(n+1)=3.011494e+001; krok(n+1)=3.793929e-007; ng(n+1)=2.401606e+002;
n=77; farx(n+1)=7.480815e+000; foe(n+1)=3.010580e+001; krok(n+1)=1.771248e-006; ng(n+1)=2.196797e+002;
n=78; farx(n+1)=7.501151e+000; foe(n+1)=3.010049e+001; krok(n+1)=1.416371e-005; ng(n+1)=6.222483e+001;
n=79; farx(n+1)=7.504628e+000; foe(n+1)=3.008676e+001; krok(n+1)=3.315087e-005; ng(n+1)=7.474291e+001;
n=80; farx(n+1)=7.487435e+000; foe(n+1)=3.008403e+001; krok(n+1)=4.506169e-005; ng(n+1)=3.003312e+001;
n=81; farx(n+1)=7.503654e+000; foe(n+1)=3.008091e+001; krok(n+1)=1.287051e-004; ng(n+1)=2.026783e+001;
n=82; farx(n+1)=7.495552e+000; foe(n+1)=3.007696e+001; krok(n+1)=6.946998e-004; ng(n+1)=1.376101e+001;
n=83; farx(n+1)=7.429968e+000; foe(n+1)=3.005868e+001; krok(n+1)=5.485375e-003; ng(n+1)=1.084751e+001;
n=84; farx(n+1)=7.406013e+000; foe(n+1)=3.005294e+001; krok(n+1)=2.156451e-003; ng(n+1)=7.354237e+000;
n=85; farx(n+1)=7.426169e+000; foe(n+1)=3.003551e+001; krok(n+1)=2.811708e-002; ng(n+1)=1.906934e+001;
n=86; farx(n+1)=7.431505e+000; foe(n+1)=3.003030e+001; krok(n+1)=3.657312e-003; ng(n+1)=1.558643e+002;
n=87; farx(n+1)=7.412049e+000; foe(n+1)=3.002161e+001; krok(n+1)=1.488820e-002; ng(n+1)=2.298238e+002;
n=88; farx(n+1)=7.332981e+000; foe(n+1)=3.001223e+001; krok(n+1)=1.346946e-002; ng(n+1)=3.251693e+002;
n=89; farx(n+1)=7.290424e+000; foe(n+1)=2.993496e+001; krok(n+1)=1.101515e-001; ng(n+1)=3.856813e+002;
n=90; farx(n+1)=7.289687e+000; foe(n+1)=2.992712e+001; krok(n+1)=1.897240e-002; ng(n+1)=1.830830e+002;
n=91; farx(n+1)=7.078191e+000; foe(n+1)=2.988935e+001; krok(n+1)=1.031518e-001; ng(n+1)=3.735198e+002;
n=92; farx(n+1)=6.986404e+000; foe(n+1)=2.985238e+001; krok(n+1)=6.589702e-002; ng(n+1)=5.166992e+002;
n=93; farx(n+1)=7.044505e+000; foe(n+1)=2.960603e+001; krok(n+1)=7.491049e-001; ng(n+1)=5.642459e+002;
n=94; farx(n+1)=7.050191e+000; foe(n+1)=2.959853e+001; krok(n+1)=5.214412e-003; ng(n+1)=8.505766e+002;
n=95; farx(n+1)=7.222274e+000; foe(n+1)=2.949763e+001; krok(n+1)=2.022243e-001; ng(n+1)=9.192392e+002;
n=96; farx(n+1)=7.192152e+000; foe(n+1)=2.917684e+001; krok(n+1)=9.730681e-001; ng(n+1)=9.484404e+002;
n=97; farx(n+1)=7.130654e+000; foe(n+1)=2.906902e+001; krok(n+1)=3.005969e-001; ng(n+1)=5.535903e+002;
n=98; farx(n+1)=7.057321e+000; foe(n+1)=2.893467e+001; krok(n+1)=4.429176e-001; ng(n+1)=7.673667e+002;
n=99; farx(n+1)=7.096327e+000; foe(n+1)=2.890447e+001; krok(n+1)=5.695627e-001; ng(n+1)=4.326376e+002;
n=100; farx(n+1)=7.081931e+000; foe(n+1)=2.887906e+001; krok(n+1)=6.762233e-001; ng(n+1)=5.713888e+002;
%odnowa zmiennej metryki
n=101; farx(n+1)=7.081161e+000; foe(n+1)=2.887427e+001; krok(n+1)=1.708080e-007; ng(n+1)=4.715889e+002;
n=102; farx(n+1)=7.075515e+000; foe(n+1)=2.887198e+001; krok(n+1)=4.273013e-006; ng(n+1)=6.659110e+001;
n=103; farx(n+1)=7.084227e+000; foe(n+1)=2.887095e+001; krok(n+1)=1.225575e-005; ng(n+1)=3.031153e+001;
n=104; farx(n+1)=7.108570e+000; foe(n+1)=2.886843e+001; krok(n+1)=5.685733e-006; ng(n+1)=6.824129e+001;
n=105; farx(n+1)=7.108499e+000; foe(n+1)=2.886595e+001; krok(n+1)=5.019773e-005; ng(n+1)=2.342667e+001;
n=106; farx(n+1)=7.106659e+000; foe(n+1)=2.886282e+001; krok(n+1)=9.975011e-004; ng(n+1)=6.887674e+000;
n=107; farx(n+1)=7.111055e+000; foe(n+1)=2.886036e+001; krok(n+1)=1.340733e-004; ng(n+1)=1.768654e+001;
n=108; farx(n+1)=7.100232e+000; foe(n+1)=2.885496e+001; krok(n+1)=9.584745e-004; ng(n+1)=1.251539e+001;
n=109; farx(n+1)=7.089179e+000; foe(n+1)=2.883471e+001; krok(n+1)=4.708571e-003; ng(n+1)=1.020472e+001;
n=110; farx(n+1)=7.136519e+000; foe(n+1)=2.882022e+001; krok(n+1)=8.346822e-003; ng(n+1)=1.011209e+002;
n=111; farx(n+1)=7.154929e+000; foe(n+1)=2.881640e+001; krok(n+1)=4.952046e-003; ng(n+1)=3.514365e+002;
n=112; farx(n+1)=7.142387e+000; foe(n+1)=2.880241e+001; krok(n+1)=2.366211e-002; ng(n+1)=4.634241e+002;
n=113; farx(n+1)=7.181804e+000; foe(n+1)=2.877619e+001; krok(n+1)=3.595477e-002; ng(n+1)=8.035677e+002;
n=114; farx(n+1)=7.247129e+000; foe(n+1)=2.876638e+001; krok(n+1)=1.617035e-002; ng(n+1)=6.041848e+002;
n=115; farx(n+1)=7.263274e+000; foe(n+1)=2.875538e+001; krok(n+1)=1.522640e-002; ng(n+1)=4.682184e+002;
n=116; farx(n+1)=7.276697e+000; foe(n+1)=2.874368e+001; krok(n+1)=3.955518e-002; ng(n+1)=4.978608e+002;
n=117; farx(n+1)=7.294626e+000; foe(n+1)=2.870899e+001; krok(n+1)=1.677703e-001; ng(n+1)=7.327075e+002;
n=118; farx(n+1)=7.343933e+000; foe(n+1)=2.868220e+001; krok(n+1)=2.012417e-001; ng(n+1)=6.192477e+002;
n=119; farx(n+1)=7.407868e+000; foe(n+1)=2.864996e+001; krok(n+1)=4.747686e-001; ng(n+1)=5.079424e+002;
n=120; farx(n+1)=7.436548e+000; foe(n+1)=2.863322e+001; krok(n+1)=6.744582e-001; ng(n+1)=2.390721e+002;
n=121; farx(n+1)=7.453077e+000; foe(n+1)=2.862117e+001; krok(n+1)=4.321320e-001; ng(n+1)=8.140160e+002;
n=122; farx(n+1)=7.486228e+000; foe(n+1)=2.861138e+001; krok(n+1)=3.556863e-001; ng(n+1)=5.012216e+002;
n=123; farx(n+1)=7.546608e+000; foe(n+1)=2.859679e+001; krok(n+1)=7.178655e-001; ng(n+1)=3.829338e+002;
n=124; farx(n+1)=7.663105e+000; foe(n+1)=2.857459e+001; krok(n+1)=7.294970e-001; ng(n+1)=5.519848e+002;
n=125; farx(n+1)=7.689948e+000; foe(n+1)=2.853897e+001; krok(n+1)=1.065517e+000; ng(n+1)=5.610659e+002;
%odnowa zmiennej metryki
n=126; farx(n+1)=7.690298e+000; foe(n+1)=2.853818e+001; krok(n+1)=1.494425e-007; ng(n+1)=2.180646e+002;
n=127; farx(n+1)=7.691019e+000; foe(n+1)=2.853791e+001; krok(n+1)=1.067446e-006; ng(n+1)=4.622695e+001;
n=128; farx(n+1)=7.691355e+000; foe(n+1)=2.853683e+001; krok(n+1)=9.309975e-006; ng(n+1)=3.334312e+001;
n=129; farx(n+1)=7.656282e+000; foe(n+1)=2.853236e+001; krok(n+1)=2.254192e-005; ng(n+1)=3.680917e+001;
n=130; farx(n+1)=7.661102e+000; foe(n+1)=2.853200e+001; krok(n+1)=2.444110e-005; ng(n+1)=1.356660e+001;
n=131; farx(n+1)=7.655292e+000; foe(n+1)=2.853136e+001; krok(n+1)=1.191612e-004; ng(n+1)=8.414530e+000;
n=132; farx(n+1)=7.639266e+000; foe(n+1)=2.852626e+001; krok(n+1)=1.088038e-003; ng(n+1)=8.340317e+000;
n=133; farx(n+1)=7.611074e+000; foe(n+1)=2.851670e+001; krok(n+1)=1.088038e-003; ng(n+1)=1.183030e+001;
n=134; farx(n+1)=7.511355e+000; foe(n+1)=2.848806e+001; krok(n+1)=7.859971e-003; ng(n+1)=3.805421e+001;
n=135; farx(n+1)=7.539668e+000; foe(n+1)=2.847912e+001; krok(n+1)=5.414124e-003; ng(n+1)=4.518923e+002;
n=136; farx(n+1)=7.551266e+000; foe(n+1)=2.847674e+001; krok(n+1)=2.100547e-003; ng(n+1)=4.776242e+002;
n=137; farx(n+1)=7.558245e+000; foe(n+1)=2.847115e+001; krok(n+1)=1.877194e-002; ng(n+1)=4.686531e+002;
n=138; farx(n+1)=7.503372e+000; foe(n+1)=2.846259e+001; krok(n+1)=1.600678e-002; ng(n+1)=4.580781e+002;
n=139; farx(n+1)=7.481163e+000; foe(n+1)=2.844374e+001; krok(n+1)=5.266571e-002; ng(n+1)=3.492548e+002;
n=140; farx(n+1)=7.450346e+000; foe(n+1)=2.843484e+001; krok(n+1)=2.788482e-002; ng(n+1)=2.306772e+002;
n=141; farx(n+1)=7.419163e+000; foe(n+1)=2.843010e+001; krok(n+1)=4.502187e-002; ng(n+1)=2.288276e+002;
n=142; farx(n+1)=7.363298e+000; foe(n+1)=2.841679e+001; krok(n+1)=4.938403e-001; ng(n+1)=1.813491e+002;
n=143; farx(n+1)=7.282366e+000; foe(n+1)=2.840194e+001; krok(n+1)=4.164005e-001; ng(n+1)=2.249503e+002;
n=144; farx(n+1)=7.234547e+000; foe(n+1)=2.839227e+001; krok(n+1)=4.126070e-001; ng(n+1)=9.987203e+001;
n=145; farx(n+1)=7.202745e+000; foe(n+1)=2.838977e+001; krok(n+1)=5.636152e-001; ng(n+1)=1.310574e+002;
n=146; farx(n+1)=7.180525e+000; foe(n+1)=2.838811e+001; krok(n+1)=1.045307e+000; ng(n+1)=5.424203e+001;
n=147; farx(n+1)=7.160855e+000; foe(n+1)=2.838650e+001; krok(n+1)=2.909207e+000; ng(n+1)=7.649594e+001;
n=148; farx(n+1)=7.155223e+000; foe(n+1)=2.838216e+001; krok(n+1)=5.159582e+000; ng(n+1)=1.222630e+002;
n=149; farx(n+1)=7.135654e+000; foe(n+1)=2.837348e+001; krok(n+1)=2.472812e+000; ng(n+1)=9.258503e+001;
n=150; farx(n+1)=7.108548e+000; foe(n+1)=2.836925e+001; krok(n+1)=7.386536e-001; ng(n+1)=1.858538e+002;
%odnowa zmiennej metryki
n=151; farx(n+1)=7.108448e+000; foe(n+1)=2.836807e+001; krok(n+1)=1.938318e-007; ng(n+1)=2.364945e+002;
n=152; farx(n+1)=7.110455e+000; foe(n+1)=2.836563e+001; krok(n+1)=1.079046e-006; ng(n+1)=1.384670e+002;
n=153; farx(n+1)=7.114677e+000; foe(n+1)=2.836531e+001; krok(n+1)=4.886358e-006; ng(n+1)=2.604671e+001;
n=154; farx(n+1)=7.118043e+000; foe(n+1)=2.836473e+001; krok(n+1)=2.282707e-005; ng(n+1)=1.613862e+001;
n=155; farx(n+1)=7.112317e+000; foe(n+1)=2.836448e+001; krok(n+1)=2.921834e-005; ng(n+1)=9.120367e+000;
n=156; farx(n+1)=7.117152e+000; foe(n+1)=2.836428e+001; krok(n+1)=1.910401e-004; ng(n+1)=3.904708e+000;
n=157; farx(n+1)=7.112432e+000; foe(n+1)=2.836317e+001; krok(n+1)=9.540187e-004; ng(n+1)=4.129239e+000;
n=158; farx(n+1)=7.124494e+000; foe(n+1)=2.836011e+001; krok(n+1)=1.463096e-003; ng(n+1)=5.116919e+000;
n=159; farx(n+1)=7.145031e+000; foe(n+1)=2.835065e+001; krok(n+1)=4.008930e-003; ng(n+1)=2.083677e+001;
n=160; farx(n+1)=7.157804e+000; foe(n+1)=2.834617e+001; krok(n+1)=5.214412e-003; ng(n+1)=2.272789e+002;
n=161; farx(n+1)=7.153166e+000; foe(n+1)=2.834512e+001; krok(n+1)=3.966104e-003; ng(n+1)=4.157810e+002;
n=162; farx(n+1)=7.139815e+000; foe(n+1)=2.834122e+001; krok(n+1)=3.911298e-002; ng(n+1)=4.368799e+002;
n=163; farx(n+1)=7.206229e+000; foe(n+1)=2.833446e+001; krok(n+1)=1.815661e-002; ng(n+1)=6.026735e+002;
n=164; farx(n+1)=7.247420e+000; foe(n+1)=2.831347e+001; krok(n+1)=6.170195e-002; ng(n+1)=7.476012e+002;
n=165; farx(n+1)=7.273694e+000; foe(n+1)=2.830511e+001; krok(n+1)=3.172883e-002; ng(n+1)=7.565133e+002;
n=166; farx(n+1)=7.300399e+000; foe(n+1)=2.829862e+001; krok(n+1)=2.096657e-002; ng(n+1)=7.382507e+002;
n=167; farx(n+1)=7.299543e+000; foe(n+1)=2.828383e+001; krok(n+1)=1.093143e-001; ng(n+1)=6.524405e+002;
n=168; farx(n+1)=7.441909e+000; foe(n+1)=2.825331e+001; krok(n+1)=7.832156e-001; ng(n+1)=5.858900e+002;
n=169; farx(n+1)=7.485168e+000; foe(n+1)=2.823970e+001; krok(n+1)=5.788951e-001; ng(n+1)=2.074692e+002;
n=170; farx(n+1)=7.507902e+000; foe(n+1)=2.821810e+001; krok(n+1)=6.293018e-001; ng(n+1)=6.402346e+002;
n=171; farx(n+1)=7.532689e+000; foe(n+1)=2.820698e+001; krok(n+1)=4.063898e-001; ng(n+1)=7.634683e+002;
n=172; farx(n+1)=7.577100e+000; foe(n+1)=2.819555e+001; krok(n+1)=3.427774e-001; ng(n+1)=2.471349e+002;
n=173; farx(n+1)=7.625084e+000; foe(n+1)=2.817614e+001; krok(n+1)=7.571530e-001; ng(n+1)=2.222341e+002;
n=174; farx(n+1)=7.637746e+000; foe(n+1)=2.816967e+001; krok(n+1)=3.550835e-001; ng(n+1)=6.941484e+002;
n=175; farx(n+1)=7.666878e+000; foe(n+1)=2.816116e+001; krok(n+1)=5.403833e-001; ng(n+1)=3.384441e+002;
%odnowa zmiennej metryki
n=176; farx(n+1)=7.667470e+000; foe(n+1)=2.815868e+001; krok(n+1)=1.244400e-007; ng(n+1)=4.058005e+002;
n=177; farx(n+1)=7.669501e+000; foe(n+1)=2.815687e+001; krok(n+1)=9.998608e-007; ng(n+1)=1.341702e+002;
n=178; farx(n+1)=7.678116e+000; foe(n+1)=2.815579e+001; krok(n+1)=2.207467e-006; ng(n+1)=6.775756e+001;
n=179; farx(n+1)=7.684484e+000; foe(n+1)=2.815541e+001; krok(n+1)=1.321254e-005; ng(n+1)=1.783247e+001;
n=180; farx(n+1)=7.677358e+000; foe(n+1)=2.815447e+001; krok(n+1)=1.784716e-004; ng(n+1)=7.542847e+000;
n=181; farx(n+1)=7.673702e+000; foe(n+1)=2.815388e+001; krok(n+1)=1.464903e-004; ng(n+1)=6.900505e+000;
n=182; farx(n+1)=7.681275e+000; foe(n+1)=2.815351e+001; krok(n+1)=1.767236e-004; ng(n+1)=5.956993e+000;
n=183; farx(n+1)=7.689051e+000; foe(n+1)=2.814496e+001; krok(n+1)=3.473809e-003; ng(n+1)=6.623888e+000;
n=184; farx(n+1)=7.689714e+000; foe(n+1)=2.813433e+001; krok(n+1)=3.910756e-003; ng(n+1)=2.648525e+001;
n=185; farx(n+1)=7.676409e+000; foe(n+1)=2.813077e+001; krok(n+1)=2.849598e-003; ng(n+1)=1.579989e+002;
n=186; farx(n+1)=7.642090e+000; foe(n+1)=2.812878e+001; krok(n+1)=5.173682e-003; ng(n+1)=2.109145e+002;
n=187; farx(n+1)=7.610880e+000; foe(n+1)=2.812540e+001; krok(n+1)=3.389322e-002; ng(n+1)=2.173033e+002;
n=188; farx(n+1)=7.619067e+000; foe(n+1)=2.812318e+001; krok(n+1)=7.505960e-003; ng(n+1)=2.240506e+002;
n=189; farx(n+1)=7.577188e+000; foe(n+1)=2.811625e+001; krok(n+1)=6.899709e-002; ng(n+1)=1.397116e+002;
n=190; farx(n+1)=7.545089e+000; foe(n+1)=2.810494e+001; krok(n+1)=4.438543e-002; ng(n+1)=1.656815e+002;
n=191; farx(n+1)=7.508270e+000; foe(n+1)=2.810166e+001; krok(n+1)=5.801457e-002; ng(n+1)=1.610264e+002;
n=192; farx(n+1)=7.415087e+000; foe(n+1)=2.808983e+001; krok(n+1)=3.307424e-001; ng(n+1)=1.400661e+002;
n=193; farx(n+1)=7.326541e+000; foe(n+1)=2.807885e+001; krok(n+1)=4.635137e-001; ng(n+1)=1.160388e+002;
n=194; farx(n+1)=7.294219e+000; foe(n+1)=2.807304e+001; krok(n+1)=3.103642e-001; ng(n+1)=1.382057e+002;
n=195; farx(n+1)=7.262688e+000; foe(n+1)=2.806785e+001; krok(n+1)=9.775646e-001; ng(n+1)=1.473469e+002;
n=196; farx(n+1)=7.196050e+000; foe(n+1)=2.806320e+001; krok(n+1)=1.601772e+000; ng(n+1)=5.105897e+001;
n=197; farx(n+1)=7.211953e+000; foe(n+1)=2.805854e+001; krok(n+1)=1.914853e+000; ng(n+1)=5.985344e+001;
n=198; farx(n+1)=7.194080e+000; foe(n+1)=2.805544e+001; krok(n+1)=1.189931e+000; ng(n+1)=1.886377e+002;
n=199; farx(n+1)=7.189797e+000; foe(n+1)=2.805390e+001; krok(n+1)=4.228884e-001; ng(n+1)=3.989762e+001;
n=200; farx(n+1)=7.214368e+000; foe(n+1)=2.805133e+001; krok(n+1)=1.033731e+000; ng(n+1)=7.020124e+001;
%odnowa zmiennej metryki
n=201; farx(n+1)=7.214147e+000; foe(n+1)=2.805093e+001; krok(n+1)=1.062719e-007; ng(n+1)=1.750927e+002;
n=202; farx(n+1)=7.213661e+000; foe(n+1)=2.805090e+001; krok(n+1)=2.233942e-006; ng(n+1)=1.157941e+001;
n=203; farx(n+1)=7.211309e+000; foe(n+1)=2.805083e+001; krok(n+1)=7.792977e-006; ng(n+1)=9.994969e+000;
n=204; farx(n+1)=7.209819e+000; foe(n+1)=2.805081e+001; krok(n+1)=4.789737e-006; ng(n+1)=6.156887e+000;
n=205; farx(n+1)=7.210889e+000; foe(n+1)=2.805071e+001; krok(n+1)=4.705554e-005; ng(n+1)=5.124106e+000;
n=206; farx(n+1)=7.214323e+000; foe(n+1)=2.805056e+001; krok(n+1)=1.714180e-004; ng(n+1)=3.077417e+000;
n=207; farx(n+1)=7.210116e+000; foe(n+1)=2.805025e+001; krok(n+1)=2.752027e-004; ng(n+1)=3.520779e+000;
n=208; farx(n+1)=7.217380e+000; foe(n+1)=2.804932e+001; krok(n+1)=1.186919e-003; ng(n+1)=3.026460e+000;
n=209; farx(n+1)=7.223018e+000; foe(n+1)=2.804766e+001; krok(n+1)=4.038010e-003; ng(n+1)=2.283838e+000;
n=210; farx(n+1)=7.229252e+000; foe(n+1)=2.804503e+001; krok(n+1)=6.480489e-003; ng(n+1)=2.036977e+001;
n=211; farx(n+1)=7.236664e+000; foe(n+1)=2.804337e+001; krok(n+1)=4.473296e-003; ng(n+1)=9.931361e+001;
n=212; farx(n+1)=7.221634e+000; foe(n+1)=2.804190e+001; krok(n+1)=3.426151e-002; ng(n+1)=1.778210e+002;
n=213; farx(n+1)=7.271053e+000; foe(n+1)=2.803835e+001; krok(n+1)=2.054980e-002; ng(n+1)=2.682588e+002;
n=214; farx(n+1)=7.296459e+000; foe(n+1)=2.803341e+001; krok(n+1)=5.983917e-002; ng(n+1)=4.672767e+002;
n=215; farx(n+1)=7.326483e+000; foe(n+1)=2.802212e+001; krok(n+1)=6.252611e-002; ng(n+1)=5.253899e+002;
n=216; farx(n+1)=7.360105e+000; foe(n+1)=2.801799e+001; krok(n+1)=4.438543e-002; ng(n+1)=4.537112e+002;
n=217; farx(n+1)=7.388659e+000; foe(n+1)=2.800549e+001; krok(n+1)=1.471339e-001; ng(n+1)=3.398925e+002;
n=218; farx(n+1)=7.437905e+000; foe(n+1)=2.799133e+001; krok(n+1)=4.352152e-001; ng(n+1)=4.297808e+002;
n=219; farx(n+1)=7.507131e+000; foe(n+1)=2.797655e+001; krok(n+1)=1.535939e+000; ng(n+1)=2.378973e+002;
n=220; farx(n+1)=7.525186e+000; foe(n+1)=2.796622e+001; krok(n+1)=4.572701e-001; ng(n+1)=5.283843e+002;
n=221; farx(n+1)=7.567386e+000; foe(n+1)=2.795898e+001; krok(n+1)=5.422916e-001; ng(n+1)=1.763618e+002;
n=222; farx(n+1)=7.603427e+000; foe(n+1)=2.795374e+001; krok(n+1)=3.923690e-001; ng(n+1)=2.967922e+002;
n=223; farx(n+1)=7.634433e+000; foe(n+1)=2.794961e+001; krok(n+1)=1.300618e-001; ng(n+1)=3.899200e+002;
n=224; farx(n+1)=7.654162e+000; foe(n+1)=2.794634e+001; krok(n+1)=3.160606e-001; ng(n+1)=9.097236e+001;
n=225; farx(n+1)=7.702023e+000; foe(n+1)=2.793718e+001; krok(n+1)=1.289896e+000; ng(n+1)=3.500934e+002;
%odnowa zmiennej metryki
n=226; farx(n+1)=7.702523e+000; foe(n+1)=2.793582e+001; krok(n+1)=9.529678e-008; ng(n+1)=3.530161e+002;
n=227; farx(n+1)=7.705202e+000; foe(n+1)=2.793465e+001; krok(n+1)=1.110414e-006; ng(n+1)=9.735570e+001;
n=228; farx(n+1)=7.712655e+000; foe(n+1)=2.793391e+001; krok(n+1)=1.623756e-006; ng(n+1)=6.276620e+001;
n=229; farx(n+1)=7.718078e+000; foe(n+1)=2.793347e+001; krok(n+1)=8.600086e-006; ng(n+1)=2.174693e+001;
n=230; farx(n+1)=7.708069e+000; foe(n+1)=2.793303e+001; krok(n+1)=7.257464e-005; ng(n+1)=7.182031e+000;
n=231; farx(n+1)=7.709567e+000; foe(n+1)=2.793293e+001; krok(n+1)=1.562643e-004; ng(n+1)=2.991415e+000;
n=232; farx(n+1)=7.708576e+000; foe(n+1)=2.793062e+001; krok(n+1)=1.830478e-003; ng(n+1)=3.918812e+000;
n=233; farx(n+1)=7.714425e+000; foe(n+1)=2.792812e+001; krok(n+1)=6.767654e-004; ng(n+1)=9.563503e+000;
n=234; farx(n+1)=7.704910e+000; foe(n+1)=2.792253e+001; krok(n+1)=2.774090e-003; ng(n+1)=2.807009e+001;
n=235; farx(n+1)=7.690168e+000; foe(n+1)=2.792058e+001; krok(n+1)=4.213325e-003; ng(n+1)=1.131410e+002;
n=236; farx(n+1)=7.681819e+000; foe(n+1)=2.792004e+001; krok(n+1)=3.966104e-003; ng(n+1)=1.165222e+002;
n=237; farx(n+1)=7.648613e+000; foe(n+1)=2.791856e+001; krok(n+1)=1.865271e-002; ng(n+1)=1.172033e+002;
n=238; farx(n+1)=7.658841e+000; foe(n+1)=2.791716e+001; krok(n+1)=9.454744e-003; ng(n+1)=1.305640e+002;
n=239; farx(n+1)=7.633958e+000; foe(n+1)=2.791486e+001; krok(n+1)=5.857529e-002; ng(n+1)=1.113838e+002;
n=240; farx(n+1)=7.565222e+000; foe(n+1)=2.789832e+001; krok(n+1)=2.144024e-001; ng(n+1)=1.590430e+002;
n=241; farx(n+1)=7.515668e+000; foe(n+1)=2.789392e+001; krok(n+1)=4.560093e-002; ng(n+1)=2.157607e+002;
n=242; farx(n+1)=7.429647e+000; foe(n+1)=2.788372e+001; krok(n+1)=3.839848e-001; ng(n+1)=1.629928e+002;
n=243; farx(n+1)=7.362501e+000; foe(n+1)=2.787699e+001; krok(n+1)=3.240219e-001; ng(n+1)=1.047633e+002;
n=244; farx(n+1)=7.322531e+000; foe(n+1)=2.787027e+001; krok(n+1)=6.954939e-001; ng(n+1)=2.465144e+002;
n=245; farx(n+1)=7.294502e+000; foe(n+1)=2.786550e+001; krok(n+1)=5.824576e-001; ng(n+1)=1.672805e+002;
n=246; farx(n+1)=7.234695e+000; foe(n+1)=2.785949e+001; krok(n+1)=1.727068e+000; ng(n+1)=4.375073e+001;
n=247; farx(n+1)=7.243124e+000; foe(n+1)=2.785692e+001; krok(n+1)=9.372047e-001; ng(n+1)=1.201877e+002;
n=248; farx(n+1)=7.229896e+000; foe(n+1)=2.785562e+001; krok(n+1)=7.113726e-001; ng(n+1)=1.591994e+002;
n=249; farx(n+1)=7.231877e+000; foe(n+1)=2.785478e+001; krok(n+1)=3.963765e-001; ng(n+1)=7.727353e+001;
n=250; farx(n+1)=7.247324e+000; foe(n+1)=2.785347e+001; krok(n+1)=1.164353e+000; ng(n+1)=6.760785e+001;
%odnowa zmiennej metryki
n=251; farx(n+1)=7.247305e+000; foe(n+1)=2.785344e+001; krok(n+1)=9.384149e-008; ng(n+1)=5.015920e+001;
n=252; farx(n+1)=7.247813e+000; foe(n+1)=2.785340e+001; krok(n+1)=1.302477e-006; ng(n+1)=1.809697e+001;
n=253; farx(n+1)=7.249598e+000; foe(n+1)=2.785335e+001; krok(n+1)=4.725500e-006; ng(n+1)=8.528933e+000;
n=254; farx(n+1)=7.249635e+000; foe(n+1)=2.785329e+001; krok(n+1)=6.370154e-006; ng(n+1)=9.187824e+000;
n=255; farx(n+1)=7.246302e+000; foe(n+1)=2.785324e+001; krok(n+1)=1.977673e-005; ng(n+1)=4.965715e+000;
n=256; farx(n+1)=7.247197e+000; foe(n+1)=2.785294e+001; krok(n+1)=1.786250e-003; ng(n+1)=1.344383e+000;
n=257; farx(n+1)=7.245901e+000; foe(n+1)=2.785264e+001; krok(n+1)=4.129769e-004; ng(n+1)=2.834431e+000;
n=258; farx(n+1)=7.251901e+000; foe(n+1)=2.785230e+001; krok(n+1)=2.135248e-004; ng(n+1)=4.489742e+000;
n=259; farx(n+1)=7.254819e+000; foe(n+1)=2.785108e+001; krok(n+1)=3.113587e-003; ng(n+1)=3.305682e+000;
n=260; farx(n+1)=7.255294e+000; foe(n+1)=2.784965e+001; krok(n+1)=5.433546e-003; ng(n+1)=2.577168e+001;
n=261; farx(n+1)=7.257768e+000; foe(n+1)=2.784916e+001; krok(n+1)=3.544455e-003; ng(n+1)=7.404834e+001;
n=262; farx(n+1)=7.251351e+000; foe(n+1)=2.784805e+001; krok(n+1)=4.842911e-002; ng(n+1)=9.751283e+001;
n=263; farx(n+1)=7.290192e+000; foe(n+1)=2.784664e+001; krok(n+1)=1.333847e-002; ng(n+1)=1.850176e+002;
n=264; farx(n+1)=7.304359e+000; foe(n+1)=2.784505e+001; krok(n+1)=3.243320e-002; ng(n+1)=2.967409e+002;
n=265; farx(n+1)=7.342240e+000; foe(n+1)=2.783861e+001; krok(n+1)=2.373843e-001; ng(n+1)=3.817162e+002;
n=266; farx(n+1)=7.337226e+000; foe(n+1)=2.783593e+001; krok(n+1)=4.048730e-002; ng(n+1)=2.156076e+002;
n=267; farx(n+1)=7.399552e+000; foe(n+1)=2.782495e+001; krok(n+1)=2.857373e-001; ng(n+1)=2.731211e+002;
n=268; farx(n+1)=7.439394e+000; foe(n+1)=2.781808e+001; krok(n+1)=3.655097e-001; ng(n+1)=2.943472e+002;
n=269; farx(n+1)=7.467164e+000; foe(n+1)=2.780698e+001; krok(n+1)=1.874409e+000; ng(n+1)=1.574114e+002;
n=270; farx(n+1)=7.476945e+000; foe(n+1)=2.780447e+001; krok(n+1)=2.658374e-001; ng(n+1)=3.739803e+002;
n=271; farx(n+1)=7.570673e+000; foe(n+1)=2.779712e+001; krok(n+1)=9.846165e-001; ng(n+1)=2.375334e+002;
n=272; farx(n+1)=7.581551e+000; foe(n+1)=2.779531e+001; krok(n+1)=1.914145e-001; ng(n+1)=3.405090e+002;
n=273; farx(n+1)=7.617598e+000; foe(n+1)=2.779052e+001; krok(n+1)=9.372047e-001; ng(n+1)=7.608327e+001;
n=274; farx(n+1)=7.649362e+000; foe(n+1)=2.778650e+001; krok(n+1)=7.310195e-001; ng(n+1)=3.523808e+002;
n=275; farx(n+1)=7.663145e+000; foe(n+1)=2.778481e+001; krok(n+1)=2.701916e-001; ng(n+1)=1.922774e+002;
%odnowa zmiennej metryki
n=276; farx(n+1)=7.663109e+000; foe(n+1)=2.778386e+001; krok(n+1)=8.192372e-008; ng(n+1)=2.952516e+002;
n=277; farx(n+1)=7.661093e+000; foe(n+1)=2.778303e+001; krok(n+1)=5.552069e-007; ng(n+1)=1.270835e+002;
n=278; farx(n+1)=7.660838e+000; foe(n+1)=2.778285e+001; krok(n+1)=3.213736e-006; ng(n+1)=2.091893e+001;
n=279; farx(n+1)=7.653318e+000; foe(n+1)=2.778266e+001; krok(n+1)=7.471755e-006; ng(n+1)=1.560568e+001;
n=280; farx(n+1)=7.647682e+000; foe(n+1)=2.778222e+001; krok(n+1)=3.999817e-005; ng(n+1)=1.129609e+001;
n=281; farx(n+1)=7.649268e+000; foe(n+1)=2.778214e+001; krok(n+1)=7.657542e-005; ng(n+1)=3.782752e+000;
n=282; farx(n+1)=7.652065e+000; foe(n+1)=2.778196e+001; krok(n+1)=3.647778e-004; ng(n+1)=3.100005e+000;
n=283; farx(n+1)=7.660719e+000; foe(n+1)=2.778100e+001; krok(n+1)=2.583925e-003; ng(n+1)=2.223415e+000;
n=284; farx(n+1)=7.665406e+000; foe(n+1)=2.777899e+001; krok(n+1)=6.042575e-003; ng(n+1)=3.219985e+000;
n=285; farx(n+1)=7.657141e+000; foe(n+1)=2.777783e+001; krok(n+1)=2.598953e-003; ng(n+1)=1.404843e+001;
n=286; farx(n+1)=7.632947e+000; foe(n+1)=2.777717e+001; krok(n+1)=5.478088e-003; ng(n+1)=2.124995e+001;
n=287; farx(n+1)=7.617874e+000; foe(n+1)=2.777622e+001; krok(n+1)=3.351559e-002; ng(n+1)=3.541772e+001;
n=288; farx(n+1)=7.604535e+000; foe(n+1)=2.777506e+001; krok(n+1)=8.403478e-003; ng(n+1)=2.078780e+001;
n=289; farx(n+1)=7.588899e+000; foe(n+1)=2.777401e+001; krok(n+1)=3.638604e-002; ng(n+1)=6.154898e+001;
n=290; farx(n+1)=7.494549e+000; foe(n+1)=2.776068e+001; krok(n+1)=7.847380e-001; ng(n+1)=1.100206e+002;
n=291; farx(n+1)=7.475582e+000; foe(n+1)=2.775745e+001; krok(n+1)=9.595947e-002; ng(n+1)=6.135178e+001;
n=292; farx(n+1)=7.408859e+000; foe(n+1)=2.775226e+001; krok(n+1)=6.345427e-002; ng(n+1)=2.066566e+002;
n=293; farx(n+1)=7.370278e+000; foe(n+1)=2.774740e+001; krok(n+1)=3.655097e-001; ng(n+1)=1.879820e+002;
n=294; farx(n+1)=7.334096e+000; foe(n+1)=2.773862e+001; krok(n+1)=7.011953e-001; ng(n+1)=2.592337e+002;
n=295; farx(n+1)=7.273088e+000; foe(n+1)=2.773376e+001; krok(n+1)=1.435731e+000; ng(n+1)=6.644644e+001;
n=296; farx(n+1)=7.261436e+000; foe(n+1)=2.773249e+001; krok(n+1)=5.556839e-001; ng(n+1)=1.485388e+002;
n=297; farx(n+1)=7.270096e+000; foe(n+1)=2.773152e+001; krok(n+1)=9.744898e-001; ng(n+1)=8.771382e+001;
n=298; farx(n+1)=7.265889e+000; foe(n+1)=2.773065e+001; krok(n+1)=1.800219e+000; ng(n+1)=1.238649e+002;
n=299; farx(n+1)=7.272737e+000; foe(n+1)=2.772997e+001; krok(n+1)=1.154088e+000; ng(n+1)=7.223061e+001;
n=300; farx(n+1)=7.328313e+000; foe(n+1)=2.772811e+001; krok(n+1)=2.446161e+000; ng(n+1)=6.200856e+001;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
