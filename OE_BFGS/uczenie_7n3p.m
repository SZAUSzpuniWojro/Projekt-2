%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.681605e+003; foe(n+1)=4.556457e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=4.224645e+003; foe(n+1)=3.871533e+003; krok(n+1)=5.897838e-004; ng(n+1)=2.919047e+003;
n=2; farx(n+1)=1.376669e+003; foe(n+1)=9.022481e+002; krok(n+1)=3.437609e-003; ng(n+1)=2.039184e+003;
n=3; farx(n+1)=1.548654e+003; foe(n+1)=7.743544e+002; krok(n+1)=2.661414e-004; ng(n+1)=2.705541e+003;
n=4; farx(n+1)=7.668711e+002; foe(n+1)=6.647373e+002; krok(n+1)=1.915637e-003; ng(n+1)=9.893008e+002;
n=5; farx(n+1)=4.497109e+002; foe(n+1)=4.836268e+002; krok(n+1)=5.623096e-004; ng(n+1)=3.875376e+003;
n=6; farx(n+1)=4.239217e+002; foe(n+1)=4.499928e+002; krok(n+1)=2.110724e-003; ng(n+1)=1.619991e+003;
n=7; farx(n+1)=3.975973e+002; foe(n+1)=4.346738e+002; krok(n+1)=1.532691e-003; ng(n+1)=1.137629e+003;
n=8; farx(n+1)=3.431899e+002; foe(n+1)=4.228082e+002; krok(n+1)=1.460933e-003; ng(n+1)=1.630200e+003;
n=9; farx(n+1)=2.669027e+002; foe(n+1)=3.958826e+002; krok(n+1)=1.291963e-003; ng(n+1)=1.377315e+003;
n=10; farx(n+1)=2.483718e+002; foe(n+1)=3.784917e+002; krok(n+1)=2.275225e-003; ng(n+1)=8.832961e+002;
n=11; farx(n+1)=2.373123e+002; foe(n+1)=3.267347e+002; krok(n+1)=2.121602e-003; ng(n+1)=1.979498e+003;
n=12; farx(n+1)=2.028822e+002; foe(n+1)=2.872766e+002; krok(n+1)=2.583925e-003; ng(n+1)=1.005012e+003;
n=13; farx(n+1)=1.225546e+002; foe(n+1)=2.551449e+002; krok(n+1)=1.242424e-003; ng(n+1)=2.437116e+003;
n=14; farx(n+1)=1.171009e+002; foe(n+1)=2.511578e+002; krok(n+1)=9.921626e-005; ng(n+1)=1.781584e+003;
n=15; farx(n+1)=9.738031e+001; foe(n+1)=2.433533e+002; krok(n+1)=7.740461e-004; ng(n+1)=1.606718e+003;
n=16; farx(n+1)=1.051169e+002; foe(n+1)=2.192813e+002; krok(n+1)=1.400894e-002; ng(n+1)=9.498186e+002;
n=17; farx(n+1)=1.089018e+002; foe(n+1)=2.147912e+002; krok(n+1)=9.064776e-004; ng(n+1)=1.104238e+003;
n=18; farx(n+1)=9.528812e+001; foe(n+1)=2.021856e+002; krok(n+1)=2.855545e-003; ng(n+1)=1.895318e+003;
n=19; farx(n+1)=8.528481e+001; foe(n+1)=1.926278e+002; krok(n+1)=2.861859e-004; ng(n+1)=3.257653e+003;
n=20; farx(n+1)=8.666808e+001; foe(n+1)=1.853375e+002; krok(n+1)=2.005784e-003; ng(n+1)=1.949586e+003;
n=21; farx(n+1)=8.043962e+001; foe(n+1)=1.674988e+002; krok(n+1)=6.256921e-003; ng(n+1)=2.179742e+003;
n=22; farx(n+1)=7.614312e+001; foe(n+1)=1.543492e+002; krok(n+1)=1.078226e-003; ng(n+1)=1.173487e+003;
n=23; farx(n+1)=7.204385e+001; foe(n+1)=1.465899e+002; krok(n+1)=2.926191e-003; ng(n+1)=1.210501e+003;
n=24; farx(n+1)=6.550555e+001; foe(n+1)=1.399496e+002; krok(n+1)=9.778245e-003; ng(n+1)=1.355035e+003;
n=25; farx(n+1)=6.373583e+001; foe(n+1)=1.283357e+002; krok(n+1)=1.175883e-002; ng(n+1)=1.674215e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=6.051784e+001; foe(n+1)=1.191675e+002; krok(n+1)=2.443486e-005; ng(n+1)=2.374830e+003;
n=27; farx(n+1)=5.619700e+001; foe(n+1)=1.088811e+002; krok(n+1)=4.673205e-005; ng(n+1)=1.776805e+003;
n=28; farx(n+1)=5.238096e+001; foe(n+1)=1.042859e+002; krok(n+1)=1.691914e-004; ng(n+1)=7.212410e+002;
n=29; farx(n+1)=4.722726e+001; foe(n+1)=9.727365e+001; krok(n+1)=5.081655e-004; ng(n+1)=7.403889e+002;
n=30; farx(n+1)=4.146225e+001; foe(n+1)=8.687596e+001; krok(n+1)=3.946779e-004; ng(n+1)=1.185442e+003;
n=31; farx(n+1)=3.734330e+001; foe(n+1)=7.949223e+001; krok(n+1)=2.414856e-003; ng(n+1)=1.132296e+003;
n=32; farx(n+1)=3.501905e+001; foe(n+1)=7.411640e+001; krok(n+1)=2.232051e-003; ng(n+1)=1.676214e+003;
n=33; farx(n+1)=3.060388e+001; foe(n+1)=6.841307e+001; krok(n+1)=3.187345e-003; ng(n+1)=1.121186e+003;
n=34; farx(n+1)=2.801350e+001; foe(n+1)=6.167158e+001; krok(n+1)=4.743099e-003; ng(n+1)=7.088559e+002;
n=35; farx(n+1)=2.496185e+001; foe(n+1)=5.793207e+001; krok(n+1)=1.739067e-003; ng(n+1)=5.609152e+002;
n=36; farx(n+1)=2.403130e+001; foe(n+1)=5.542136e+001; krok(n+1)=5.557598e-003; ng(n+1)=1.731897e+003;
n=37; farx(n+1)=2.192919e+001; foe(n+1)=5.188902e+001; krok(n+1)=8.076024e-003; ng(n+1)=5.243390e+002;
n=38; farx(n+1)=1.938833e+001; foe(n+1)=4.614424e+001; krok(n+1)=1.428969e-002; ng(n+1)=6.632385e+002;
n=39; farx(n+1)=1.796436e+001; foe(n+1)=4.411294e+001; krok(n+1)=3.187345e-003; ng(n+1)=8.243626e+002;
n=40; farx(n+1)=1.730300e+001; foe(n+1)=4.289942e+001; krok(n+1)=3.671029e-003; ng(n+1)=6.251805e+002;
n=41; farx(n+1)=1.591881e+001; foe(n+1)=4.033955e+001; krok(n+1)=1.725161e-002; ng(n+1)=6.265131e+002;
n=42; farx(n+1)=1.612404e+001; foe(n+1)=3.915660e+001; krok(n+1)=1.518514e-002; ng(n+1)=1.224269e+003;
n=43; farx(n+1)=1.601223e+001; foe(n+1)=3.590896e+001; krok(n+1)=1.214548e-002; ng(n+1)=1.494720e+003;
n=44; farx(n+1)=1.534679e+001; foe(n+1)=3.450709e+001; krok(n+1)=4.625835e-003; ng(n+1)=3.828106e+002;
n=45; farx(n+1)=1.537399e+001; foe(n+1)=3.387794e+001; krok(n+1)=4.312903e-003; ng(n+1)=9.525169e+002;
n=46; farx(n+1)=1.484372e+001; foe(n+1)=3.307971e+001; krok(n+1)=1.675011e-002; ng(n+1)=4.448101e+002;
n=47; farx(n+1)=1.419563e+001; foe(n+1)=3.167550e+001; krok(n+1)=7.557177e-003; ng(n+1)=3.650970e+002;
n=48; farx(n+1)=1.390044e+001; foe(n+1)=3.123230e+001; krok(n+1)=1.007231e-002; ng(n+1)=4.409769e+002;
n=49; farx(n+1)=1.342216e+001; foe(n+1)=3.054496e+001; krok(n+1)=1.356166e-002; ng(n+1)=9.289333e+002;
n=50; farx(n+1)=1.218418e+001; foe(n+1)=2.919116e+001; krok(n+1)=3.469665e-002; ng(n+1)=8.565297e+002;
%odnowa zmiennej metryki
n=51; farx(n+1)=1.204454e+001; foe(n+1)=2.842169e+001; krok(n+1)=7.258975e-006; ng(n+1)=1.204430e+003;
n=52; farx(n+1)=1.194314e+001; foe(n+1)=2.828209e+001; krok(n+1)=3.594329e-005; ng(n+1)=2.466366e+002;
n=53; farx(n+1)=1.203567e+001; foe(n+1)=2.815027e+001; krok(n+1)=1.278030e-004; ng(n+1)=1.248563e+002;
n=54; farx(n+1)=1.167950e+001; foe(n+1)=2.754727e+001; krok(n+1)=2.647552e-003; ng(n+1)=7.074354e+001;
n=55; farx(n+1)=1.088409e+001; foe(n+1)=2.661613e+001; krok(n+1)=2.257553e-003; ng(n+1)=1.027057e+002;
n=56; farx(n+1)=1.039357e+001; foe(n+1)=2.542386e+001; krok(n+1)=1.324250e-003; ng(n+1)=1.647511e+002;
n=57; farx(n+1)=9.597145e+000; foe(n+1)=2.456584e+001; krok(n+1)=2.041615e-003; ng(n+1)=2.603638e+002;
n=58; farx(n+1)=9.159968e+000; foe(n+1)=2.403488e+001; krok(n+1)=5.608389e-003; ng(n+1)=5.474775e+002;
n=59; farx(n+1)=9.288959e+000; foe(n+1)=2.363242e+001; krok(n+1)=7.152527e-003; ng(n+1)=8.895815e+002;
n=60; farx(n+1)=9.033843e+000; foe(n+1)=2.323226e+001; krok(n+1)=6.193383e-003; ng(n+1)=7.385456e+002;
n=61; farx(n+1)=9.038483e+000; foe(n+1)=2.245970e+001; krok(n+1)=1.462925e-002; ng(n+1)=8.640903e+002;
n=62; farx(n+1)=9.114329e+000; foe(n+1)=2.230131e+001; krok(n+1)=3.271540e-003; ng(n+1)=4.571111e+002;
n=63; farx(n+1)=8.728952e+000; foe(n+1)=2.195849e+001; krok(n+1)=1.097075e-002; ng(n+1)=3.524685e+002;
n=64; farx(n+1)=8.470806e+000; foe(n+1)=2.179293e+001; krok(n+1)=1.321526e-002; ng(n+1)=2.037757e+002;
n=65; farx(n+1)=8.430847e+000; foe(n+1)=2.152647e+001; krok(n+1)=3.627640e-002; ng(n+1)=1.765779e+002;
n=66; farx(n+1)=8.276317e+000; foe(n+1)=2.096145e+001; krok(n+1)=4.254539e-002; ng(n+1)=1.448008e+002;
n=67; farx(n+1)=8.172250e+000; foe(n+1)=2.063810e+001; krok(n+1)=2.854877e-002; ng(n+1)=2.981123e+002;
n=68; farx(n+1)=8.172508e+000; foe(n+1)=2.016892e+001; krok(n+1)=3.006981e-002; ng(n+1)=6.569121e+002;
n=69; farx(n+1)=7.581351e+000; foe(n+1)=1.896436e+001; krok(n+1)=5.695205e-002; ng(n+1)=2.469510e+002;
n=70; farx(n+1)=7.419280e+000; foe(n+1)=1.855630e+001; krok(n+1)=2.515521e-002; ng(n+1)=1.236930e+002;
n=71; farx(n+1)=7.255434e+000; foe(n+1)=1.840560e+001; krok(n+1)=2.191235e-002; ng(n+1)=2.864487e+002;
n=72; farx(n+1)=6.987244e+000; foe(n+1)=1.813073e+001; krok(n+1)=3.149998e-002; ng(n+1)=2.794698e+002;
n=73; farx(n+1)=6.791873e+000; foe(n+1)=1.739400e+001; krok(n+1)=6.306675e-002; ng(n+1)=2.171591e+002;
n=74; farx(n+1)=6.705215e+000; foe(n+1)=1.710170e+001; krok(n+1)=4.346837e-002; ng(n+1)=1.923149e+002;
n=75; farx(n+1)=6.704059e+000; foe(n+1)=1.659042e+001; krok(n+1)=1.532184e-001; ng(n+1)=1.719465e+002;
%odnowa zmiennej metryki
n=76; farx(n+1)=6.616697e+000; foe(n+1)=1.634082e+001; krok(n+1)=3.740683e-005; ng(n+1)=3.476111e+002;
n=77; farx(n+1)=6.618732e+000; foe(n+1)=1.631267e+001; krok(n+1)=1.308609e-005; ng(n+1)=1.671614e+002;
n=78; farx(n+1)=6.652119e+000; foe(n+1)=1.626152e+001; krok(n+1)=4.702476e-005; ng(n+1)=1.451984e+002;
n=79; farx(n+1)=6.422846e+000; foe(n+1)=1.603050e+001; krok(n+1)=1.358387e-003; ng(n+1)=5.521656e+001;
n=80; farx(n+1)=6.356308e+000; foe(n+1)=1.585497e+001; krok(n+1)=9.352282e-004; ng(n+1)=6.260487e+001;
n=81; farx(n+1)=6.243533e+000; foe(n+1)=1.572043e+001; krok(n+1)=2.133311e-003; ng(n+1)=3.846649e+001;
n=82; farx(n+1)=5.924816e+000; foe(n+1)=1.531346e+001; krok(n+1)=6.800237e-003; ng(n+1)=7.220121e+001;
n=83; farx(n+1)=5.864534e+000; foe(n+1)=1.511985e+001; krok(n+1)=4.236653e-003; ng(n+1)=1.836749e+002;
n=84; farx(n+1)=5.770540e+000; foe(n+1)=1.497967e+001; krok(n+1)=3.334617e-003; ng(n+1)=2.928279e+002;
n=85; farx(n+1)=5.538580e+000; foe(n+1)=1.475973e+001; krok(n+1)=1.238677e-002; ng(n+1)=4.307396e+002;
n=86; farx(n+1)=4.816041e+000; foe(n+1)=1.366783e+001; krok(n+1)=1.012568e-002; ng(n+1)=5.406781e+002;
n=87; farx(n+1)=4.699123e+000; foe(n+1)=1.356853e+001; krok(n+1)=2.041615e-003; ng(n+1)=4.508604e+002;
n=88; farx(n+1)=4.649760e+000; foe(n+1)=1.339821e+001; krok(n+1)=2.024365e-002; ng(n+1)=2.666084e+002;
n=89; farx(n+1)=4.518132e+000; foe(n+1)=1.290091e+001; krok(n+1)=3.262119e-002; ng(n+1)=7.206158e+002;
n=90; farx(n+1)=4.417848e+000; foe(n+1)=1.269204e+001; krok(n+1)=6.274942e-003; ng(n+1)=5.995086e+002;
n=91; farx(n+1)=4.253161e+000; foe(n+1)=1.248220e+001; krok(n+1)=8.996953e-003; ng(n+1)=4.452106e+002;
n=92; farx(n+1)=4.111113e+000; foe(n+1)=1.219582e+001; krok(n+1)=1.506159e-002; ng(n+1)=3.856931e+002;
n=93; farx(n+1)=4.058693e+000; foe(n+1)=1.210258e+001; krok(n+1)=9.045235e-003; ng(n+1)=5.021812e+002;
n=94; farx(n+1)=4.003344e+000; foe(n+1)=1.197386e+001; krok(n+1)=1.153579e-002; ng(n+1)=4.270402e+002;
n=95; farx(n+1)=3.824807e+000; foe(n+1)=1.172585e+001; krok(n+1)=5.730945e-002; ng(n+1)=7.680043e+002;
n=96; farx(n+1)=3.667720e+000; foe(n+1)=1.143362e+001; krok(n+1)=3.971457e-002; ng(n+1)=4.600803e+002;
n=97; farx(n+1)=3.024047e+000; foe(n+1)=1.061586e+001; krok(n+1)=7.636705e-002; ng(n+1)=5.433022e+002;
n=98; farx(n+1)=2.549764e+000; foe(n+1)=1.006918e+001; krok(n+1)=2.351648e-002; ng(n+1)=1.087200e+003;
n=99; farx(n+1)=2.450483e+000; foe(n+1)=8.566172e+000; krok(n+1)=1.269085e-001; ng(n+1)=6.002151e+002;
n=100; farx(n+1)=2.575428e+000; foe(n+1)=8.149956e+000; krok(n+1)=4.254539e-002; ng(n+1)=4.971149e+002;
%odnowa zmiennej metryki
n=101; farx(n+1)=2.572117e+000; foe(n+1)=7.997473e+000; krok(n+1)=7.578315e-006; ng(n+1)=5.200555e+002;
n=102; farx(n+1)=2.540411e+000; foe(n+1)=7.862596e+000; krok(n+1)=6.070286e-006; ng(n+1)=5.569939e+002;
n=103; farx(n+1)=2.572318e+000; foe(n+1)=7.702808e+000; krok(n+1)=3.616555e-005; ng(n+1)=2.671026e+002;
n=104; farx(n+1)=2.552444e+000; foe(n+1)=7.420360e+000; krok(n+1)=2.224690e-004; ng(n+1)=1.589392e+002;
n=105; farx(n+1)=2.525770e+000; foe(n+1)=7.191729e+000; krok(n+1)=9.854180e-004; ng(n+1)=6.654306e+001;
n=106; farx(n+1)=2.535405e+000; foe(n+1)=6.968089e+000; krok(n+1)=1.120856e-003; ng(n+1)=6.539979e+001;
n=107; farx(n+1)=2.604887e+000; foe(n+1)=6.457667e+000; krok(n+1)=3.715831e-003; ng(n+1)=1.076172e+002;
n=108; farx(n+1)=2.601939e+000; foe(n+1)=6.269250e+000; krok(n+1)=2.187783e-003; ng(n+1)=4.378922e+002;
n=109; farx(n+1)=2.612739e+000; foe(n+1)=6.172341e+000; krok(n+1)=2.677949e-003; ng(n+1)=7.532539e+002;
n=110; farx(n+1)=2.591426e+000; foe(n+1)=6.079392e+000; krok(n+1)=2.707062e-003; ng(n+1)=6.475135e+002;
n=111; farx(n+1)=2.555114e+000; foe(n+1)=5.845303e+000; krok(n+1)=9.424758e-003; ng(n+1)=5.900908e+002;
n=112; farx(n+1)=2.640300e+000; foe(n+1)=5.678640e+000; krok(n+1)=8.579668e-003; ng(n+1)=1.987669e+002;
n=113; farx(n+1)=2.652706e+000; foe(n+1)=5.519286e+000; krok(n+1)=1.366559e-002; ng(n+1)=5.472088e+002;
n=114; farx(n+1)=2.544945e+000; foe(n+1)=5.434625e+000; krok(n+1)=1.423591e-002; ng(n+1)=1.697628e+002;
n=115; farx(n+1)=2.475978e+000; foe(n+1)=5.331076e+000; krok(n+1)=3.153693e-002; ng(n+1)=4.388373e+002;
n=116; farx(n+1)=2.394576e+000; foe(n+1)=5.221637e+000; krok(n+1)=2.667694e-002; ng(n+1)=4.800403e+002;
n=117; farx(n+1)=2.346866e+000; foe(n+1)=5.096714e+000; krok(n+1)=2.191235e-002; ng(n+1)=4.078943e+002;
n=118; farx(n+1)=2.338853e+000; foe(n+1)=5.038031e+000; krok(n+1)=2.876670e-002; ng(n+1)=2.751719e+002;
n=119; farx(n+1)=2.344724e+000; foe(n+1)=4.911422e+000; krok(n+1)=4.568872e-002; ng(n+1)=5.712027e+002;
n=120; farx(n+1)=2.164623e+000; foe(n+1)=4.675252e+000; krok(n+1)=1.036878e-001; ng(n+1)=3.212690e+002;
n=121; farx(n+1)=2.043474e+000; foe(n+1)=4.591525e+000; krok(n+1)=2.667694e-002; ng(n+1)=6.029017e+002;
n=122; farx(n+1)=1.829036e+000; foe(n+1)=4.491260e+000; krok(n+1)=9.363811e-002; ng(n+1)=5.604926e+002;
n=123; farx(n+1)=1.681047e+000; foe(n+1)=4.142933e+000; krok(n+1)=2.401907e-001; ng(n+1)=3.408558e+002;
n=124; farx(n+1)=1.844183e+000; foe(n+1)=4.013149e+000; krok(n+1)=6.240867e-002; ng(n+1)=4.629894e+002;
n=125; farx(n+1)=1.665493e+000; foe(n+1)=3.656986e+000; krok(n+1)=2.268451e-001; ng(n+1)=2.399326e+002;
%odnowa zmiennej metryki
n=126; farx(n+1)=1.652751e+000; foe(n+1)=3.536708e+000; krok(n+1)=1.364051e-006; ng(n+1)=9.617482e+002;
n=127; farx(n+1)=1.651933e+000; foe(n+1)=3.532461e+000; krok(n+1)=8.431791e-006; ng(n+1)=8.864810e+001;
n=128; farx(n+1)=1.641033e+000; foe(n+1)=3.504059e+000; krok(n+1)=1.239407e-004; ng(n+1)=5.640823e+001;
n=129; farx(n+1)=1.641284e+000; foe(n+1)=3.466398e+000; krok(n+1)=5.561724e-005; ng(n+1)=1.079555e+002;
n=130; farx(n+1)=1.636093e+000; foe(n+1)=3.447283e+000; krok(n+1)=2.921277e-004; ng(n+1)=3.293965e+001;
n=131; farx(n+1)=1.637421e+000; foe(n+1)=3.433808e+000; krok(n+1)=4.356528e-004; ng(n+1)=2.911947e+001;
n=132; farx(n+1)=1.610332e+000; foe(n+1)=3.401877e+000; krok(n+1)=1.667309e-003; ng(n+1)=2.151328e+001;
n=133; farx(n+1)=1.559988e+000; foe(n+1)=3.318244e+000; krok(n+1)=1.841613e-003; ng(n+1)=3.373981e+001;
n=134; farx(n+1)=1.505507e+000; foe(n+1)=3.212987e+000; krok(n+1)=3.485331e-003; ng(n+1)=5.385673e+001;
n=135; farx(n+1)=1.456323e+000; foe(n+1)=3.167883e+000; krok(n+1)=7.444102e-003; ng(n+1)=2.495124e+002;
n=136; farx(n+1)=1.412727e+000; foe(n+1)=3.128829e+000; krok(n+1)=1.376948e-002; ng(n+1)=3.566616e+002;
n=137; farx(n+1)=1.337664e+000; foe(n+1)=3.080095e+000; krok(n+1)=1.290032e-002; ng(n+1)=5.816955e+002;
n=138; farx(n+1)=1.265858e+000; foe(n+1)=3.056214e+000; krok(n+1)=1.290032e-002; ng(n+1)=6.142899e+002;
n=139; farx(n+1)=1.222497e+000; foe(n+1)=3.024526e+000; krok(n+1)=2.643052e-002; ng(n+1)=5.844783e+002;
n=140; farx(n+1)=1.178966e+000; foe(n+1)=3.006331e+000; krok(n+1)=1.023953e-002; ng(n+1)=5.620461e+002;
n=141; farx(n+1)=1.111837e+000; foe(n+1)=2.945684e+000; krok(n+1)=3.353295e-002; ng(n+1)=5.424162e+002;
n=142; farx(n+1)=1.120205e+000; foe(n+1)=2.914325e+000; krok(n+1)=3.481722e-002; ng(n+1)=2.998302e+002;
n=143; farx(n+1)=1.088273e+000; foe(n+1)=2.860917e+000; krok(n+1)=2.447549e-002; ng(n+1)=4.342955e+002;
n=144; farx(n+1)=1.053306e+000; foe(n+1)=2.801027e+000; krok(n+1)=5.514490e-002; ng(n+1)=2.219117e+002;
n=145; farx(n+1)=9.963697e-001; foe(n+1)=2.721816e+000; krok(n+1)=6.414287e-002; ng(n+1)=3.528561e+002;
n=146; farx(n+1)=1.032991e+000; foe(n+1)=2.651503e+000; krok(n+1)=3.178877e-002; ng(n+1)=2.171383e+002;
n=147; farx(n+1)=1.032097e+000; foe(n+1)=2.594617e+000; krok(n+1)=2.812841e-002; ng(n+1)=2.344393e+002;
n=148; farx(n+1)=8.853465e-001; foe(n+1)=2.470595e+000; krok(n+1)=1.710829e-001; ng(n+1)=3.821736e+002;
n=149; farx(n+1)=8.265815e-001; foe(n+1)=2.344888e+000; krok(n+1)=1.191056e-001; ng(n+1)=3.593652e+002;
n=150; farx(n+1)=8.848061e-001; foe(n+1)=2.245456e+000; krok(n+1)=6.287977e-002; ng(n+1)=3.648062e+002;
%odnowa zmiennej metryki
n=151; farx(n+1)=8.841560e-001; foe(n+1)=2.213395e+000; krok(n+1)=2.998147e-006; ng(n+1)=4.091508e+002;
n=152; farx(n+1)=8.843812e-001; foe(n+1)=2.211880e+000; krok(n+1)=2.296189e-006; ng(n+1)=1.008834e+002;
n=153; farx(n+1)=8.859656e-001; foe(n+1)=2.205491e+000; krok(n+1)=8.106422e-006; ng(n+1)=8.605778e+001;
n=154; farx(n+1)=8.876725e-001; foe(n+1)=2.190804e+000; krok(n+1)=1.054119e-004; ng(n+1)=4.251971e+001;
n=155; farx(n+1)=8.855170e-001; foe(n+1)=2.184362e+000; krok(n+1)=4.532388e-004; ng(n+1)=1.825907e+001;
n=156; farx(n+1)=8.803027e-001; foe(n+1)=2.152773e+000; krok(n+1)=6.767654e-004; ng(n+1)=2.791172e+001;
n=157; farx(n+1)=8.742889e-001; foe(n+1)=2.098384e+000; krok(n+1)=7.964319e-004; ng(n+1)=3.908536e+001;
n=158; farx(n+1)=8.702824e-001; foe(n+1)=2.079878e+000; krok(n+1)=6.459813e-004; ng(n+1)=1.307749e+002;
n=159; farx(n+1)=8.681838e-001; foe(n+1)=2.055074e+000; krok(n+1)=2.106031e-003; ng(n+1)=1.905993e+002;
n=160; farx(n+1)=8.686041e-001; foe(n+1)=2.045506e+000; krok(n+1)=2.595225e-003; ng(n+1)=2.928901e+002;
n=161; farx(n+1)=9.013112e-001; foe(n+1)=1.995146e+000; krok(n+1)=2.047906e-002; ng(n+1)=3.282573e+002;
n=162; farx(n+1)=9.123027e-001; foe(n+1)=1.959985e+000; krok(n+1)=3.427028e-003; ng(n+1)=4.422603e+002;
n=163; farx(n+1)=8.672797e-001; foe(n+1)=1.913477e+000; krok(n+1)=1.169933e-002; ng(n+1)=2.723498e+002;
n=164; farx(n+1)=8.230301e-001; foe(n+1)=1.895242e+000; krok(n+1)=1.520419e-002; ng(n+1)=1.988526e+002;
n=165; farx(n+1)=8.189482e-001; foe(n+1)=1.887533e+000; krok(n+1)=1.526430e-002; ng(n+1)=1.852033e+002;
n=166; farx(n+1)=8.221489e-001; foe(n+1)=1.869299e+000; krok(n+1)=3.389322e-002; ng(n+1)=1.736486e+002;
n=167; farx(n+1)=7.803867e-001; foe(n+1)=1.797799e+000; krok(n+1)=5.730115e-002; ng(n+1)=2.311838e+002;
n=168; farx(n+1)=7.628330e-001; foe(n+1)=1.773823e+000; krok(n+1)=9.100900e-003; ng(n+1)=3.838850e+002;
n=169; farx(n+1)=7.473045e-001; foe(n+1)=1.761302e+000; krok(n+1)=2.079162e-002; ng(n+1)=3.026240e+002;
n=170; farx(n+1)=7.290037e-001; foe(n+1)=1.735904e+000; krok(n+1)=3.407387e-002; ng(n+1)=4.917439e+002;
n=171; farx(n+1)=6.976708e-001; foe(n+1)=1.608192e+000; krok(n+1)=1.400181e-001; ng(n+1)=3.868472e+002;
n=172; farx(n+1)=6.718596e-001; foe(n+1)=1.478176e+000; krok(n+1)=5.507577e-002; ng(n+1)=7.580782e+002;
n=173; farx(n+1)=6.605194e-001; foe(n+1)=1.431745e+000; krok(n+1)=7.822596e-002; ng(n+1)=2.999915e+002;
n=174; farx(n+1)=6.197965e-001; foe(n+1)=1.374168e+000; krok(n+1)=4.235232e-002; ng(n+1)=3.125665e+002;
n=175; farx(n+1)=6.210742e-001; foe(n+1)=1.342186e+000; krok(n+1)=2.262062e-002; ng(n+1)=5.467178e+002;
%odnowa zmiennej metryki
n=176; farx(n+1)=6.209638e-001; foe(n+1)=1.338064e+000; krok(n+1)=1.341315e-006; ng(n+1)=2.097451e+002;
n=177; farx(n+1)=6.207570e-001; foe(n+1)=1.336897e+000; krok(n+1)=1.192618e-006; ng(n+1)=1.139202e+002;
n=178; farx(n+1)=6.202370e-001; foe(n+1)=1.333933e+000; krok(n+1)=4.269784e-006; ng(n+1)=1.237900e+002;
n=179; farx(n+1)=6.214219e-001; foe(n+1)=1.328060e+000; krok(n+1)=1.431236e-004; ng(n+1)=2.576330e+001;
n=180; farx(n+1)=6.215794e-001; foe(n+1)=1.323043e+000; krok(n+1)=1.629504e-004; ng(n+1)=2.604590e+001;
n=181; farx(n+1)=6.218758e-001; foe(n+1)=1.313812e+000; krok(n+1)=2.942857e-004; ng(n+1)=2.934503e+001;
n=182; farx(n+1)=6.248181e-001; foe(n+1)=1.300576e+000; krok(n+1)=3.657739e-004; ng(n+1)=2.661197e+001;
n=183; farx(n+1)=6.277830e-001; foe(n+1)=1.288002e+000; krok(n+1)=4.532388e-004; ng(n+1)=4.812685e+001;
n=184; farx(n+1)=6.301154e-001; foe(n+1)=1.276763e+000; krok(n+1)=1.442683e-003; ng(n+1)=7.970101e+001;
n=185; farx(n+1)=6.316265e-001; foe(n+1)=1.259980e+000; krok(n+1)=4.494346e-003; ng(n+1)=1.232130e+002;
n=186; farx(n+1)=6.280300e-001; foe(n+1)=1.245480e+000; krok(n+1)=6.227174e-003; ng(n+1)=1.703381e+002;
n=187; farx(n+1)=6.167710e-001; foe(n+1)=1.231384e+000; krok(n+1)=1.660195e-002; ng(n+1)=3.285326e+002;
n=188; farx(n+1)=6.160587e-001; foe(n+1)=1.198931e+000; krok(n+1)=2.753789e-002; ng(n+1)=6.387630e+002;
n=189; farx(n+1)=6.100336e-001; foe(n+1)=1.173781e+000; krok(n+1)=9.683310e-003; ng(n+1)=6.802686e+002;
n=190; farx(n+1)=6.062131e-001; foe(n+1)=1.168322e+000; krok(n+1)=9.409526e-003; ng(n+1)=3.189018e+002;
n=191; farx(n+1)=5.998212e-001; foe(n+1)=1.157626e+000; krok(n+1)=1.056599e-002; ng(n+1)=3.328779e+002;
n=192; farx(n+1)=6.000299e-001; foe(n+1)=1.146617e+000; krok(n+1)=2.085765e-002; ng(n+1)=2.373938e+002;
n=193; farx(n+1)=5.937801e-001; foe(n+1)=1.137627e+000; krok(n+1)=2.351648e-002; ng(n+1)=2.417461e+002;
n=194; farx(n+1)=5.987641e-001; foe(n+1)=1.128723e+000; krok(n+1)=6.288276e-002; ng(n+1)=3.009637e+002;
n=195; farx(n+1)=5.732238e-001; foe(n+1)=1.105406e+000; krok(n+1)=1.081188e-001; ng(n+1)=1.935483e+002;
n=196; farx(n+1)=5.734607e-001; foe(n+1)=1.098746e+000; krok(n+1)=1.952078e-002; ng(n+1)=2.214451e+002;
n=197; farx(n+1)=5.785214e-001; foe(n+1)=1.078781e+000; krok(n+1)=1.076512e-001; ng(n+1)=8.740062e+001;
n=198; farx(n+1)=5.684309e-001; foe(n+1)=1.058440e+000; krok(n+1)=1.301252e-002; ng(n+1)=2.595985e+002;
n=199; farx(n+1)=5.593189e-001; foe(n+1)=1.042984e+000; krok(n+1)=2.788482e-002; ng(n+1)=2.462272e+002;
n=200; farx(n+1)=5.519258e-001; foe(n+1)=1.016841e+000; krok(n+1)=5.801457e-002; ng(n+1)=4.729601e+002;
%odnowa zmiennej metryki
n=201; farx(n+1)=5.519040e-001; foe(n+1)=1.011853e+000; krok(n+1)=6.104072e-007; ng(n+1)=3.463003e+002;
n=202; farx(n+1)=5.518679e-001; foe(n+1)=1.010029e+000; krok(n+1)=7.440534e-007; ng(n+1)=2.064318e+002;
n=203; farx(n+1)=5.514519e-001; foe(n+1)=1.004015e+000; krok(n+1)=1.083629e-005; ng(n+1)=8.886696e+001;
n=204; farx(n+1)=5.518791e-001; foe(n+1)=9.937380e-001; krok(n+1)=2.609816e-005; ng(n+1)=8.865456e+001;
n=205; farx(n+1)=5.524401e-001; foe(n+1)=9.894599e-001; krok(n+1)=1.314970e-004; ng(n+1)=2.724505e+001;
n=206; farx(n+1)=5.526491e-001; foe(n+1)=9.819382e-001; krok(n+1)=2.996980e-004; ng(n+1)=2.243180e+001;
n=207; farx(n+1)=5.539044e-001; foe(n+1)=9.747207e-001; krok(n+1)=2.808967e-004; ng(n+1)=2.720706e+001;
n=208; farx(n+1)=5.544378e-001; foe(n+1)=9.715114e-001; krok(n+1)=4.264342e-004; ng(n+1)=2.080538e+001;
n=209; farx(n+1)=5.554893e-001; foe(n+1)=9.680220e-001; krok(n+1)=1.064566e-003; ng(n+1)=2.111551e+001;
n=210; farx(n+1)=5.601939e-001; foe(n+1)=9.437064e-001; krok(n+1)=8.928204e-003; ng(n+1)=3.237918e+001;
n=211; farx(n+1)=5.562611e-001; foe(n+1)=9.366999e-001; krok(n+1)=6.947618e-003; ng(n+1)=2.862416e+002;
n=212; farx(n+1)=5.539196e-001; foe(n+1)=9.299742e-001; krok(n+1)=9.219686e-003; ng(n+1)=4.693570e+002;
n=213; farx(n+1)=5.549737e-001; foe(n+1)=9.203424e-001; krok(n+1)=7.152527e-003; ng(n+1)=5.198689e+002;
n=214; farx(n+1)=5.560503e-001; foe(n+1)=9.156584e-001; krok(n+1)=9.949946e-003; ng(n+1)=4.752850e+002;
n=215; farx(n+1)=5.607474e-001; foe(n+1)=9.119298e-001; krok(n+1)=2.721013e-002; ng(n+1)=2.330883e+002;
n=216; farx(n+1)=5.579190e-001; foe(n+1)=9.069447e-001; krok(n+1)=2.896961e-002; ng(n+1)=2.735482e+002;
n=217; farx(n+1)=5.568290e-001; foe(n+1)=9.018838e-001; krok(n+1)=2.780492e-002; ng(n+1)=1.941770e+002;
n=218; farx(n+1)=5.620448e-001; foe(n+1)=8.972920e-001; krok(n+1)=1.767909e-002; ng(n+1)=1.910616e+002;
n=219; farx(n+1)=5.598914e-001; foe(n+1)=8.905690e-001; krok(n+1)=5.042484e-002; ng(n+1)=1.087559e+002;
n=220; farx(n+1)=5.568572e-001; foe(n+1)=8.804311e-001; krok(n+1)=4.890625e-002; ng(n+1)=1.944873e+002;
n=221; farx(n+1)=5.486155e-001; foe(n+1)=8.642603e-001; krok(n+1)=1.171506e-001; ng(n+1)=1.803174e+002;
n=222; farx(n+1)=5.454600e-001; foe(n+1)=8.488066e-001; krok(n+1)=2.868906e-002; ng(n+1)=1.620505e+002;
n=223; farx(n+1)=5.471696e-001; foe(n+1)=8.384989e-001; krok(n+1)=3.408723e-002; ng(n+1)=2.377653e+002;
n=224; farx(n+1)=5.477138e-001; foe(n+1)=8.174266e-001; krok(n+1)=9.759878e-002; ng(n+1)=3.226130e+002;
n=225; farx(n+1)=5.489173e-001; foe(n+1)=8.028349e-001; krok(n+1)=9.137743e-002; ng(n+1)=2.613543e+002;
%odnowa zmiennej metryki
n=226; farx(n+1)=5.488909e-001; foe(n+1)=8.004365e-001; krok(n+1)=2.933371e-007; ng(n+1)=3.204255e+002;
n=227; farx(n+1)=5.487991e-001; foe(n+1)=7.997778e-001; krok(n+1)=4.136318e-006; ng(n+1)=5.652075e+001;
n=228; farx(n+1)=5.486475e-001; foe(n+1)=7.986903e-001; krok(n+1)=2.208309e-006; ng(n+1)=9.341530e+001;
n=229; farx(n+1)=5.485590e-001; foe(n+1)=7.947897e-001; krok(n+1)=2.832743e-005; ng(n+1)=5.609525e+001;
n=230; farx(n+1)=5.483264e-001; foe(n+1)=7.874487e-001; krok(n+1)=1.781286e-004; ng(n+1)=2.620810e+001;
n=231; farx(n+1)=5.488631e-001; foe(n+1)=7.816578e-001; krok(n+1)=1.486449e-004; ng(n+1)=2.514717e+001;
n=232; farx(n+1)=5.490429e-001; foe(n+1)=7.775141e-001; krok(n+1)=1.802468e-004; ng(n+1)=2.726552e+001;
n=233; farx(n+1)=5.492752e-001; foe(n+1)=7.741986e-001; krok(n+1)=3.984182e-004; ng(n+1)=3.246204e+001;
n=234; farx(n+1)=5.493798e-001; foe(n+1)=7.699983e-001; krok(n+1)=8.747958e-004; ng(n+1)=3.818281e+001;
n=235; farx(n+1)=5.466899e-001; foe(n+1)=7.515623e-001; krok(n+1)=5.843731e-003; ng(n+1)=4.881004e+001;
n=236; farx(n+1)=5.469030e-001; foe(n+1)=7.490218e-001; krok(n+1)=3.009584e-003; ng(n+1)=2.118127e+002;
n=237; farx(n+1)=5.474769e-001; foe(n+1)=7.457076e-001; krok(n+1)=8.839546e-003; ng(n+1)=2.990108e+002;
n=238; farx(n+1)=5.453058e-001; foe(n+1)=7.392327e-001; krok(n+1)=5.918252e-002; ng(n+1)=4.288805e+002;
n=239; farx(n+1)=5.482780e-001; foe(n+1)=7.357799e-001; krok(n+1)=1.419977e-002; ng(n+1)=2.748903e+002;
n=240; farx(n+1)=5.468748e-001; foe(n+1)=7.343232e-001; krok(n+1)=1.417782e-002; ng(n+1)=1.731871e+002;
n=241; farx(n+1)=5.461736e-001; foe(n+1)=7.294132e-001; krok(n+1)=2.667694e-002; ng(n+1)=2.516654e+002;
n=242; farx(n+1)=5.435215e-001; foe(n+1)=7.251172e-001; krok(n+1)=1.966098e-002; ng(n+1)=1.723928e+002;
n=243; farx(n+1)=5.410237e-001; foe(n+1)=7.084716e-001; krok(n+1)=3.350022e-002; ng(n+1)=1.280759e+002;
n=244; farx(n+1)=5.381029e-001; foe(n+1)=6.956173e-001; krok(n+1)=1.251384e-002; ng(n+1)=3.174638e+002;
n=245; farx(n+1)=5.348850e-001; foe(n+1)=6.888812e-001; krok(n+1)=3.377159e-002; ng(n+1)=2.000686e+002;
n=246; farx(n+1)=5.351792e-001; foe(n+1)=6.645480e-001; krok(n+1)=2.725909e-001; ng(n+1)=1.931604e+002;
n=247; farx(n+1)=5.365808e-001; foe(n+1)=6.533964e-001; krok(n+1)=3.526755e-002; ng(n+1)=5.033107e+002;
n=248; farx(n+1)=5.369986e-001; foe(n+1)=6.511060e-001; krok(n+1)=9.267996e-003; ng(n+1)=4.522558e+002;
n=249; farx(n+1)=5.341761e-001; foe(n+1)=6.387679e-001; krok(n+1)=1.196783e-001; ng(n+1)=3.808583e+002;
n=250; farx(n+1)=5.298370e-001; foe(n+1)=6.102390e-001; krok(n+1)=7.557928e-002; ng(n+1)=5.641372e+002;
%odnowa zmiennej metryki
n=251; farx(n+1)=5.298059e-001; foe(n+1)=6.086303e-001; krok(n+1)=2.533257e-007; ng(n+1)=3.124837e+002;
n=252; farx(n+1)=5.296504e-001; foe(n+1)=6.054537e-001; krok(n+1)=4.782160e-006; ng(n+1)=1.033544e+002;
n=253; farx(n+1)=5.296334e-001; foe(n+1)=6.049839e-001; krok(n+1)=1.052955e-006; ng(n+1)=8.745827e+001;
n=254; farx(n+1)=5.296840e-001; foe(n+1)=6.016381e-001; krok(n+1)=2.122479e-005; ng(n+1)=5.982022e+001;
n=255; farx(n+1)=5.301621e-001; foe(n+1)=5.985113e-001; krok(n+1)=3.748417e-004; ng(n+1)=1.340229e+001;
n=256; farx(n+1)=5.303697e-001; foe(n+1)=5.969260e-001; krok(n+1)=1.300680e-004; ng(n+1)=1.580451e+001;
n=257; farx(n+1)=5.309845e-001; foe(n+1)=5.938068e-001; krok(n+1)=2.911565e-004; ng(n+1)=1.492938e+001;
n=258; farx(n+1)=5.309792e-001; foe(n+1)=5.929511e-001; krok(n+1)=2.015675e-004; ng(n+1)=1.389004e+001;
n=259; farx(n+1)=5.308942e-001; foe(n+1)=5.877143e-001; krok(n+1)=3.677109e-003; ng(n+1)=1.099799e+001;
n=260; farx(n+1)=5.311356e-001; foe(n+1)=5.803544e-001; krok(n+1)=1.093891e-003; ng(n+1)=4.627577e+001;
n=261; farx(n+1)=5.317040e-001; foe(n+1)=5.784972e-001; krok(n+1)=3.881663e-003; ng(n+1)=1.842764e+002;
n=262; farx(n+1)=5.335091e-001; foe(n+1)=5.717135e-001; krok(n+1)=3.990004e-003; ng(n+1)=2.469416e+002;
n=263; farx(n+1)=5.315306e-001; foe(n+1)=5.639676e-001; krok(n+1)=1.080437e-002; ng(n+1)=4.817573e+002;
n=264; farx(n+1)=5.325458e-001; foe(n+1)=5.594732e-001; krok(n+1)=6.256921e-003; ng(n+1)=5.260208e+002;
n=265; farx(n+1)=5.347602e-001; foe(n+1)=5.544264e-001; krok(n+1)=1.521910e-002; ng(n+1)=5.582098e+002;
n=266; farx(n+1)=5.371438e-001; foe(n+1)=5.482695e-001; krok(n+1)=3.627640e-002; ng(n+1)=3.382290e+002;
n=267; farx(n+1)=5.421096e-001; foe(n+1)=5.431440e-001; krok(n+1)=1.789318e-002; ng(n+1)=1.662214e+002;
n=268; farx(n+1)=5.430105e-001; foe(n+1)=5.399911e-001; krok(n+1)=1.199953e-002; ng(n+1)=2.024976e+002;
n=269; farx(n+1)=5.439493e-001; foe(n+1)=5.347079e-001; krok(n+1)=3.123706e-002; ng(n+1)=2.537793e+002;
n=270; farx(n+1)=5.439142e-001; foe(n+1)=5.296212e-001; krok(n+1)=6.078596e-002; ng(n+1)=2.065598e+002;
n=271; farx(n+1)=5.417229e-001; foe(n+1)=5.250128e-001; krok(n+1)=4.614317e-002; ng(n+1)=1.726370e+002;
n=272; farx(n+1)=5.407299e-001; foe(n+1)=5.206231e-001; krok(n+1)=1.871779e-002; ng(n+1)=3.751275e+002;
n=273; farx(n+1)=5.376763e-001; foe(n+1)=5.136957e-001; krok(n+1)=9.676700e-002; ng(n+1)=1.630944e+002;
n=274; farx(n+1)=5.357468e-001; foe(n+1)=5.082009e-001; krok(n+1)=7.165516e-002; ng(n+1)=1.905601e+002;
n=275; farx(n+1)=5.360724e-001; foe(n+1)=4.993970e-001; krok(n+1)=1.191056e-001; ng(n+1)=1.041588e+002;
%odnowa zmiennej metryki
n=276; farx(n+1)=5.360608e-001; foe(n+1)=4.964536e-001; krok(n+1)=3.646499e-007; ng(n+1)=3.396838e+002;
n=277; farx(n+1)=5.360490e-001; foe(n+1)=4.956643e-001; krok(n+1)=3.565875e-007; ng(n+1)=1.848873e+002;
n=278; farx(n+1)=5.360073e-001; foe(n+1)=4.951975e-001; krok(n+1)=4.182313e-006; ng(n+1)=4.267447e+001;
n=279; farx(n+1)=5.360479e-001; foe(n+1)=4.944616e-001; krok(n+1)=6.946355e-005; ng(n+1)=1.224419e+001;
n=280; farx(n+1)=5.360822e-001; foe(n+1)=4.937948e-001; krok(n+1)=7.598256e-005; ng(n+1)=1.396671e+001;
n=281; farx(n+1)=5.361187e-001; foe(n+1)=4.929862e-001; krok(n+1)=2.431917e-004; ng(n+1)=9.029476e+000;
n=282; farx(n+1)=5.361460e-001; foe(n+1)=4.927045e-001; krok(n+1)=7.856157e-005; ng(n+1)=8.500827e+000;
n=283; farx(n+1)=5.363149e-001; foe(n+1)=4.923848e-001; krok(n+1)=5.045766e-004; ng(n+1)=6.174574e+000;
n=284; farx(n+1)=5.363287e-001; foe(n+1)=4.920870e-001; krok(n+1)=6.650130e-004; ng(n+1)=8.221250e+000;
n=285; farx(n+1)=5.362651e-001; foe(n+1)=4.901946e-001; krok(n+1)=3.473809e-003; ng(n+1)=8.456589e+000;
n=286; farx(n+1)=5.358947e-001; foe(n+1)=4.868740e-001; krok(n+1)=1.742666e-003; ng(n+1)=2.150309e+001;
n=287; farx(n+1)=5.389025e-001; foe(n+1)=4.818576e-001; krok(n+1)=7.932207e-003; ng(n+1)=6.382214e+001;
n=288; farx(n+1)=5.400066e-001; foe(n+1)=4.798197e-001; krok(n+1)=8.685477e-003; ng(n+1)=2.443084e+002;
n=289; farx(n+1)=5.402707e-001; foe(n+1)=4.776313e-001; krok(n+1)=1.835014e-002; ng(n+1)=3.752841e+002;
n=290; farx(n+1)=5.424413e-001; foe(n+1)=4.752071e-001; krok(n+1)=3.301014e-002; ng(n+1)=3.773834e+002;
n=291; farx(n+1)=5.442703e-001; foe(n+1)=4.740827e-001; krok(n+1)=1.260621e-002; ng(n+1)=1.850414e+002;
n=292; farx(n+1)=5.438786e-001; foe(n+1)=4.723685e-001; krok(n+1)=5.985461e-002; ng(n+1)=1.569649e+002;
n=293; farx(n+1)=5.445254e-001; foe(n+1)=4.694601e-001; krok(n+1)=3.535818e-002; ng(n+1)=2.420344e+002;
n=294; farx(n+1)=5.432260e-001; foe(n+1)=4.658352e-001; krok(n+1)=8.509079e-002; ng(n+1)=2.041568e+002;
n=295; farx(n+1)=5.408560e-001; foe(n+1)=4.617866e-001; krok(n+1)=1.364675e-001; ng(n+1)=2.594455e+002;
n=296; farx(n+1)=5.415512e-001; foe(n+1)=4.591023e-001; krok(n+1)=3.338729e-002; ng(n+1)=1.558132e+002;
n=297; farx(n+1)=5.432849e-001; foe(n+1)=4.564223e-001; krok(n+1)=3.096333e-002; ng(n+1)=1.792447e+002;
n=298; farx(n+1)=5.447545e-001; foe(n+1)=4.513928e-001; krok(n+1)=1.467726e-001; ng(n+1)=1.580628e+002;
n=299; farx(n+1)=5.553091e-001; foe(n+1)=4.427596e-001; krok(n+1)=2.010759e-001; ng(n+1)=3.894186e+002;
n=300; farx(n+1)=5.564782e-001; foe(n+1)=4.369781e-001; krok(n+1)=1.075010e-001; ng(n+1)=1.828518e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)