%uczenie predyktora arx
clear all;
n=0; farx(n+1)=4.391345e+003; foe(n+1)=4.351952e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.793937e+003; foe(n+1)=3.953925e+003; krok(n+1)=4.674935e-004; ng(n+1)=3.402939e+003;
n=2; farx(n+1)=2.462988e+003; foe(n+1)=8.290012e+003; krok(n+1)=3.749852e-004; ng(n+1)=3.679217e+003;
n=3; farx(n+1)=1.718568e+003; foe(n+1)=1.410759e+004; krok(n+1)=4.168272e-004; ng(n+1)=6.692469e+003;
n=4; farx(n+1)=1.367454e+003; foe(n+1)=1.213837e+004; krok(n+1)=3.111744e-003; ng(n+1)=3.434241e+003;
n=5; farx(n+1)=1.184029e+003; foe(n+1)=1.117227e+004; krok(n+1)=4.168272e-004; ng(n+1)=3.105033e+003;
n=6; farx(n+1)=6.579270e+002; foe(n+1)=8.210664e+003; krok(n+1)=1.640492e-002; ng(n+1)=2.617495e+003;
n=7; farx(n+1)=4.750983e+002; foe(n+1)=7.224495e+003; krok(n+1)=2.133516e-003; ng(n+1)=2.009356e+003;
n=8; farx(n+1)=3.861024e+002; foe(n+1)=7.296568e+003; krok(n+1)=3.032635e-003; ng(n+1)=1.417818e+003;
n=9; farx(n+1)=1.521165e+002; foe(n+1)=4.151771e+003; krok(n+1)=6.518600e-003; ng(n+1)=1.355275e+003;
n=10; farx(n+1)=1.302311e+002; foe(n+1)=3.886759e+003; krok(n+1)=9.696269e-004; ng(n+1)=1.317698e+003;
n=11; farx(n+1)=1.268383e+002; foe(n+1)=3.860505e+003; krok(n+1)=1.694223e-003; ng(n+1)=1.306333e+003;
n=12; farx(n+1)=9.259379e+001; foe(n+1)=3.320386e+003; krok(n+1)=2.337492e-002; ng(n+1)=1.422676e+003;
n=13; farx(n+1)=8.417642e+001; foe(n+1)=3.086416e+003; krok(n+1)=2.147336e-002; ng(n+1)=8.541455e+002;
n=14; farx(n+1)=7.336612e+001; foe(n+1)=1.427521e+003; krok(n+1)=4.438543e-002; ng(n+1)=1.019256e+003;
n=15; farx(n+1)=4.987192e+001; foe(n+1)=8.454363e+003; krok(n+1)=1.362955e-001; ng(n+1)=8.292483e+002;
n=16; farx(n+1)=3.705778e+001; foe(n+1)=1.971663e+003; krok(n+1)=2.567403e-001; ng(n+1)=1.251048e+003;
n=17; farx(n+1)=2.258209e+001; foe(n+1)=3.056735e+003; krok(n+1)=1.150668e-001; ng(n+1)=7.775017e+002;
n=18; farx(n+1)=1.728644e+001; foe(n+1)=7.390646e+002; krok(n+1)=3.310405e-001; ng(n+1)=3.467958e+002;
n=19; farx(n+1)=1.203678e+001; foe(n+1)=5.338120e+002; krok(n+1)=1.037862e+000; ng(n+1)=3.491962e+002;
n=20; farx(n+1)=1.014478e+001; foe(n+1)=8.519898e+002; krok(n+1)=3.355407e-001; ng(n+1)=1.424957e+002;
n=21; farx(n+1)=7.394642e+000; foe(n+1)=1.519304e+003; krok(n+1)=1.000418e+000; ng(n+1)=1.408466e+002;
n=22; farx(n+1)=6.818115e+000; foe(n+1)=9.730432e+002; krok(n+1)=2.715719e-001; ng(n+1)=1.993895e+002;
n=23; farx(n+1)=6.096405e+000; foe(n+1)=5.436960e+003; krok(n+1)=5.461275e-001; ng(n+1)=1.308284e+002;
n=24; farx(n+1)=5.502814e+000; foe(n+1)=2.774996e+002; krok(n+1)=1.111368e+000; ng(n+1)=9.567768e+001;
n=25; farx(n+1)=4.761332e+000; foe(n+1)=2.545192e+002; krok(n+1)=2.087756e+000; ng(n+1)=8.991007e+001;
%odnowa zmiennej metryki
n=26; farx(n+1)=4.734150e+000; foe(n+1)=2.326882e+002; krok(n+1)=1.966794e-004; ng(n+1)=3.714180e+001;
n=27; farx(n+1)=4.681884e+000; foe(n+1)=2.199448e+002; krok(n+1)=9.854180e-004; ng(n+1)=2.566833e+001;
n=28; farx(n+1)=4.507165e+000; foe(n+1)=2.283678e+002; krok(n+1)=6.502194e-004; ng(n+1)=4.562861e+001;
n=29; farx(n+1)=4.090704e+000; foe(n+1)=2.455211e+002; krok(n+1)=2.343845e-003; ng(n+1)=4.438346e+001;
n=30; farx(n+1)=3.774032e+000; foe(n+1)=2.658552e+002; krok(n+1)=2.899547e-003; ng(n+1)=6.239052e+001;
n=31; farx(n+1)=3.624019e+000; foe(n+1)=2.121128e+002; krok(n+1)=1.571994e-002; ng(n+1)=9.228671e+001;
n=32; farx(n+1)=3.423311e+000; foe(n+1)=2.072698e+002; krok(n+1)=2.541206e-002; ng(n+1)=1.383566e+002;
n=33; farx(n+1)=3.309296e+000; foe(n+1)=2.106822e+002; krok(n+1)=3.766673e-002; ng(n+1)=1.162976e+002;
n=34; farx(n+1)=2.836954e+000; foe(n+1)=1.573291e+002; krok(n+1)=2.427388e-001; ng(n+1)=1.116313e+002;
n=35; farx(n+1)=2.660178e+000; foe(n+1)=1.400249e+002; krok(n+1)=1.036878e-001; ng(n+1)=9.935612e+001;
n=36; farx(n+1)=2.343891e+000; foe(n+1)=1.262007e+002; krok(n+1)=3.655097e-001; ng(n+1)=9.590247e+001;
n=37; farx(n+1)=1.991152e+000; foe(n+1)=1.935608e+002; krok(n+1)=7.763546e-001; ng(n+1)=9.882676e+001;
n=38; farx(n+1)=1.597981e+000; foe(n+1)=1.965305e+002; krok(n+1)=6.268073e-001; ng(n+1)=1.203020e+002;
n=39; farx(n+1)=1.432340e+000; foe(n+1)=1.488545e+002; krok(n+1)=3.329123e-001; ng(n+1)=6.114632e+001;
n=40; farx(n+1)=1.302000e+000; foe(n+1)=1.012102e+002; krok(n+1)=7.228765e-001; ng(n+1)=7.560763e+001;
n=41; farx(n+1)=1.210976e+000; foe(n+1)=1.310202e+002; krok(n+1)=7.847380e-001; ng(n+1)=3.648309e+001;
n=42; farx(n+1)=1.159613e+000; foe(n+1)=1.731850e+002; krok(n+1)=3.806077e-001; ng(n+1)=5.780837e+001;
n=43; farx(n+1)=1.118861e+000; foe(n+1)=1.687144e+002; krok(n+1)=2.688700e-001; ng(n+1)=1.614573e+001;
n=44; farx(n+1)=1.070999e+000; foe(n+1)=1.198370e+002; krok(n+1)=7.487569e-001; ng(n+1)=2.588785e+001;
n=45; farx(n+1)=1.001672e+000; foe(n+1)=6.739131e+001; krok(n+1)=7.807902e-001; ng(n+1)=3.309308e+001;
n=46; farx(n+1)=9.133940e-001; foe(n+1)=6.633694e+001; krok(n+1)=9.471081e-001; ng(n+1)=4.327698e+001;
n=47; farx(n+1)=8.675991e-001; foe(n+1)=5.250895e+001; krok(n+1)=6.212565e-001; ng(n+1)=3.344542e+001;
n=48; farx(n+1)=8.577653e-001; foe(n+1)=5.551103e+001; krok(n+1)=5.183452e-001; ng(n+1)=2.108484e+001;
n=49; farx(n+1)=8.533010e-001; foe(n+1)=5.565699e+001; krok(n+1)=1.328156e-001; ng(n+1)=2.256343e+001;
n=50; farx(n+1)=8.419132e-001; foe(n+1)=5.216721e+001; krok(n+1)=9.282331e-001; ng(n+1)=1.228804e+001;
%odnowa zmiennej metryki
n=51; farx(n+1)=8.399772e-001; foe(n+1)=5.356187e+001; krok(n+1)=2.439194e-004; ng(n+1)=9.746430e+000;
n=52; farx(n+1)=8.381058e-001; foe(n+1)=5.576291e+001; krok(n+1)=2.969451e-004; ng(n+1)=9.541826e+000;
n=53; farx(n+1)=8.376793e-001; foe(n+1)=5.355062e+001; krok(n+1)=1.711903e-004; ng(n+1)=4.824860e+000;
n=54; farx(n+1)=8.360735e-001; foe(n+1)=5.012869e+001; krok(n+1)=3.660956e-003; ng(n+1)=2.371612e+000;
n=55; farx(n+1)=8.181772e-001; foe(n+1)=5.332802e+001; krok(n+1)=8.172294e-002; ng(n+1)=1.673437e+000;
n=56; farx(n+1)=8.159208e-001; foe(n+1)=4.810886e+001; krok(n+1)=2.038266e-002; ng(n+1)=9.996463e+000;
n=57; farx(n+1)=8.133807e-001; foe(n+1)=5.077477e+001; krok(n+1)=2.080702e-002; ng(n+1)=1.231022e+001;
n=58; farx(n+1)=7.963724e-001; foe(n+1)=4.857884e+001; krok(n+1)=4.547265e-001; ng(n+1)=1.525164e+001;
n=59; farx(n+1)=7.879661e-001; foe(n+1)=4.852217e+001; krok(n+1)=2.833785e-001; ng(n+1)=1.924299e+001;
n=60; farx(n+1)=7.787587e-001; foe(n+1)=4.405021e+001; krok(n+1)=2.155114e-001; ng(n+1)=2.788771e+001;
n=61; farx(n+1)=7.694844e-001; foe(n+1)=4.819255e+001; krok(n+1)=1.902702e-001; ng(n+1)=2.022038e+001;
n=62; farx(n+1)=7.589995e-001; foe(n+1)=4.388118e+001; krok(n+1)=3.819012e-001; ng(n+1)=1.115189e+001;
n=63; farx(n+1)=7.497164e-001; foe(n+1)=4.217498e+001; krok(n+1)=2.775080e-001; ng(n+1)=1.150003e+001;
n=64; farx(n+1)=7.430128e-001; foe(n+1)=4.397019e+001; krok(n+1)=2.288809e-001; ng(n+1)=1.268971e+001;
n=65; farx(n+1)=7.330576e-001; foe(n+1)=4.235179e+001; krok(n+1)=6.121507e-001; ng(n+1)=7.028445e+000;
n=66; farx(n+1)=7.272384e-001; foe(n+1)=4.133765e+001; krok(n+1)=8.131859e-001; ng(n+1)=9.507737e+000;
n=67; farx(n+1)=7.254419e-001; foe(n+1)=4.246715e+001; krok(n+1)=3.200292e-001; ng(n+1)=9.352102e+000;
n=68; farx(n+1)=7.200961e-001; foe(n+1)=4.275863e+001; krok(n+1)=1.078288e+000; ng(n+1)=6.993244e+000;
n=69; farx(n+1)=7.163113e-001; foe(n+1)=4.472711e+001; krok(n+1)=1.111741e+000; ng(n+1)=4.382802e+000;
n=70; farx(n+1)=7.143405e-001; foe(n+1)=4.120256e+001; krok(n+1)=4.188878e-001; ng(n+1)=4.475234e+000;
n=71; farx(n+1)=7.127833e-001; foe(n+1)=4.260427e+001; krok(n+1)=9.673187e-001; ng(n+1)=2.943362e+000;
n=72; farx(n+1)=7.116175e-001; foe(n+1)=4.295197e+001; krok(n+1)=1.056270e+000; ng(n+1)=3.877678e+000;
n=73; farx(n+1)=7.112947e-001; foe(n+1)=4.201712e+001; krok(n+1)=1.015268e+000; ng(n+1)=3.225560e+000;
n=74; farx(n+1)=7.111910e-001; foe(n+1)=4.196665e+001; krok(n+1)=2.002310e+000; ng(n+1)=1.368815e+000;
n=75; farx(n+1)=7.111474e-001; foe(n+1)=4.203315e+001; krok(n+1)=1.361947e+000; ng(n+1)=1.510587e+000;
%odnowa zmiennej metryki
n=76; farx(n+1)=7.111409e-001; foe(n+1)=4.227260e+001; krok(n+1)=1.381179e-004; ng(n+1)=6.635283e-001;
n=77; farx(n+1)=7.111403e-001; foe(n+1)=4.229262e+001; krok(n+1)=2.108238e-004; ng(n+1)=1.916737e-001;
n=78; farx(n+1)=7.111392e-001; foe(n+1)=4.235271e+001; krok(n+1)=6.710424e-004; ng(n+1)=1.622897e-001;
n=79; farx(n+1)=7.111332e-001; foe(n+1)=4.272860e+001; krok(n+1)=1.321526e-002; ng(n+1)=6.890668e-002;
n=80; farx(n+1)=7.111047e-001; foe(n+1)=4.293641e+001; krok(n+1)=2.941687e-002; ng(n+1)=9.475980e-002;
n=81; farx(n+1)=7.110306e-001; foe(n+1)=4.229334e+001; krok(n+1)=3.542231e-002; ng(n+1)=1.613519e-001;
n=82; farx(n+1)=7.109840e-001; foe(n+1)=4.242175e+001; krok(n+1)=4.102628e-002; ng(n+1)=1.669667e-001;
n=83; farx(n+1)=7.108558e-001; foe(n+1)=4.269924e+001; krok(n+1)=5.343541e-001; ng(n+1)=2.591042e-001;
n=84; farx(n+1)=7.105850e-001; foe(n+1)=4.244122e+001; krok(n+1)=4.164005e-001; ng(n+1)=8.996361e-001;
n=85; farx(n+1)=7.103069e-001; foe(n+1)=4.199821e+001; krok(n+1)=3.240219e-001; ng(n+1)=3.188030e+000;
n=86; farx(n+1)=7.097635e-001; foe(n+1)=4.251646e+001; krok(n+1)=4.543926e-001; ng(n+1)=5.405554e+000;
n=87; farx(n+1)=7.094174e-001; foe(n+1)=4.296507e+001; krok(n+1)=4.268310e-001; ng(n+1)=4.979671e+000;
n=88; farx(n+1)=7.090812e-001; foe(n+1)=4.294424e+001; krok(n+1)=8.728449e-001; ng(n+1)=2.087576e+000;
n=89; farx(n+1)=7.087133e-001; foe(n+1)=4.241707e+001; krok(n+1)=3.745525e-001; ng(n+1)=1.162306e+000;
n=90; farx(n+1)=7.080100e-001; foe(n+1)=4.197421e+001; krok(n+1)=9.001093e-001; ng(n+1)=5.091357e+000;
n=91; farx(n+1)=7.077737e-001; foe(n+1)=4.336544e+001; krok(n+1)=2.936023e-001; ng(n+1)=5.726864e+000;
n=92; farx(n+1)=7.073594e-001; foe(n+1)=4.371937e+001; krok(n+1)=4.801921e-001; ng(n+1)=4.876693e+000;
n=93; farx(n+1)=7.068489e-001; foe(n+1)=4.282348e+001; krok(n+1)=5.231840e-001; ng(n+1)=6.794249e-001;
n=94; farx(n+1)=7.064718e-001; foe(n+1)=4.185904e+001; krok(n+1)=6.665509e-001; ng(n+1)=3.400174e+000;
n=95; farx(n+1)=7.060386e-001; foe(n+1)=4.211880e+001; krok(n+1)=9.372047e-001; ng(n+1)=5.286968e+000;
n=96; farx(n+1)=7.055049e-001; foe(n+1)=4.272087e+001; krok(n+1)=9.073805e-001; ng(n+1)=1.286990e+000;
n=97; farx(n+1)=7.051410e-001; foe(n+1)=4.341925e+001; krok(n+1)=5.189311e-001; ng(n+1)=4.644736e+000;
n=98; farx(n+1)=7.046773e-001; foe(n+1)=4.324744e+001; krok(n+1)=1.273593e+000; ng(n+1)=4.977279e+000;
n=99; farx(n+1)=7.044011e-001; foe(n+1)=4.297254e+001; krok(n+1)=4.268310e-001; ng(n+1)=1.466911e+000;
n=100; farx(n+1)=7.040216e-001; foe(n+1)=4.374602e+001; krok(n+1)=8.773192e-001; ng(n+1)=3.584518e+000;
%odnowa zmiennej metryki
n=101; farx(n+1)=7.038532e-001; foe(n+1)=4.328184e+001; krok(n+1)=1.323954e-004; ng(n+1)=3.978300e+000;
n=102; farx(n+1)=7.038465e-001; foe(n+1)=4.351191e+001; krok(n+1)=1.443891e-004; ng(n+1)=7.117862e-001;
n=103; farx(n+1)=7.038252e-001; foe(n+1)=4.364495e+001; krok(n+1)=5.173682e-003; ng(n+1)=2.427628e-001;
n=104; farx(n+1)=7.038125e-001; foe(n+1)=4.400453e+001; krok(n+1)=1.181843e-003; ng(n+1)=3.983537e-001;
n=105; farx(n+1)=7.037966e-001; foe(n+1)=4.442386e+001; krok(n+1)=1.600678e-002; ng(n+1)=1.115898e-001;
n=106; farx(n+1)=7.037637e-001; foe(n+1)=4.375179e+001; krok(n+1)=1.048565e-002; ng(n+1)=1.777272e-001;
n=107; farx(n+1)=7.035575e-001; foe(n+1)=4.348237e+001; krok(n+1)=1.926134e-001; ng(n+1)=1.189746e-001;
n=108; farx(n+1)=7.033790e-001; foe(n+1)=4.327129e+001; krok(n+1)=3.589369e-001; ng(n+1)=3.852564e-001;
n=109; farx(n+1)=7.032147e-001; foe(n+1)=4.305064e+001; krok(n+1)=4.770341e-001; ng(n+1)=1.169630e+000;
n=110; farx(n+1)=7.030673e-001; foe(n+1)=4.258636e+001; krok(n+1)=3.005969e-001; ng(n+1)=2.099358e+000;
n=111; farx(n+1)=7.028868e-001; foe(n+1)=4.287587e+001; krok(n+1)=5.362494e-001; ng(n+1)=3.045389e+000;
n=112; farx(n+1)=7.026747e-001; foe(n+1)=4.367610e+001; krok(n+1)=4.461224e-001; ng(n+1)=1.648199e+000;
n=113; farx(n+1)=7.023350e-001; foe(n+1)=4.286891e+001; krok(n+1)=6.277333e-001; ng(n+1)=7.343121e-001;
n=114; farx(n+1)=7.020069e-001; foe(n+1)=4.371666e+001; krok(n+1)=7.970783e-001; ng(n+1)=3.240476e+000;
n=115; farx(n+1)=7.017016e-001; foe(n+1)=4.505133e+001; krok(n+1)=1.295594e+000; ng(n+1)=4.604052e+000;
n=116; farx(n+1)=7.011802e-001; foe(n+1)=4.336835e+001; krok(n+1)=7.382907e-001; ng(n+1)=1.391252e+000;
n=117; farx(n+1)=7.008499e-001; foe(n+1)=4.324397e+001; krok(n+1)=5.403454e-001; ng(n+1)=4.366352e+000;
n=118; farx(n+1)=7.003769e-001; foe(n+1)=4.349944e+001; krok(n+1)=5.758050e-001; ng(n+1)=6.921741e+000;
n=119; farx(n+1)=6.997282e-001; foe(n+1)=4.359205e+001; krok(n+1)=1.244898e+000; ng(n+1)=2.751283e+000;
n=120; farx(n+1)=6.994997e-001; foe(n+1)=4.221747e+001; krok(n+1)=4.686024e-001; ng(n+1)=4.758942e+000;
n=121; farx(n+1)=6.988316e-001; foe(n+1)=4.304521e+001; krok(n+1)=7.273018e-001; ng(n+1)=6.394708e+000;
n=122; farx(n+1)=6.981938e-001; foe(n+1)=4.265227e+001; krok(n+1)=5.491644e-001; ng(n+1)=1.300680e+000;
n=123; farx(n+1)=6.974188e-001; foe(n+1)=4.222815e+001; krok(n+1)=7.310195e-001; ng(n+1)=5.630423e+000;
n=124; farx(n+1)=6.960012e-001; foe(n+1)=4.375658e+001; krok(n+1)=9.512526e-001; ng(n+1)=8.408940e+000;
n=125; farx(n+1)=6.948875e-001; foe(n+1)=4.128403e+001; krok(n+1)=6.553300e-001; ng(n+1)=4.228268e+000;
%odnowa zmiennej metryki
n=126; farx(n+1)=6.940137e-001; foe(n+1)=4.359770e+001; krok(n+1)=1.231773e-004; ng(n+1)=9.830538e+000;
n=127; farx(n+1)=6.939899e-001; foe(n+1)=4.344547e+001; krok(n+1)=1.355338e-004; ng(n+1)=1.655450e+000;
n=128; farx(n+1)=6.939093e-001; foe(n+1)=4.420779e+001; krok(n+1)=2.561707e-003; ng(n+1)=6.176388e-001;
n=129; farx(n+1)=6.938659e-001; foe(n+1)=4.506042e+001; krok(n+1)=1.345435e-003; ng(n+1)=6.881066e-001;
n=130; farx(n+1)=6.936344e-001; foe(n+1)=4.418544e+001; krok(n+1)=1.356166e-002; ng(n+1)=4.987288e-001;
n=131; farx(n+1)=6.934957e-001; foe(n+1)=4.506604e+001; krok(n+1)=4.095813e-002; ng(n+1)=1.948972e-001;
n=132; farx(n+1)=6.928492e-001; foe(n+1)=4.378367e+001; krok(n+1)=5.005537e-002; ng(n+1)=4.197671e-001;
n=133; farx(n+1)=6.922970e-001; foe(n+1)=4.306665e+001; krok(n+1)=2.857025e-001; ng(n+1)=8.919139e-001;
n=134; farx(n+1)=6.916847e-001; foe(n+1)=4.257855e+001; krok(n+1)=3.923690e-001; ng(n+1)=2.237059e+000;
n=135; farx(n+1)=6.908102e-001; foe(n+1)=4.182560e+001; krok(n+1)=4.079802e-001; ng(n+1)=4.307004e+000;
n=136; farx(n+1)=6.890377e-001; foe(n+1)=4.385367e+001; krok(n+1)=5.391440e-001; ng(n+1)=1.031969e+001;
n=137; farx(n+1)=6.863895e-001; foe(n+1)=4.439758e+001; krok(n+1)=1.003200e+000; ng(n+1)=1.272839e+001;
n=138; farx(n+1)=6.847004e-001; foe(n+1)=4.081906e+001; krok(n+1)=2.711458e-001; ng(n+1)=1.739840e+000;
n=139; farx(n+1)=6.821231e-001; foe(n+1)=4.040269e+001; krok(n+1)=2.250273e-001; ng(n+1)=8.393045e+000;
n=140; farx(n+1)=6.766510e-001; foe(n+1)=4.507570e+001; krok(n+1)=6.071167e-001; ng(n+1)=1.866799e+001;
n=141; farx(n+1)=6.705496e-001; foe(n+1)=3.935925e+001; krok(n+1)=3.706964e-001; ng(n+1)=2.028406e+001;
n=142; farx(n+1)=6.592216e-001; foe(n+1)=4.330022e+001; krok(n+1)=4.543926e-001; ng(n+1)=2.061804e+001;
n=143; farx(n+1)=6.477577e-001; foe(n+1)=3.787776e+001; krok(n+1)=6.400584e-001; ng(n+1)=8.268035e+000;
n=144; farx(n+1)=6.320807e-001; foe(n+1)=3.643852e+001; krok(n+1)=2.044970e-001; ng(n+1)=2.937659e+001;
n=145; farx(n+1)=6.146265e-001; foe(n+1)=2.498341e+001; krok(n+1)=9.790078e-002; ng(n+1)=2.516258e+001;
n=146; farx(n+1)=6.096629e-001; foe(n+1)=2.587863e+001; krok(n+1)=1.022485e-001; ng(n+1)=9.252824e+000;
n=147; farx(n+1)=5.966685e-001; foe(n+1)=2.306224e+001; krok(n+1)=2.230542e-001; ng(n+1)=1.034889e+001;
n=148; farx(n+1)=5.827814e-001; foe(n+1)=2.035131e+001; krok(n+1)=3.691454e-001; ng(n+1)=2.141448e+001;
n=149; farx(n+1)=5.642104e-001; foe(n+1)=2.368463e+001; krok(n+1)=6.329112e-001; ng(n+1)=1.770407e+001;
n=150; farx(n+1)=5.495703e-001; foe(n+1)=1.500241e+001; krok(n+1)=6.155707e-001; ng(n+1)=1.333798e+001;
%odnowa zmiennej metryki
n=151; farx(n+1)=5.454667e-001; foe(n+1)=1.509195e+001; krok(n+1)=7.818173e-005; ng(n+1)=2.726986e+001;
n=152; farx(n+1)=5.452368e-001; foe(n+1)=1.542272e+001; krok(n+1)=1.554679e-004; ng(n+1)=4.437353e+000;
n=153; farx(n+1)=5.443374e-001; foe(n+1)=1.367883e+001; krok(n+1)=1.100811e-003; ng(n+1)=3.396472e+000;
n=154; farx(n+1)=5.439473e-001; foe(n+1)=1.322940e+001; krok(n+1)=1.025539e-003; ng(n+1)=2.334034e+000;
n=155; farx(n+1)=5.414624e-001; foe(n+1)=1.289531e+001; krok(n+1)=2.812841e-002; ng(n+1)=1.020884e+000;
n=156; farx(n+1)=5.405884e-001; foe(n+1)=1.228131e+001; krok(n+1)=4.904612e-002; ng(n+1)=5.286109e-001;
n=157; farx(n+1)=5.398789e-001; foe(n+1)=1.217814e+001; krok(n+1)=3.696362e-002; ng(n+1)=7.491873e-001;
n=158; farx(n+1)=5.370825e-001; foe(n+1)=1.085007e+001; krok(n+1)=2.250273e-001; ng(n+1)=1.118000e+000;
n=159; farx(n+1)=5.356348e-001; foe(n+1)=9.652361e+000; krok(n+1)=2.926858e-001; ng(n+1)=9.039231e+000;
n=160; farx(n+1)=5.343590e-001; foe(n+1)=9.315643e+000; krok(n+1)=4.500546e-001; ng(n+1)=1.174770e+001;
n=161; farx(n+1)=5.322001e-001; foe(n+1)=9.402766e+000; krok(n+1)=9.001093e-001; ng(n+1)=9.499628e+000;
n=162; farx(n+1)=5.304416e-001; foe(n+1)=1.008514e+001; krok(n+1)=2.957090e-001; ng(n+1)=1.088430e+001;
n=163; farx(n+1)=5.286164e-001; foe(n+1)=9.077710e+000; krok(n+1)=7.783923e-001; ng(n+1)=2.538167e+000;
n=164; farx(n+1)=5.270025e-001; foe(n+1)=8.328752e+000; krok(n+1)=5.312625e-001; ng(n+1)=7.246489e+000;
n=165; farx(n+1)=5.257964e-001; foe(n+1)=8.225335e+000; krok(n+1)=2.401907e-001; ng(n+1)=9.342151e+000;
n=166; farx(n+1)=5.246901e-001; foe(n+1)=7.670456e+000; krok(n+1)=3.749470e-001; ng(n+1)=5.864865e+000;
n=167; farx(n+1)=5.234589e-001; foe(n+1)=8.005420e+000; krok(n+1)=5.231840e-001; ng(n+1)=1.045268e+001;
n=168; farx(n+1)=5.227342e-001; foe(n+1)=8.403505e+000; krok(n+1)=9.830700e-001; ng(n+1)=4.385069e+000;
n=169; farx(n+1)=5.221939e-001; foe(n+1)=8.213635e+000; krok(n+1)=1.885812e+000; ng(n+1)=2.913630e+000;
n=170; farx(n+1)=5.219226e-001; foe(n+1)=8.538322e+000; krok(n+1)=6.871474e-001; ng(n+1)=3.645497e+000;
n=171; farx(n+1)=5.214924e-001; foe(n+1)=8.156572e+000; krok(n+1)=8.295027e-001; ng(n+1)=3.067782e+000;
n=172; farx(n+1)=5.212073e-001; foe(n+1)=8.350493e+000; krok(n+1)=1.086288e+000; ng(n+1)=1.796975e+000;
n=173; farx(n+1)=5.207421e-001; foe(n+1)=8.147401e+000; krok(n+1)=9.853052e-001; ng(n+1)=3.596038e+000;
n=174; farx(n+1)=5.204627e-001; foe(n+1)=8.392095e+000; krok(n+1)=8.812470e-001; ng(n+1)=1.603611e+000;
n=175; farx(n+1)=5.202161e-001; foe(n+1)=8.570882e+000; krok(n+1)=9.453922e-001; ng(n+1)=4.271345e+000;
%odnowa zmiennej metryki
n=176; farx(n+1)=5.201912e-001; foe(n+1)=8.436470e+000; krok(n+1)=1.024331e-004; ng(n+1)=1.812322e+000;
n=177; farx(n+1)=5.201759e-001; foe(n+1)=8.398818e+000; krok(n+1)=5.873414e-004; ng(n+1)=6.431598e-001;
n=178; farx(n+1)=5.201635e-001; foe(n+1)=8.389608e+000; krok(n+1)=1.251495e-004; ng(n+1)=1.201523e+000;
n=179; farx(n+1)=5.201475e-001; foe(n+1)=8.229019e+000; krok(n+1)=1.259039e-003; ng(n+1)=4.431315e-001;
n=180; farx(n+1)=5.201123e-001; foe(n+1)=8.334511e+000; krok(n+1)=1.121665e-002; ng(n+1)=1.923150e-001;
n=181; farx(n+1)=5.199797e-001; foe(n+1)=8.239148e+000; krok(n+1)=1.827549e-001; ng(n+1)=9.816864e-002;
n=182; farx(n+1)=5.199239e-001; foe(n+1)=8.217536e+000; krok(n+1)=4.509055e-002; ng(n+1)=2.876959e-001;
n=183; farx(n+1)=5.197430e-001; foe(n+1)=7.996929e+000; krok(n+1)=3.381117e-001; ng(n+1)=7.218528e-001;
n=184; farx(n+1)=5.195180e-001; foe(n+1)=7.821477e+000; krok(n+1)=1.186316e-001; ng(n+1)=1.951307e+000;
n=185; farx(n+1)=5.190889e-001; foe(n+1)=7.610672e+000; krok(n+1)=8.773192e-001; ng(n+1)=5.018697e+000;
n=186; farx(n+1)=5.186269e-001; foe(n+1)=7.622379e+000; krok(n+1)=7.011953e-001; ng(n+1)=7.085105e+000;
n=187; farx(n+1)=5.182804e-001; foe(n+1)=8.167645e+000; krok(n+1)=4.079802e-001; ng(n+1)=6.491344e+000;
n=188; farx(n+1)=5.173276e-001; foe(n+1)=8.439942e+000; krok(n+1)=5.856145e-001; ng(n+1)=4.640953e+000;
n=189; farx(n+1)=5.169037e-001; foe(n+1)=7.830550e+000; krok(n+1)=2.676098e-001; ng(n+1)=5.333789e+000;
n=190; farx(n+1)=5.160788e-001; foe(n+1)=7.264136e+000; krok(n+1)=6.643777e-001; ng(n+1)=9.653817e+000;
n=191; farx(n+1)=5.153263e-001; foe(n+1)=7.238397e+000; krok(n+1)=6.674447e-001; ng(n+1)=4.618072e+000;
n=192; farx(n+1)=5.147191e-001; foe(n+1)=7.970376e+000; krok(n+1)=5.695627e-001; ng(n+1)=6.490344e+000;
n=193; farx(n+1)=5.143364e-001; foe(n+1)=7.700131e+000; krok(n+1)=1.162023e+000; ng(n+1)=3.655162e+000;
n=194; farx(n+1)=5.141077e-001; foe(n+1)=7.868786e+000; krok(n+1)=7.273018e-001; ng(n+1)=6.058519e+000;
n=195; farx(n+1)=5.137837e-001; foe(n+1)=7.435577e+000; krok(n+1)=3.262230e-001; ng(n+1)=2.404015e+000;
n=196; farx(n+1)=5.134121e-001; foe(n+1)=7.192831e+000; krok(n+1)=1.112649e+000; ng(n+1)=4.300386e+000;
n=197; farx(n+1)=5.132419e-001; foe(n+1)=7.158546e+000; krok(n+1)=7.382907e-001; ng(n+1)=4.221042e+000;
n=198; farx(n+1)=5.131672e-001; foe(n+1)=7.407167e+000; krok(n+1)=2.320583e-001; ng(n+1)=1.721574e+000;
n=199; farx(n+1)=5.130638e-001; foe(n+1)=7.223840e+000; krok(n+1)=1.056270e+000; ng(n+1)=1.258109e+000;
n=200; farx(n+1)=5.127476e-001; foe(n+1)=7.368880e+000; krok(n+1)=3.644931e+000; ng(n+1)=2.935268e+000;
%odnowa zmiennej metryki
n=201; farx(n+1)=5.127116e-001; foe(n+1)=7.351484e+000; krok(n+1)=7.524177e-005; ng(n+1)=2.359672e+000;
n=202; farx(n+1)=5.126704e-001; foe(n+1)=7.323535e+000; krok(n+1)=2.302017e-004; ng(n+1)=1.626287e+000;
n=203; farx(n+1)=5.126463e-001; foe(n+1)=7.450815e+000; krok(n+1)=1.530291e-004; ng(n+1)=1.379681e+000;
n=204; farx(n+1)=5.126385e-001; foe(n+1)=7.526309e+000; krok(n+1)=3.904632e-003; ng(n+1)=1.558841e-001;
n=205; farx(n+1)=5.125924e-001; foe(n+1)=7.365984e+000; krok(n+1)=2.147336e-002; ng(n+1)=1.484650e-001;
n=206; farx(n+1)=5.125281e-001; foe(n+1)=7.330082e+000; krok(n+1)=3.411686e-002; ng(n+1)=1.601811e-001;
n=207; farx(n+1)=5.124837e-001; foe(n+1)=7.276513e+000; krok(n+1)=4.438543e-002; ng(n+1)=1.219717e-001;
n=208; farx(n+1)=5.124426e-001; foe(n+1)=7.140678e+000; krok(n+1)=7.301859e-002; ng(n+1)=2.701789e-001;
n=209; farx(n+1)=5.122880e-001; foe(n+1)=6.772918e+000; krok(n+1)=4.376495e-001; ng(n+1)=3.994800e-001;
n=210; farx(n+1)=5.120068e-001; foe(n+1)=6.661695e+000; krok(n+1)=1.289896e+000; ng(n+1)=1.757343e+000;
n=211; farx(n+1)=5.118060e-001; foe(n+1)=6.814159e+000; krok(n+1)=3.372291e-001; ng(n+1)=5.079796e+000;
n=212; farx(n+1)=5.114128e-001; foe(n+1)=7.441344e+000; krok(n+1)=8.162976e-001; ng(n+1)=6.594825e+000;
n=213; farx(n+1)=5.110981e-001; foe(n+1)=7.363928e+000; krok(n+1)=1.304892e+000; ng(n+1)=1.554371e+000;
n=214; farx(n+1)=5.109542e-001; foe(n+1)=7.216401e+000; krok(n+1)=4.862877e-001; ng(n+1)=2.972533e+000;
n=215; farx(n+1)=5.108591e-001; foe(n+1)=7.341167e+000; krok(n+1)=3.214238e-001; ng(n+1)=4.880995e+000;
n=216; farx(n+1)=5.106747e-001; foe(n+1)=7.295285e+000; krok(n+1)=4.128102e-001; ng(n+1)=2.883390e+000;
n=217; farx(n+1)=5.105595e-001; foe(n+1)=7.105258e+000; krok(n+1)=9.685429e-001; ng(n+1)=3.117078e+000;
n=218; farx(n+1)=5.103772e-001; foe(n+1)=7.554219e+000; krok(n+1)=9.639418e-001; ng(n+1)=2.223529e+000;
n=219; farx(n+1)=5.100843e-001; foe(n+1)=7.805684e+000; krok(n+1)=1.707324e+000; ng(n+1)=3.587295e+000;
n=220; farx(n+1)=5.097481e-001; foe(n+1)=7.403838e+000; krok(n+1)=2.640676e-001; ng(n+1)=3.416412e+000;
n=221; farx(n+1)=5.094980e-001; foe(n+1)=7.043919e+000; krok(n+1)=4.884576e-001; ng(n+1)=6.992510e+000;
n=222; farx(n+1)=5.090328e-001; foe(n+1)=6.990768e+000; krok(n+1)=1.000418e+000; ng(n+1)=7.358727e+000;
n=223; farx(n+1)=5.087922e-001; foe(n+1)=6.994273e+000; krok(n+1)=4.322425e-001; ng(n+1)=2.945530e+000;
n=224; farx(n+1)=5.085620e-001; foe(n+1)=7.022215e+000; krok(n+1)=1.280117e+000; ng(n+1)=4.709089e+000;
n=225; farx(n+1)=5.083951e-001; foe(n+1)=7.248455e+000; krok(n+1)=1.507294e+000; ng(n+1)=5.376699e+000;
%odnowa zmiennej metryki
n=226; farx(n+1)=5.083941e-001; foe(n+1)=7.253176e+000; krok(n+1)=1.395201e-004; ng(n+1)=3.173733e-001;
n=227; farx(n+1)=5.083862e-001; foe(n+1)=7.342259e+000; krok(n+1)=4.526501e-004; ng(n+1)=4.749512e-001;
n=228; farx(n+1)=5.083819e-001; foe(n+1)=7.376623e+000; krok(n+1)=1.291963e-003; ng(n+1)=2.065547e-001;
n=229; farx(n+1)=5.083811e-001; foe(n+1)=7.360280e+000; krok(n+1)=2.626087e-004; ng(n+1)=2.064302e-001;
n=230; farx(n+1)=5.083612e-001; foe(n+1)=7.394765e+000; krok(n+1)=3.499183e-003; ng(n+1)=2.974638e-001;
n=231; farx(n+1)=5.083126e-001; foe(n+1)=7.317387e+000; krok(n+1)=6.173004e-002; ng(n+1)=1.056085e-001;
n=232; farx(n+1)=5.082395e-001; foe(n+1)=7.115171e+000; krok(n+1)=4.486711e-002; ng(n+1)=2.228561e-001;
n=233; farx(n+1)=5.081570e-001; foe(n+1)=6.936456e+000; krok(n+1)=3.655097e-001; ng(n+1)=4.154590e-001;
n=234; farx(n+1)=5.081301e-001; foe(n+1)=6.893441e+000; krok(n+1)=8.268560e-002; ng(n+1)=1.004332e+000;
n=235; farx(n+1)=5.080133e-001; foe(n+1)=6.745497e+000; krok(n+1)=4.756263e-001; ng(n+1)=1.215745e+000;
n=236; farx(n+1)=5.078819e-001; foe(n+1)=6.763685e+000; krok(n+1)=5.312625e-001; ng(n+1)=3.006060e+000;
n=237; farx(n+1)=5.077381e-001; foe(n+1)=6.868338e+000; krok(n+1)=4.555491e-001; ng(n+1)=4.970208e+000;
n=238; farx(n+1)=5.075702e-001; foe(n+1)=7.128348e+000; krok(n+1)=5.945315e-001; ng(n+1)=4.985398e+000;
n=239; farx(n+1)=5.073707e-001; foe(n+1)=6.886734e+000; krok(n+1)=1.065517e+000; ng(n+1)=2.640812e+000;
n=240; farx(n+1)=5.072021e-001; foe(n+1)=6.888401e+000; krok(n+1)=6.002260e-001; ng(n+1)=4.563780e+000;
n=241; farx(n+1)=5.069153e-001; foe(n+1)=6.916635e+000; krok(n+1)=1.691553e+000; ng(n+1)=4.246306e+000;
n=242; farx(n+1)=5.068196e-001; foe(n+1)=6.970828e+000; krok(n+1)=4.854776e-001; ng(n+1)=3.177733e+000;
n=243; farx(n+1)=5.065711e-001; foe(n+1)=6.806482e+000; krok(n+1)=6.121507e-001; ng(n+1)=4.088255e+000;
n=244; farx(n+1)=5.064186e-001; foe(n+1)=7.074069e+000; krok(n+1)=4.268310e-001; ng(n+1)=3.438556e+000;
n=245; farx(n+1)=5.060738e-001; foe(n+1)=7.363308e+000; krok(n+1)=2.196657e+000; ng(n+1)=5.167155e+000;
n=246; farx(n+1)=5.057933e-001; foe(n+1)=7.019873e+000; krok(n+1)=2.818076e-001; ng(n+1)=2.427820e+000;
n=247; farx(n+1)=5.055256e-001; foe(n+1)=6.793974e+000; krok(n+1)=8.922169e-001; ng(n+1)=3.631241e+000;
n=248; farx(n+1)=5.051333e-001; foe(n+1)=6.870503e+000; krok(n+1)=7.847380e-001; ng(n+1)=7.519817e+000;
n=249; farx(n+1)=5.048624e-001; foe(n+1)=6.959874e+000; krok(n+1)=7.238398e-001; ng(n+1)=1.505992e+000;
n=250; farx(n+1)=5.045111e-001; foe(n+1)=6.879083e+000; krok(n+1)=6.780158e-001; ng(n+1)=7.277003e+000;
%odnowa zmiennej metryki
n=251; farx(n+1)=5.044243e-001; foe(n+1)=6.960189e+000; krok(n+1)=3.639457e-005; ng(n+1)=5.658013e+000;
n=252; farx(n+1)=5.044065e-001; foe(n+1)=7.020469e+000; krok(n+1)=7.813213e-005; ng(n+1)=1.837641e+000;
n=253; farx(n+1)=5.043519e-001; foe(n+1)=6.982728e+000; krok(n+1)=4.375565e-003; ng(n+1)=4.158070e-001;
n=254; farx(n+1)=5.043105e-001; foe(n+1)=7.013741e+000; krok(n+1)=1.564230e-003; ng(n+1)=5.316703e-001;
n=255; farx(n+1)=5.042836e-001; foe(n+1)=7.175889e+000; krok(n+1)=1.579878e-003; ng(n+1)=4.710088e-001;
n=256; farx(n+1)=5.042557e-001; foe(n+1)=7.165862e+000; krok(n+1)=3.084406e-002; ng(n+1)=1.475066e-001;
n=257; farx(n+1)=5.041434e-001; foe(n+1)=6.906459e+000; krok(n+1)=1.183885e-001; ng(n+1)=2.254259e-001;
n=258; farx(n+1)=5.040758e-001; foe(n+1)=6.778399e+000; krok(n+1)=6.864555e-002; ng(n+1)=3.738114e-001;
n=259; farx(n+1)=5.038619e-001; foe(n+1)=6.442708e+000; krok(n+1)=4.041502e-001; ng(n+1)=7.693327e-001;
n=260; farx(n+1)=5.036967e-001; foe(n+1)=6.303644e+000; krok(n+1)=3.354651e-001; ng(n+1)=2.757771e+000;
n=261; farx(n+1)=5.035197e-001; foe(n+1)=6.463627e+000; krok(n+1)=3.214238e-001; ng(n+1)=5.739884e+000;
n=262; farx(n+1)=5.031901e-001; foe(n+1)=6.794579e+000; krok(n+1)=7.926233e-001; ng(n+1)=8.181324e+000;
n=263; farx(n+1)=5.029741e-001; foe(n+1)=6.796192e+000; krok(n+1)=7.215816e-001; ng(n+1)=3.694578e+000;
n=264; farx(n+1)=5.027396e-001; foe(n+1)=6.528030e+000; krok(n+1)=5.792353e-001; ng(n+1)=3.178876e+000;
n=265; farx(n+1)=5.025038e-001; foe(n+1)=6.548522e+000; krok(n+1)=5.727234e-001; ng(n+1)=6.540208e+000;
n=266; farx(n+1)=5.023513e-001; foe(n+1)=6.531629e+000; krok(n+1)=6.155707e-001; ng(n+1)=5.753548e+000;
n=267; farx(n+1)=5.020989e-001; foe(n+1)=6.364234e+000; krok(n+1)=9.270274e-001; ng(n+1)=2.741520e+000;
n=268; farx(n+1)=5.019292e-001; foe(n+1)=6.553524e+000; krok(n+1)=4.863560e-001; ng(n+1)=5.358202e+000;
n=269; farx(n+1)=5.016812e-001; foe(n+1)=6.630748e+000; krok(n+1)=8.620454e-001; ng(n+1)=3.109045e+000;
n=270; farx(n+1)=5.014322e-001; foe(n+1)=6.642374e+000; krok(n+1)=1.221873e+000; ng(n+1)=4.519314e+000;
n=271; farx(n+1)=5.012288e-001; foe(n+1)=6.918703e+000; krok(n+1)=6.622313e-001; ng(n+1)=6.434443e+000;
n=272; farx(n+1)=5.006008e-001; foe(n+1)=6.442115e+000; krok(n+1)=1.236644e+000; ng(n+1)=1.620613e+000;
n=273; farx(n+1)=4.999564e-001; foe(n+1)=6.104451e+000; krok(n+1)=1.342163e+000; ng(n+1)=7.017469e+000;
n=274; farx(n+1)=4.995500e-001; foe(n+1)=6.378908e+000; krok(n+1)=4.329486e-001; ng(n+1)=7.183956e+000;
n=275; farx(n+1)=4.986541e-001; foe(n+1)=6.517767e+000; krok(n+1)=4.228884e-001; ng(n+1)=6.170742e+000;
%odnowa zmiennej metryki
n=276; farx(n+1)=4.983443e-001; foe(n+1)=6.577246e+000; krok(n+1)=2.394830e-005; ng(n+1)=1.133616e+001;
n=277; farx(n+1)=4.982909e-001; foe(n+1)=6.641916e+000; krok(n+1)=5.716520e-004; ng(n+1)=1.220214e+000;
n=278; farx(n+1)=4.982792e-001; foe(n+1)=6.594530e+000; krok(n+1)=1.070946e-004; ng(n+1)=1.245314e+000;
n=279; farx(n+1)=4.982251e-001; foe(n+1)=6.546365e+000; krok(n+1)=5.282995e-003; ng(n+1)=3.314182e-001;
n=280; farx(n+1)=4.981857e-001; foe(n+1)=6.723395e+000; krok(n+1)=2.385598e-003; ng(n+1)=4.770488e-001;
n=281; farx(n+1)=4.980850e-001; foe(n+1)=6.601827e+000; krok(n+1)=7.233428e-002; ng(n+1)=1.844527e-001;
n=282; farx(n+1)=4.980052e-001; foe(n+1)=6.598921e+000; krok(n+1)=5.671128e-002; ng(n+1)=2.058571e-001;
n=283; farx(n+1)=4.979892e-001; foe(n+1)=6.519587e+000; krok(n+1)=2.355148e-002; ng(n+1)=4.673470e-001;
n=284; farx(n+1)=4.979335e-001; foe(n+1)=6.327653e+000; krok(n+1)=1.961845e-001; ng(n+1)=4.375053e-001;
n=285; farx(n+1)=4.978264e-001; foe(n+1)=6.152896e+000; krok(n+1)=3.337224e-001; ng(n+1)=4.007296e-001;
n=286; farx(n+1)=4.975501e-001; foe(n+1)=5.884602e+000; krok(n+1)=6.058278e-001; ng(n+1)=9.094632e-001;
n=287; farx(n+1)=4.970794e-001; foe(n+1)=5.920752e+000; krok(n+1)=8.181022e-001; ng(n+1)=5.265052e+000;
n=288; farx(n+1)=4.966994e-001; foe(n+1)=6.198902e+000; krok(n+1)=5.596492e-001; ng(n+1)=1.011860e+001;
n=289; farx(n+1)=4.959512e-001; foe(n+1)=6.069395e+000; krok(n+1)=1.449678e+000; ng(n+1)=5.598562e+000;
n=290; farx(n+1)=4.953900e-001; foe(n+1)=5.920906e+000; krok(n+1)=1.040494e+000; ng(n+1)=7.887395e+000;
n=291; farx(n+1)=4.951482e-001; foe(n+1)=6.063517e+000; krok(n+1)=1.655578e-001; ng(n+1)=8.524624e+000;
n=292; farx(n+1)=4.940247e-001; foe(n+1)=5.825916e+000; krok(n+1)=1.505055e+000; ng(n+1)=5.734324e+000;
n=293; farx(n+1)=4.931011e-001; foe(n+1)=5.998449e+000; krok(n+1)=6.744582e-001; ng(n+1)=4.841309e+000;
n=294; farx(n+1)=4.923559e-001; foe(n+1)=6.393890e+000; krok(n+1)=2.701916e-001; ng(n+1)=8.383384e+000;
n=295; farx(n+1)=4.919369e-001; foe(n+1)=6.367236e+000; krok(n+1)=3.304177e-001; ng(n+1)=6.105625e+000;
n=296; farx(n+1)=4.913590e-001; foe(n+1)=6.098676e+000; krok(n+1)=4.014974e-001; ng(n+1)=4.600473e+000;
n=297; farx(n+1)=4.907384e-001; foe(n+1)=5.681004e+000; krok(n+1)=8.773192e-001; ng(n+1)=2.926955e+000;
n=298; farx(n+1)=4.903648e-001; foe(n+1)=5.449556e+000; krok(n+1)=6.614848e-001; ng(n+1)=4.171120e+000;
n=299; farx(n+1)=4.900255e-001; foe(n+1)=5.414316e+000; krok(n+1)=1.068393e+000; ng(n+1)=5.157042e+000;
n=300; farx(n+1)=4.895816e-001; foe(n+1)=5.488045e+000; krok(n+1)=8.127796e-001; ng(n+1)=3.994729e+000;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora ARX');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
