%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.475612e+003; foe(n+1)=4.571199e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=3.761043e+003; foe(n+1)=3.907792e+003; krok(n+1)=5.053235e-004; ng(n+1)=3.945823e+003;
n=2; farx(n+1)=9.276067e+002; foe(n+1)=1.177501e+003; krok(n+1)=3.065383e-003; ng(n+1)=2.338939e+003;
n=3; farx(n+1)=5.171372e+002; foe(n+1)=8.623928e+002; krok(n+1)=1.120093e-004; ng(n+1)=7.463035e+003;
n=4; farx(n+1)=5.412733e+002; foe(n+1)=8.288214e+002; krok(n+1)=2.034749e-004; ng(n+1)=2.136019e+003;
n=5; farx(n+1)=4.673794e+002; foe(n+1)=8.057553e+002; krok(n+1)=1.070946e-004; ng(n+1)=2.391647e+003;
n=6; farx(n+1)=4.881077e+002; foe(n+1)=7.868236e+002; krok(n+1)=1.569452e-004; ng(n+1)=1.506832e+003;
n=7; farx(n+1)=4.275667e+002; foe(n+1)=7.673818e+002; krok(n+1)=1.007837e-004; ng(n+1)=2.066015e+003;
n=8; farx(n+1)=4.449175e+002; foe(n+1)=7.494870e+002; krok(n+1)=1.515450e-004; ng(n+1)=1.576888e+003;
n=9; farx(n+1)=3.897082e+002; foe(n+1)=7.305139e+002; krok(n+1)=9.458811e-005; ng(n+1)=2.087342e+003;
n=10; farx(n+1)=4.038854e+002; foe(n+1)=7.127807e+002; krok(n+1)=1.445012e-004; ng(n+1)=1.612377e+003;
n=11; farx(n+1)=3.535435e+002; foe(n+1)=6.938286e+002; krok(n+1)=8.836179e-005; ng(n+1)=2.149728e+003;
n=12; farx(n+1)=3.651773e+002; foe(n+1)=6.759765e+002; krok(n+1)=1.383684e-004; ng(n+1)=1.655843e+003;
n=13; farx(n+1)=3.191746e+002; foe(n+1)=6.568851e+002; krok(n+1)=8.233132e-005; ng(n+1)=2.244732e+003;
n=14; farx(n+1)=3.288520e+002; foe(n+1)=6.390702e+002; krok(n+1)=1.297027e-004; ng(n+1)=1.709447e+003;
n=15; farx(n+1)=2.875044e+002; foe(n+1)=6.201009e+002; krok(n+1)=7.651456e-005; ng(n+1)=2.326846e+003;
n=16; farx(n+1)=2.954278e+002; foe(n+1)=6.023696e+002; krok(n+1)=1.220198e-004; ng(n+1)=1.749400e+003;
n=17; farx(n+1)=2.584158e+002; foe(n+1)=5.836873e+002; krok(n+1)=7.103937e-005; ng(n+1)=2.409198e+003;
n=18; farx(n+1)=2.650242e+002; foe(n+1)=5.663018e+002; krok(n+1)=1.133097e-004; ng(n+1)=1.791407e+003;
n=19; farx(n+1)=2.320787e+002; foe(n+1)=5.481564e+002; krok(n+1)=6.615676e-005; ng(n+1)=2.478764e+003;
n=20; farx(n+1)=2.377018e+002; foe(n+1)=5.313797e+002; krok(n+1)=1.032442e-004; ng(n+1)=1.842773e+003;
n=21; farx(n+1)=2.087567e+002; foe(n+1)=5.140605e+002; krok(n+1)=6.158863e-005; ng(n+1)=2.527368e+003;
n=22; farx(n+1)=2.135295e+002; foe(n+1)=4.978954e+002; krok(n+1)=9.579321e-005; ng(n+1)=1.879364e+003;
n=23; farx(n+1)=1.878726e+002; foe(n+1)=4.813322e+002; krok(n+1)=5.714551e-005; ng(n+1)=2.596092e+003;
n=24; farx(n+1)=1.921094e+002; foe(n+1)=4.659463e+002; krok(n+1)=8.878196e-005; ng(n+1)=1.923321e+003;
n=25; farx(n+1)=1.693883e+002; foe(n+1)=4.501410e+002; krok(n+1)=5.261508e-005; ng(n+1)=2.672378e+003;
n=26; farx(n+1)=1.729905e+002; foe(n+1)=4.354719e+002; krok(n+1)=8.271938e-005; ng(n+1)=1.955001e+003;
n=27; farx(n+1)=1.528170e+002; foe(n+1)=4.205278e+002; krok(n+1)=4.888219e-005; ng(n+1)=2.738029e+003;
n=28; farx(n+1)=1.560684e+002; foe(n+1)=4.067768e+002; krok(n+1)=7.571652e-005; ng(n+1)=2.011575e+003;
n=29; farx(n+1)=1.382163e+002; foe(n+1)=3.927769e+002; krok(n+1)=4.554340e-005; ng(n+1)=2.783063e+003;
n=30; farx(n+1)=1.411264e+002; foe(n+1)=3.799150e+002; krok(n+1)=6.946355e-005; ng(n+1)=2.057838e+003;
n=31; farx(n+1)=1.253546e+002; foe(n+1)=3.668748e+002; krok(n+1)=4.244958e-005; ng(n+1)=2.813560e+003;
n=32; farx(n+1)=1.279369e+002; foe(n+1)=3.548442e+002; krok(n+1)=6.447294e-005; ng(n+1)=2.086700e+003;
n=33; farx(n+1)=1.139396e+002; foe(n+1)=3.427012e+002; krok(n+1)=3.953838e-005; ng(n+1)=2.843110e+003;
n=34; farx(n+1)=1.162265e+002; foe(n+1)=3.314864e+002; krok(n+1)=5.990465e-005; ng(n+1)=2.105858e+003;
n=35; farx(n+1)=1.038041e+002; foe(n+1)=3.202401e+002; krok(n+1)=3.693259e-005; ng(n+1)=2.861368e+003;
n=36; farx(n+1)=1.058497e+002; foe(n+1)=3.098160e+002; krok(n+1)=5.618446e-005; ng(n+1)=2.114469e+003;
n=37; farx(n+1)=9.479125e+001; foe(n+1)=2.993908e+002; krok(n+1)=3.437390e-005; ng(n+1)=2.884234e+003;
n=38; farx(n+1)=9.658496e+001; foe(n+1)=2.896960e+002; krok(n+1)=5.314479e-005; ng(n+1)=2.104441e+003;
n=39; farx(n+1)=8.672608e+001; foe(n+1)=2.800614e+002; krok(n+1)=3.204087e-005; ng(n+1)=2.897671e+003;
n=40; farx(n+1)=8.830698e+001; foe(n+1)=2.710705e+002; krok(n+1)=5.059212e-005; ng(n+1)=2.085789e+003;
n=41; farx(n+1)=7.947585e+001; foe(n+1)=2.621733e+002; krok(n+1)=2.993183e-005; ng(n+1)=2.906928e+003;
n=42; farx(n+1)=8.090686e+001; foe(n+1)=2.538835e+002; krok(n+1)=4.837788e-005; ng(n+1)=2.069546e+003;
n=43; farx(n+1)=7.299478e+001; foe(n+1)=2.456556e+002; krok(n+1)=2.780862e-005; ng(n+1)=2.916251e+003;
n=44; farx(n+1)=7.422416e+001; foe(n+1)=2.379705e+002; krok(n+1)=4.661880e-005; ng(n+1)=2.033137e+003;
n=45; farx(n+1)=6.718073e+001; foe(n+1)=2.304097e+002; krok(n+1)=2.581106e-005; ng(n+1)=2.907974e+003;
n=46; farx(n+1)=6.820407e+001; foe(n+1)=2.232313e+002; krok(n+1)=4.614756e-005; ng(n+1)=1.977278e+003;
n=47; farx(n+1)=6.184141e+001; foe(n+1)=2.162076e+002; krok(n+1)=2.394830e-005; ng(n+1)=2.923743e+003;
n=48; farx(n+1)=6.271447e+001; foe(n+1)=2.095809e+002; krok(n+1)=4.461789e-005; ng(n+1)=1.933019e+003;
n=49; farx(n+1)=5.698221e+001; foe(n+1)=2.031639e+002; krok(n+1)=2.262551e-005; ng(n+1)=2.900768e+003;
n=50; farx(n+1)=5.781181e+001; foe(n+1)=1.971651e+002; krok(n+1)=4.244958e-005; ng(n+1)=1.907100e+003;
n=51; farx(n+1)=5.266243e+001; foe(n+1)=1.913090e+002; krok(n+1)=2.132008e-005; ng(n+1)=2.870065e+003;
n=52; farx(n+1)=5.339876e+001; foe(n+1)=1.858234e+002; krok(n+1)=4.068304e-005; ng(n+1)=1.864263e+003;
n=53; farx(n+1)=4.877986e+001; foe(n+1)=1.805135e+002; krok(n+1)=2.020969e-005; ng(n+1)=2.825229e+003;
n=54; farx(n+1)=4.944837e+001; foe(n+1)=1.755115e+002; krok(n+1)=3.909086e-005; ng(n+1)=1.826077e+003;
n=55; farx(n+1)=4.528337e+001; foe(n+1)=1.706835e+002; krok(n+1)=1.920477e-005; ng(n+1)=2.782719e+003;
n=56; farx(n+1)=4.590370e+001; foe(n+1)=1.661396e+002; krok(n+1)=3.746225e-005; ng(n+1)=1.788595e+003;
n=57; farx(n+1)=4.214156e+001; foe(n+1)=1.617535e+002; krok(n+1)=1.831129e-005; ng(n+1)=2.733658e+003;
n=58; farx(n+1)=4.272292e+001; foe(n+1)=1.576303e+002; krok(n+1)=3.592136e-005; ng(n+1)=1.750011e+003;
n=59; farx(n+1)=3.932181e+001; foe(n+1)=1.536438e+002; krok(n+1)=1.747693e-005; ng(n+1)=2.681030e+003;
n=60; farx(n+1)=3.986072e+001; foe(n+1)=1.498895e+002; krok(n+1)=3.463091e-005; ng(n+1)=1.707655e+003;
n=61; farx(n+1)=3.678555e+001; foe(n+1)=1.462608e+002; krok(n+1)=1.666649e-005; ng(n+1)=2.627233e+003;
n=62; farx(n+1)=3.727380e+001; foe(n+1)=1.428203e+002; krok(n+1)=3.371798e-005; ng(n+1)=1.661057e+003;
n=63; farx(n+1)=3.447716e+001; foe(n+1)=1.395067e+002; krok(n+1)=1.595012e-005; ng(n+1)=2.575819e+003;
n=64; farx(n+1)=3.493289e+001; foe(n+1)=1.363719e+002; krok(n+1)=3.261438e-005; ng(n+1)=1.618171e+003;
n=65; farx(n+1)=3.238320e+001; foe(n+1)=1.333508e+002; krok(n+1)=1.535151e-005; ng(n+1)=2.518245e+003;
n=66; farx(n+1)=3.281693e+001; foe(n+1)=1.305087e+002; krok(n+1)=3.121480e-005; ng(n+1)=1.578681e+003;
n=67; farx(n+1)=3.050198e+001; foe(n+1)=1.277700e+002; krok(n+1)=1.482745e-005; ng(n+1)=2.448446e+003;
n=68; farx(n+1)=3.091118e+001; foe(n+1)=1.251782e+002; krok(n+1)=3.023762e-005; ng(n+1)=1.537061e+003;
n=69; farx(n+1)=2.881475e+001; foe(n+1)=1.226755e+002; krok(n+1)=1.416371e-005; ng(n+1)=2.390388e+003;
n=70; farx(n+1)=2.917196e+001; foe(n+1)=1.202509e+002; krok(n+1)=3.059974e-005; ng(n+1)=1.485071e+003;
n=71; farx(n+1)=2.722622e+001; foe(n+1)=1.179133e+002; krok(n+1)=1.357488e-005; ng(n+1)=2.366648e+003;
n=72; farx(n+1)=2.756410e+001; foe(n+1)=1.157043e+002; krok(n+1)=2.952483e-005; ng(n+1)=1.443630e+003;
n=73; farx(n+1)=2.579112e+001; foe(n+1)=1.135785e+002; krok(n+1)=1.315377e-005; ng(n+1)=2.296851e+003;
n=74; farx(n+1)=2.610635e+001; foe(n+1)=1.115523e+002; krok(n+1)=2.877521e-005; ng(n+1)=1.400961e+003;
n=75; farx(n+1)=2.447622e+001; foe(n+1)=1.096074e+002; krok(n+1)=1.280413e-005; ng(n+1)=2.236848e+003;
n=76; farx(n+1)=2.477963e+001; foe(n+1)=1.077673e+002; krok(n+1)=2.746916e-005; ng(n+1)=1.362858e+003;
n=77; farx(n+1)=2.328615e+001; foe(n+1)=1.060076e+002; krok(n+1)=1.259177e-005; ng(n+1)=2.156006e+003;
n=78; farx(n+1)=2.358781e+001; foe(n+1)=1.043429e+002; krok(n+1)=2.624535e-005; ng(n+1)=1.328294e+003;
n=79; farx(n+1)=2.221887e+001; foe(n+1)=1.027368e+002; krok(n+1)=1.229246e-005; ng(n+1)=2.087944e+003;
n=80; farx(n+1)=2.250441e+001; foe(n+1)=1.012014e+002; krok(n+1)=2.580710e-005; ng(n+1)=1.288658e+003;
n=81; farx(n+1)=2.123475e+001; foe(n+1)=9.971446e+001; krok(n+1)=1.197513e-005; ng(n+1)=2.042481e+003;
n=82; farx(n+1)=2.150744e+001; foe(n+1)=9.830190e+001; krok(n+1)=2.518354e-005; ng(n+1)=1.251286e+003;
n=83; farx(n+1)=2.033181e+001; foe(n+1)=9.693079e+001; krok(n+1)=1.171380e-005; ng(n+1)=1.988516e+003;
n=84; farx(n+1)=2.059094e+001; foe(n+1)=9.562981e+001; krok(n+1)=2.443486e-005; ng(n+1)=1.215281e+003;
n=85; farx(n+1)=1.950587e+001; foe(n+1)=9.437173e+001; krok(n+1)=1.149974e-005; ng(n+1)=1.925120e+003;
n=86; farx(n+1)=1.975651e+001; foe(n+1)=9.316914e+001; krok(n+1)=2.418894e-005; ng(n+1)=1.181069e+003;
n=87; farx(n+1)=1.873646e+001; foe(n+1)=9.199111e+001; krok(n+1)=1.126542e-005; ng(n+1)=1.892205e+003;
n=88; farx(n+1)=1.897738e+001; foe(n+1)=9.088902e+001; krok(n+1)=2.299949e-005; ng(n+1)=1.150351e+003;
n=89; farx(n+1)=1.803409e+001; foe(n+1)=8.982056e+001; krok(n+1)=1.126542e-005; ng(n+1)=1.809015e+003;
n=90; farx(n+1)=1.827514e+001; foe(n+1)=8.881628e+001; krok(n+1)=2.190752e-005; ng(n+1)=1.124097e+003;
n=91; farx(n+1)=1.740809e+001; foe(n+1)=8.783396e+001; krok(n+1)=1.104522e-005; ng(n+1)=1.744594e+003;
n=92; farx(n+1)=1.762967e+001; foe(n+1)=8.688337e+001; krok(n+1)=2.212034e-005; ng(n+1)=1.088195e+003;
n=93; farx(n+1)=1.682078e+001; foe(n+1)=8.595863e+001; krok(n+1)=1.069939e-005; ng(n+1)=1.715583e+003;
n=94; farx(n+1)=1.702458e+001; foe(n+1)=8.505252e+001; krok(n+1)=2.277170e-005; ng(n+1)=1.052379e+003;
n=95; farx(n+1)=1.625190e+001; foe(n+1)=8.416866e+001; krok(n+1)=1.046074e-005; ng(n+1)=1.706396e+003;
n=96; farx(n+1)=1.644771e+001; foe(n+1)=8.333299e+001; krok(n+1)=2.170937e-005; ng(n+1)=1.024960e+003;
n=97; farx(n+1)=1.573397e+001; foe(n+1)=8.252776e+001; krok(n+1)=1.044207e-005; ng(n+1)=1.629955e+003;
n=98; farx(n+1)=1.592586e+001; foe(n+1)=8.175163e+001; krok(n+1)=2.132008e-005; ng(n+1)=9.980489e+002;
n=99; farx(n+1)=1.525753e+001; foe(n+1)=8.099656e+001; krok(n+1)=1.023989e-005; ng(n+1)=1.592175e+003;
n=100; farx(n+1)=1.544131e+001; foe(n+1)=8.026297e+001; krok(n+1)=2.167257e-005; ng(n+1)=9.683895e+002;
n=101; farx(n+1)=1.480689e+001; foe(n+1)=7.954003e+001; krok(n+1)=9.955536e-006; ng(n+1)=1.580556e+003;
n=102; farx(n+1)=1.497529e+001; foe(n+1)=7.884443e+001; krok(n+1)=2.139878e-005; ng(n+1)=9.389180e+002;
n=103; farx(n+1)=1.438317e+001; foe(n+1)=7.817519e+001; krok(n+1)=9.888364e-006; ng(n+1)=1.525750e+003;
n=104; farx(n+1)=1.454697e+001; foe(n+1)=7.752432e+001; krok(n+1)=2.114892e-005; ng(n+1)=9.135914e+002;
n=105; farx(n+1)=1.398542e+001; foe(n+1)=7.689419e+001; krok(n+1)=9.809075e-006; ng(n+1)=1.491489e+003;
n=106; farx(n+1)=1.414667e+001; foe(n+1)=7.629033e+001; krok(n+1)=2.047978e-005; ng(n+1)=8.911879e+002;
n=107; farx(n+1)=1.362155e+001; foe(n+1)=7.570346e+001; krok(n+1)=9.716744e-006; ng(n+1)=1.443557e+003;
n=108; farx(n+1)=1.377690e+001; foe(n+1)=7.513201e+001; krok(n+1)=2.058283e-005; ng(n+1)=8.662983e+002;
n=109; farx(n+1)=1.327738e+001; foe(n+1)=7.457295e+001; krok(n+1)=9.554760e-006; ng(n+1)=1.422051e+003;
n=110; farx(n+1)=1.342652e+001; foe(n+1)=7.403324e+001; krok(n+1)=2.044242e-005; ng(n+1)=8.426832e+002;
n=111; farx(n+1)=1.295469e+001; foe(n+1)=7.350569e+001; krok(n+1)=9.434846e-006; ng(n+1)=1.389186e+003;
n=112; farx(n+1)=1.309744e+001; foe(n+1)=7.299431e+001; krok(n+1)=2.036880e-005; ng(n+1)=8.195296e+002;
n=113; farx(n+1)=1.264949e+001; foe(n+1)=7.249551e+001; krok(n+1)=9.351708e-006; ng(n+1)=1.358344e+003;
n=114; farx(n+1)=1.278949e+001; foe(n+1)=7.201469e+001; krok(n+1)=2.018691e-005; ng(n+1)=7.987394e+002;
n=115; farx(n+1)=1.236466e+001; foe(n+1)=7.154134e+001; krok(n+1)=9.225986e-006; ng(n+1)=1.331059e+003;
n=116; farx(n+1)=1.249869e+001; foe(n+1)=7.108407e+001; krok(n+1)=2.018691e-005; ng(n+1)=7.769437e+002;
n=117; farx(n+1)=1.209529e+001; foe(n+1)=7.063493e+001; krok(n+1)=9.119360e-006; ng(n+1)=1.303496e+003;
n=118; farx(n+1)=1.222476e+001; foe(n+1)=7.020076e+001; krok(n+1)=2.018691e-005; ng(n+1)=7.565691e+002;
n=119; farx(n+1)=1.184107e+001; foe(n+1)=6.977330e+001; krok(n+1)=9.003600e-006; ng(n+1)=1.278824e+003;
n=120; farx(n+1)=1.196540e+001; foe(n+1)=6.935990e+001; krok(n+1)=2.018691e-005; ng(n+1)=7.365898e+002;
n=121; farx(n+1)=1.160125e+001; foe(n+1)=6.895371e+001; krok(n+1)=8.887598e-006; ng(n+1)=1.252434e+003;
n=122; farx(n+1)=1.172006e+001; foe(n+1)=6.855812e+001; krok(n+1)=2.031518e-005; ng(n+1)=7.167397e+002;
n=123; farx(n+1)=1.137221e+001; foe(n+1)=6.817087e+001; krok(n+1)=8.804989e-006; ng(n+1)=1.229181e+003;
n=124; farx(n+1)=1.148807e+001; foe(n+1)=6.779627e+001; krok(n+1)=2.018691e-005; ng(n+1)=6.988905e+002;
n=125; farx(n+1)=1.115626e+001; foe(n+1)=6.742808e+001; krok(n+1)=8.739066e-006; ng(n+1)=1.203275e+003;
n=126; farx(n+1)=1.126834e+001; foe(n+1)=6.707311e+001; krok(n+1)=1.991107e-005; ng(n+1)=6.817079e+002;
n=127; farx(n+1)=1.095290e+001; foe(n+1)=6.672574e+001; krok(n+1)=8.718949e-006; ng(n+1)=1.170490e+003;
n=128; farx(n+1)=1.106209e+001; foe(n+1)=6.639038e+001; krok(n+1)=1.954543e-005; ng(n+1)=6.657169e+002;
n=129; farx(n+1)=1.076181e+001; foe(n+1)=6.606275e+001; krok(n+1)=8.718949e-006; ng(n+1)=1.137578e+003;
n=130; farx(n+1)=1.086902e+001; foe(n+1)=6.574648e+001; krok(n+1)=1.916970e-005; ng(n+1)=6.507819e+002;
n=131; farx(n+1)=1.058341e+001; foe(n+1)=6.543656e+001; krok(n+1)=8.686387e-006; ng(n+1)=1.107469e+003;
n=132; farx(n+1)=1.068773e+001; foe(n+1)=6.513531e+001; krok(n+1)=1.908189e-005; ng(n+1)=6.354163e+002;
n=133; farx(n+1)=1.041393e+001; foe(n+1)=6.483972e+001; krok(n+1)=8.648025e-006; ng(n+1)=1.084631e+003;
n=134; farx(n+1)=1.051570e+001; foe(n+1)=6.455393e+001; krok(n+1)=1.881044e-005; ng(n+1)=6.211364e+002;
n=135; farx(n+1)=1.025373e+001; foe(n+1)=6.427354e+001; krok(n+1)=8.648025e-006; ng(n+1)=1.057676e+003;
n=136; farx(n+1)=1.035342e+001; foe(n+1)=6.400313e+001; krok(n+1)=1.843594e-005; ng(n+1)=6.078598e+002;
n=137; farx(n+1)=1.010338e+001; foe(n+1)=6.373764e+001; krok(n+1)=8.648025e-006; ng(n+1)=1.029621e+003;
n=138; farx(n+1)=1.020074e+001; foe(n+1)=6.348053e+001; krok(n+1)=1.819728e-005; ng(n+1)=5.946883e+002;
n=139; farx(n+1)=9.960949e+000; foe(n+1)=6.322799e+001; krok(n+1)=8.648025e-006; ng(n+1)=1.005749e+003;
n=140; farx(n+1)=1.005622e+001; foe(n+1)=6.298398e+001; krok(n+1)=1.787154e-005; ng(n+1)=5.823368e+002;
n=141; farx(n+1)=9.826851e+000; foe(n+1)=6.274422e+001; krok(n+1)=8.648025e-006; ng(n+1)=9.802190e+002;
n=142; farx(n+1)=9.919775e+000; foe(n+1)=6.251160e+001; krok(n+1)=1.765974e-005; ng(n+1)=5.700785e+002;
n=143; farx(n+1)=9.699550e+000; foe(n+1)=6.228316e+001; krok(n+1)=8.648025e-006; ng(n+1)=9.578713e+002;
n=144; farx(n+1)=9.790628e+000; foe(n+1)=6.206166e+001; krok(n+1)=1.743790e-005; ng(n+1)=5.586097e+002;
n=145; farx(n+1)=9.579660e+000; foe(n+1)=6.184361e+001; krok(n+1)=8.610243e-006; ng(n+1)=9.371323e+002;
n=146; farx(n+1)=9.667994e+000; foe(n+1)=6.163072e+001; krok(n+1)=1.744002e-005; ng(n+1)=5.463064e+002;
n=147; farx(n+1)=9.465380e+000; foe(n+1)=6.142163e+001; krok(n+1)=8.546025e-006; ng(n+1)=9.204896e+002;
n=148; farx(n+1)=9.551316e+000; foe(n+1)=6.121633e+001; krok(n+1)=1.760998e-005; ng(n+1)=5.340301e+002;
n=149; farx(n+1)=9.355304e+000; foe(n+1)=6.101436e+001; krok(n+1)=8.470861e-006; ng(n+1)=9.098753e+002;
n=150; farx(n+1)=9.438872e+000; foe(n+1)=6.081705e+001; krok(n+1)=1.761549e-005; ng(n+1)=5.226325e+002;
n=151; farx(n+1)=9.249570e+000; foe(n+1)=6.062313e+001; krok(n+1)=8.445281e-006; ng(n+1)=8.943129e+002;
n=152; farx(n+1)=9.331450e+000; foe(n+1)=6.043456e+001; krok(n+1)=1.747813e-005; ng(n+1)=5.124618e+002;
n=153; farx(n+1)=9.149917e+000; foe(n+1)=6.024876e+001; krok(n+1)=8.379583e-006; ng(n+1)=8.771644e+002;
n=154; farx(n+1)=9.229225e+000; foe(n+1)=6.006598e+001; krok(n+1)=1.773106e-005; ng(n+1)=5.011748e+002;
n=155; farx(n+1)=9.053194e+000; foe(n+1)=5.988620e+001; krok(n+1)=8.311102e-006; ng(n+1)=8.678381e+002;
n=156; farx(n+1)=9.130628e+000; foe(n+1)=5.971047e+001; krok(n+1)=1.775984e-005; ng(n+1)=4.909722e+002;
n=157; farx(n+1)=8.960361e+000; foe(n+1)=5.953738e+001; krok(n+1)=8.272636e-006; ng(n+1)=8.546226e+002;
n=158; farx(n+1)=9.035938e+000; foe(n+1)=5.936882e+001; krok(n+1)=1.767686e-005; ng(n+1)=4.812783e+002;
n=159; farx(n+1)=8.871339e+000; foe(n+1)=5.920283e+001; krok(n+1)=8.269595e-006; ng(n+1)=8.385161e+002;
n=160; farx(n+1)=8.945365e+000; foe(n+1)=5.904189e+001; krok(n+1)=1.743790e-005; ng(n+1)=4.723644e+002;
n=161; farx(n+1)=8.787356e+000; foe(n+1)=5.888333e+001; krok(n+1)=8.245594e-006; ng(n+1)=8.196825e+002;
n=162; farx(n+1)=8.859561e+000; foe(n+1)=5.872796e+001; krok(n+1)=1.755604e-005; ng(n+1)=4.627357e+002;
n=163; farx(n+1)=8.706155e+000; foe(n+1)=5.857476e+001; krok(n+1)=8.205260e-006; ng(n+1)=8.091402e+002;
n=164; farx(n+1)=8.776369e+000; foe(n+1)=5.842563e+001; krok(n+1)=1.738167e-005; ng(n+1)=4.538926e+002;
n=165; farx(n+1)=8.628713e+000; foe(n+1)=5.827930e+001; krok(n+1)=8.206521e-006; ng(n+1)=7.905789e+002;
n=166; farx(n+1)=8.697484e+000; foe(n+1)=5.813603e+001; krok(n+1)=1.733605e-005; ng(n+1)=4.452740e+002;
n=167; farx(n+1)=8.554503e+000; foe(n+1)=5.799517e+001; krok(n+1)=8.182517e-006; ng(n+1)=7.775204e+002;
n=168; farx(n+1)=8.621734e+000; foe(n+1)=5.785747e+001; krok(n+1)=1.726474e-005; ng(n+1)=4.368707e+002;
n=169; farx(n+1)=8.483387e+000; foe(n+1)=5.772207e+001; krok(n+1)=8.166804e-006; ng(n+1)=7.635672e+002;
n=170; farx(n+1)=8.549135e+000; foe(n+1)=5.758969e+001; krok(n+1)=1.718695e-005; ng(n+1)=4.287133e+002;
n=171; farx(n+1)=8.415052e+000; foe(n+1)=5.745952e+001; krok(n+1)=8.164371e-006; ng(n+1)=7.497954e+002;
n=172; farx(n+1)=8.479726e+000; foe(n+1)=5.733261e+001; krok(n+1)=1.707914e-005; ng(n+1)=4.210506e+002;
n=173; farx(n+1)=8.349700e+000; foe(n+1)=5.720726e+001; krok(n+1)=8.151193e-006; ng(n+1)=7.371590e+002;
n=174; farx(n+1)=8.412677e+000; foe(n+1)=5.708529e+001; krok(n+1)=1.689172e-005; ng(n+1)=4.134367e+002;
n=175; farx(n+1)=8.287353e+000; foe(n+1)=5.696539e+001; krok(n+1)=8.155065e-006; ng(n+1)=7.202606e+002;
n=176; farx(n+1)=8.349198e+000; foe(n+1)=5.684780e+001; krok(n+1)=1.693219e-005; ng(n+1)=4.059732e+002;
n=177; farx(n+1)=8.227356e+000; foe(n+1)=5.673175e+001; krok(n+1)=8.112714e-006; ng(n+1)=7.110447e+002;
n=178; farx(n+1)=8.287748e+000; foe(n+1)=5.661819e+001; krok(n+1)=1.693774e-005; ng(n+1)=3.985321e+002;
n=179; farx(n+1)=8.169946e+000; foe(n+1)=5.650626e+001; krok(n+1)=8.058704e-006; ng(n+1)=6.996879e+002;
n=180; farx(n+1)=8.228478e+000; foe(n+1)=5.639576e+001; krok(n+1)=1.709205e-005; ng(n+1)=3.905526e+002;
n=181; farx(n+1)=8.113567e+000; foe(n+1)=5.628757e+001; krok(n+1)=8.059117e-006; ng(n+1)=6.897149e+002;
n=182; farx(n+1)=8.171294e+000; foe(n+1)=5.618183e+001; krok(n+1)=1.689172e-005; ng(n+1)=3.843814e+002;
n=183; farx(n+1)=8.059238e+000; foe(n+1)=5.607776e+001; krok(n+1)=8.106422e-006; ng(n+1)=6.768670e+002;
n=184; farx(n+1)=8.115947e+000; foe(n+1)=5.597724e+001; krok(n+1)=1.637351e-005; ng(n+1)=3.791413e+002;
n=185; farx(n+1)=8.008206e+000; foe(n+1)=5.587832e+001; krok(n+1)=8.164371e-006; ng(n+1)=6.566918e+002;
n=186; farx(n+1)=8.063958e+000; foe(n+1)=5.578164e+001; krok(n+1)=1.625776e-005; ng(n+1)=3.730671e+002;
n=187; farx(n+1)=7.959344e+000; foe(n+1)=5.568623e+001; krok(n+1)=8.153596e-006; ng(n+1)=6.454730e+002;
n=188; farx(n+1)=8.013539e+000; foe(n+1)=5.559283e+001; krok(n+1)=1.611741e-005; ng(n+1)=3.667465e+002;
n=189; farx(n+1)=7.911920e+000; foe(n+1)=5.550130e+001; krok(n+1)=8.205260e-006; ng(n+1)=6.312062e+002;
n=190; farx(n+1)=7.965658e+000; foe(n+1)=5.541186e+001; krok(n+1)=1.592832e-005; ng(n+1)=3.615594e+002;
n=191; farx(n+1)=7.866958e+000; foe(n+1)=5.532343e+001; krok(n+1)=8.182517e-006; ng(n+1)=6.210508e+002;
n=192; farx(n+1)=7.919526e+000; foe(n+1)=5.523674e+001; krok(n+1)=1.596631e-005; ng(n+1)=3.553097e+002;
n=193; farx(n+1)=7.823313e+000; foe(n+1)=5.515110e+001; krok(n+1)=8.164371e-006; ng(n+1)=6.123794e+002;
n=194; farx(n+1)=7.874699e+000; foe(n+1)=5.506731e+001; krok(n+1)=1.588017e-005; ng(n+1)=3.495058e+002;
n=195; farx(n+1)=7.781229e+000; foe(n+1)=5.498474e+001; krok(n+1)=8.166804e-006; ng(n+1)=6.012291e+002;
n=196; farx(n+1)=7.831773e+000; foe(n+1)=5.490373e+001; krok(n+1)=1.587668e-005; ng(n+1)=3.440929e+002;
n=197; farx(n+1)=7.740737e+000; foe(n+1)=5.482366e+001; krok(n+1)=8.129247e-006; ng(n+1)=5.934289e+002;
n=198; farx(n+1)=7.790241e+000; foe(n+1)=5.474498e+001; krok(n+1)=1.596631e-005; ng(n+1)=3.384942e+002;
n=199; farx(n+1)=7.701750e+000; foe(n+1)=5.466722e+001; krok(n+1)=8.058704e-006; ng(n+1)=5.866730e+002;
n=200; farx(n+1)=7.750231e+000; foe(n+1)=5.459023e+001; krok(n+1)=1.630239e-005; ng(n+1)=3.323826e+002;
n=201; farx(n+1)=7.663297e+000; foe(n+1)=5.451409e+001; krok(n+1)=7.989188e-006; ng(n+1)=5.844815e+002;
n=202; farx(n+1)=7.710694e+000; foe(n+1)=5.443934e+001; krok(n+1)=1.633595e-005; ng(n+1)=3.270547e+002;
n=203; farx(n+1)=7.625849e+000; foe(n+1)=5.436559e+001; krok(n+1)=7.983154e-006; ng(n+1)=5.762412e+002;
n=204; farx(n+1)=7.672405e+000; foe(n+1)=5.429338e+001; krok(n+1)=1.625776e-005; ng(n+1)=3.222805e+002;
n=205; farx(n+1)=7.589706e+000; foe(n+1)=5.422206e+001; krok(n+1)=7.983154e-006; ng(n+1)=5.672267e+002;
n=206; farx(n+1)=7.635504e+000; foe(n+1)=5.415224e+001; krok(n+1)=1.621284e-005; ng(n+1)=3.175935e+002;
n=207; farx(n+1)=7.555116e+000; foe(n+1)=5.408312e+001; krok(n+1)=7.940086e-006; ng(n+1)=5.593066e+002;
n=208; farx(n+1)=7.599897e+000; foe(n+1)=5.401493e+001; krok(n+1)=1.641052e-005; ng(n+1)=3.123237e+002;
n=209; farx(n+1)=7.521096e+000; foe(n+1)=5.394757e+001; krok(n+1)=7.909320e-006; ng(n+1)=5.542785e+002;
n=210; farx(n+1)=7.564961e+000; foe(n+1)=5.388153e+001; krok(n+1)=1.634484e-005; ng(n+1)=3.077226e+002;
n=211; farx(n+1)=7.488120e+000; foe(n+1)=5.381640e+001; krok(n+1)=7.918258e-006; ng(n+1)=5.452502e+002;
n=212; farx(n+1)=7.531276e+000; foe(n+1)=5.375257e+001; krok(n+1)=1.625776e-005; ng(n+1)=3.034774e+002;
n=213; farx(n+1)=7.456437e+000; foe(n+1)=5.368950e+001; krok(n+1)=7.902593e-006; ng(n+1)=5.371549e+002;
n=214; farx(n+1)=7.498815e+000; foe(n+1)=5.362743e+001; krok(n+1)=1.632874e-005; ng(n+1)=2.988963e+002;
n=215; farx(n+1)=7.425693e+000; foe(n+1)=5.356608e+001; krok(n+1)=7.866658e-006; ng(n+1)=5.314517e+002;
n=216; farx(n+1)=7.467296e+000; foe(n+1)=5.350563e+001; krok(n+1)=1.643953e-005; ng(n+1)=2.944619e+002;
n=217; farx(n+1)=7.395703e+000; foe(n+1)=5.344586e+001; krok(n+1)=7.835092e-006; ng(n+1)=5.263971e+002;
n=218; farx(n+1)=7.436482e+000; foe(n+1)=5.338711e+001; krok(n+1)=1.645337e-005; ng(n+1)=2.902793e+002;
n=219; farx(n+1)=7.366607e+000; foe(n+1)=5.332909e+001; krok(n+1)=7.817314e-006; ng(n+1)=5.194580e+002;
n=220; farx(n+1)=7.406749e+000; foe(n+1)=5.327195e+001; krok(n+1)=1.653919e-005; ng(n+1)=2.862435e+002;
n=221; farx(n+1)=7.338218e+000; foe(n+1)=5.321539e+001; krok(n+1)=7.794783e-006; ng(n+1)=5.146973e+002;
n=222; farx(n+1)=7.377625e+000; foe(n+1)=5.315996e+001; krok(n+1)=1.649003e-005; ng(n+1)=2.825266e+002;
n=223; farx(n+1)=7.310913e+000; foe(n+1)=5.310513e+001; krok(n+1)=7.768137e-006; ng(n+1)=5.072790e+002;
n=224; farx(n+1)=7.349589e+000; foe(n+1)=5.305094e+001; krok(n+1)=1.666809e-005; ng(n+1)=2.784983e+002;
n=225; farx(n+1)=7.284400e+000; foe(n+1)=5.299738e+001; krok(n+1)=7.698578e-006; ng(n+1)=5.033759e+002;
n=226; farx(n+1)=7.322466e+000; foe(n+1)=5.294409e+001; krok(n+1)=1.707914e-005; ng(n+1)=2.743097e+002;
n=227; farx(n+1)=7.258119e+000; foe(n+1)=5.289127e+001; krok(n+1)=7.626235e-006; ng(n+1)=5.034545e+002;
n=228; farx(n+1)=7.295195e+000; foe(n+1)=5.283924e+001; krok(n+1)=1.709205e-005; ng(n+1)=2.706112e+002;
n=229; farx(n+1)=7.232385e+000; foe(n+1)=5.278812e+001; krok(n+1)=7.638207e-006; ng(n+1)=4.953760e+002;
n=230; farx(n+1)=7.268883e+000; foe(n+1)=5.273773e+001; krok(n+1)=1.700060e-005; ng(n+1)=2.674159e+002;
n=231; farx(n+1)=7.207505e+000; foe(n+1)=5.268816e+001; krok(n+1)=7.640629e-006; ng(n+1)=4.882480e+002;
n=232; farx(n+1)=7.243717e+000; foe(n+1)=5.263926e+001; krok(n+1)=1.707914e-005; ng(n+1)=2.642177e+002;
n=233; farx(n+1)=7.183496e+000; foe(n+1)=5.259070e+001; krok(n+1)=7.585237e-006; ng(n+1)=4.852809e+002;
n=234; farx(n+1)=7.218894e+000; foe(n+1)=5.254279e+001; krok(n+1)=1.722049e-005; ng(n+1)=2.606076e+002;
n=235; farx(n+1)=7.159962e+000; foe(n+1)=5.249547e+001; krok(n+1)=7.559984e-006; ng(n+1)=4.800151e+002;
n=236; farx(n+1)=7.194794e+000; foe(n+1)=5.244871e+001; krok(n+1)=1.733471e-005; ng(n+1)=2.572581e+002;
n=237; farx(n+1)=7.136976e+000; foe(n+1)=5.240249e+001; krok(n+1)=7.532364e-006; ng(n+1)=4.758053e+002;
n=238; farx(n+1)=7.171322e+000; foe(n+1)=5.235690e+001; krok(n+1)=1.743790e-005; ng(n+1)=2.540287e+002;
n=239; farx(n+1)=7.114670e+000; foe(n+1)=5.231171e+001; krok(n+1)=7.488082e-006; ng(n+1)=4.719143e+002;
n=240; farx(n+1)=7.148387e+000; foe(n+1)=5.226700e+001; krok(n+1)=1.763687e-005; ng(n+1)=2.506116e+002;
n=241; farx(n+1)=7.092808e+000; foe(n+1)=5.222278e+001; krok(n+1)=7.447416e-006; ng(n+1)=4.683212e+002;
n=242; farx(n+1)=7.125887e+000; foe(n+1)=5.217900e+001; krok(n+1)=1.777520e-005; ng(n+1)=2.473169e+002;
n=243; farx(n+1)=7.071369e+000; foe(n+1)=5.213583e+001; krok(n+1)=7.426417e-006; ng(n+1)=4.638256e+002;
n=244; farx(n+1)=7.103936e+000; foe(n+1)=5.209314e+001; krok(n+1)=1.782892e-005; ng(n+1)=2.442672e+002;
n=245; farx(n+1)=7.050672e+000; foe(n+1)=5.205102e+001; krok(n+1)=7.381362e-006; ng(n+1)=4.590830e+002;
n=246; farx(n+1)=7.082775e+000; foe(n+1)=5.200895e+001; krok(n+1)=1.821997e-005; ng(n+1)=2.408594e+002;
n=247; farx(n+1)=7.030088e+000; foe(n+1)=5.196737e+001; krok(n+1)=7.336129e-006; ng(n+1)=4.586527e+002;
n=248; farx(n+1)=7.061697e+000; foe(n+1)=5.192637e+001; krok(n+1)=1.825798e-005; ng(n+1)=2.379128e+002;
n=249; farx(n+1)=7.009962e+000; foe(n+1)=5.188581e+001; krok(n+1)=7.334139e-006; ng(n+1)=4.537838e+002;
n=250; farx(n+1)=7.041088e+000; foe(n+1)=5.184603e+001; krok(n+1)=1.814366e-005; ng(n+1)=2.352101e+002;
n=251; farx(n+1)=6.990582e+000; foe(n+1)=5.180665e+001; krok(n+1)=7.325937e-006; ng(n+1)=4.472760e+002;
n=252; farx(n+1)=7.021202e+000; foe(n+1)=5.176775e+001; krok(n+1)=1.823411e-005; ng(n+1)=2.322910e+002;
n=253; farx(n+1)=6.971668e+000; foe(n+1)=5.172928e+001; krok(n+1)=7.303802e-006; ng(n+1)=4.429764e+002;
n=254; farx(n+1)=7.001836e+000; foe(n+1)=5.169123e+001; krok(n+1)=1.834698e-005; ng(n+1)=2.294245e+002;
n=255; farx(n+1)=6.953275e+000; foe(n+1)=5.165358e+001; krok(n+1)=7.265431e-006; ng(n+1)=4.392701e+002;
n=256; farx(n+1)=6.982994e+000; foe(n+1)=5.161614e+001; krok(n+1)=1.861995e-005; ng(n+1)=2.264121e+002;
n=257; farx(n+1)=6.934940e+000; foe(n+1)=5.157911e+001; krok(n+1)=7.258975e-006; ng(n+1)=4.372548e+002;
n=258; farx(n+1)=6.964207e+000; foe(n+1)=5.154294e+001; krok(n+1)=1.834698e-005; ng(n+1)=2.240639e+002;
n=259; farx(n+1)=6.917187e+000; foe(n+1)=5.150714e+001; krok(n+1)=7.304702e-006; ng(n+1)=4.294159e+002;
n=260; farx(n+1)=6.945959e+000; foe(n+1)=5.147232e+001; krok(n+1)=1.792006e-005; ng(n+1)=2.219405e+002;
n=261; farx(n+1)=6.900323e+000; foe(n+1)=5.143789e+001; krok(n+1)=7.336129e-006; ng(n+1)=4.197690e+002;
n=262; farx(n+1)=6.928661e+000; foe(n+1)=5.140393e+001; krok(n+1)=1.788662e-005; ng(n+1)=2.197159e+002;
n=263; farx(n+1)=6.883904e+000; foe(n+1)=5.137037e+001; krok(n+1)=7.336129e-006; ng(n+1)=4.147705e+002;
n=264; farx(n+1)=6.911879e+000; foe(n+1)=5.133724e+001; krok(n+1)=1.787641e-005; ng(n+1)=2.175778e+002;
n=265; farx(n+1)=6.867890e+000; foe(n+1)=5.130445e+001; krok(n+1)=7.334139e-006; ng(n+1)=4.105890e+002;
n=266; farx(n+1)=6.895500e+000; foe(n+1)=5.127216e+001; krok(n+1)=1.782680e-005; ng(n+1)=2.155222e+002;
n=267; farx(n+1)=6.852317e+000; foe(n+1)=5.124015e+001; krok(n+1)=7.334139e-006; ng(n+1)=4.060250e+002;
n=268; farx(n+1)=6.879474e+000; foe(n+1)=5.120866e+001; krok(n+1)=1.773994e-005; ng(n+1)=2.135894e+002;
n=269; farx(n+1)=6.837085e+000; foe(n+1)=5.117750e+001; krok(n+1)=7.359141e-006; ng(n+1)=4.005091e+002;
n=270; farx(n+1)=6.863938e+000; foe(n+1)=5.114696e+001; krok(n+1)=1.755604e-005; ng(n+1)=2.119554e+002;
n=271; farx(n+1)=6.822567e+000; foe(n+1)=5.111663e+001; krok(n+1)=7.338903e-006; ng(n+1)=3.950929e+002;
n=272; farx(n+1)=6.848960e+000; foe(n+1)=5.108647e+001; krok(n+1)=1.777520e-005; ng(n+1)=2.097698e+002;
n=273; farx(n+1)=6.808091e+000; foe(n+1)=5.105666e+001; krok(n+1)=7.334139e-006; ng(n+1)=3.926714e+002;
n=274; farx(n+1)=6.834146e+000; foe(n+1)=5.102732e+001; krok(n+1)=1.763713e-005; ng(n+1)=2.080966e+002;
n=275; farx(n+1)=6.793995e+000; foe(n+1)=5.099828e+001; krok(n+1)=7.359141e-006; ng(n+1)=3.875521e+002;
n=276; farx(n+1)=6.819653e+000; foe(n+1)=5.096979e+001; krok(n+1)=1.738167e-005; ng(n+1)=2.065457e+002;
n=277; farx(n+1)=6.780216e+000; foe(n+1)=5.094163e+001; krok(n+1)=7.425417e-006; ng(n+1)=3.809225e+002;
n=278; farx(n+1)=6.805534e+000; foe(n+1)=5.091421e+001; krok(n+1)=1.689172e-005; ng(n+1)=2.053600e+002;
n=279; farx(n+1)=6.767179e+000; foe(n+1)=5.088706e+001; krok(n+1)=7.470954e-006; ng(n+1)=3.724495e+002;
n=280; farx(n+1)=6.792076e+000; foe(n+1)=5.086027e+001; krok(n+1)=1.675917e-005; ng(n+1)=2.037262e+002;
n=281; farx(n+1)=6.754516e+000; foe(n+1)=5.083385e+001; krok(n+1)=7.488082e-006; ng(n+1)=3.669904e+002;
n=282; farx(n+1)=6.779291e+000; foe(n+1)=5.080763e+001; krok(n+1)=1.684728e-005; ng(n+1)=2.020920e+002;
n=283; farx(n+1)=6.742077e+000; foe(n+1)=5.078155e+001; krok(n+1)=7.471755e-006; ng(n+1)=3.658553e+002;
n=284; farx(n+1)=6.766409e+000; foe(n+1)=5.075594e+001; krok(n+1)=1.666649e-005; ng(n+1)=2.006766e+002;
n=285; farx(n+1)=6.730060e+000; foe(n+1)=5.073061e+001; krok(n+1)=7.488082e-006; ng(n+1)=3.598489e+002;
n=286; farx(n+1)=6.754166e+000; foe(n+1)=5.070546e+001; krok(n+1)=1.675717e-005; ng(n+1)=1.990298e+002;
n=287; farx(n+1)=6.718250e+000; foe(n+1)=5.068051e+001; krok(n+1)=7.470954e-006; ng(n+1)=3.580281e+002;
n=288; farx(n+1)=6.742080e+000; foe(n+1)=5.065586e+001; krok(n+1)=1.675717e-005; ng(n+1)=1.975097e+002;
n=289; farx(n+1)=6.706677e+000; foe(n+1)=5.063138e+001; krok(n+1)=7.471755e-006; ng(n+1)=3.549998e+002;
n=290; farx(n+1)=6.730150e+000; foe(n+1)=5.060727e+001; krok(n+1)=1.664416e-005; ng(n+1)=1.961136e+002;
n=291; farx(n+1)=6.695384e+000; foe(n+1)=5.058338e+001; krok(n+1)=7.493898e-006; ng(n+1)=3.503368e+002;
n=292; farx(n+1)=6.718621e+000; foe(n+1)=5.055980e+001; krok(n+1)=1.656922e-005; ng(n+1)=1.948202e+002;
n=293; farx(n+1)=6.684446e+000; foe(n+1)=5.053639e+001; krok(n+1)=7.486185e-006; ng(n+1)=3.469432e+002;
n=294; farx(n+1)=6.707416e+000; foe(n+1)=5.051320e+001; krok(n+1)=1.663763e-005; ng(n+1)=1.932303e+002;
n=295; farx(n+1)=6.673656e+000; foe(n+1)=5.049018e+001; krok(n+1)=7.482921e-006; ng(n+1)=3.446598e+002;
n=296; farx(n+1)=6.696367e+000; foe(n+1)=5.046747e+001; krok(n+1)=1.656922e-005; ng(n+1)=1.919224e+002;
n=297; farx(n+1)=6.663149e+000; foe(n+1)=5.044492e+001; krok(n+1)=7.486185e-006; ng(n+1)=3.411990e+002;
n=298; farx(n+1)=6.685589e+000; foe(n+1)=5.042264e+001; krok(n+1)=1.654527e-005; ng(n+1)=1.905309e+002;
n=299; farx(n+1)=6.652922e+000; foe(n+1)=5.040051e+001; krok(n+1)=7.475688e-006; ng(n+1)=3.380743e+002;
n=300; farx(n+1)=6.675118e+000; foe(n+1)=5.037856e+001; krok(n+1)=1.663763e-005; ng(n+1)=1.889730e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
