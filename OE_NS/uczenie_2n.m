%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.375264e+003; foe(n+1)=4.369417e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=3.973497e+003; foe(n+1)=3.985082e+003; krok(n+1)=5.333277e-004; ng(n+1)=1.350474e+003;
n=2; farx(n+1)=1.130327e+003; foe(n+1)=1.187399e+003; krok(n+1)=6.607631e-003; ng(n+1)=6.597537e+002;
n=3; farx(n+1)=6.999424e+002; foe(n+1)=8.974786e+002; krok(n+1)=1.389822e-004; ng(n+1)=3.550302e+003;
n=4; farx(n+1)=7.292628e+002; foe(n+1)=8.196427e+002; krok(n+1)=3.885492e-004; ng(n+1)=1.112350e+003;
n=5; farx(n+1)=5.881784e+002; foe(n+1)=7.658883e+002; krok(n+1)=1.163141e-004; ng(n+1)=1.794551e+003;
n=6; farx(n+1)=5.894646e+002; foe(n+1)=7.222397e+002; krok(n+1)=2.952315e-004; ng(n+1)=7.478293e+002;
n=7; farx(n+1)=4.933431e+002; foe(n+1)=6.814955e+002; krok(n+1)=9.968957e-005; ng(n+1)=1.555557e+003;
n=8; farx(n+1)=4.947750e+002; foe(n+1)=6.455882e+002; krok(n+1)=2.695564e-004; ng(n+1)=7.835204e+002;
n=9; farx(n+1)=4.194907e+002; foe(n+1)=6.104493e+002; krok(n+1)=8.683747e-005; ng(n+1)=1.523349e+003;
n=10; farx(n+1)=4.157262e+002; foe(n+1)=5.783125e+002; krok(n+1)=2.596664e-004; ng(n+1)=7.664275e+002;
n=11; farx(n+1)=3.561260e+002; foe(n+1)=5.473927e+002; krok(n+1)=7.639254e-005; ng(n+1)=1.525793e+003;
n=12; farx(n+1)=3.513729e+002; foe(n+1)=5.201619e+002; krok(n+1)=2.370496e-004; ng(n+1)=7.649877e+002;
n=13; farx(n+1)=3.051120e+002; foe(n+1)=4.940949e+002; krok(n+1)=6.831655e-005; ng(n+1)=1.478755e+003;
n=14; farx(n+1)=3.008393e+002; foe(n+1)=4.715971e+002; krok(n+1)=2.114007e-004; ng(n+1)=7.582763e+002;
n=15; farx(n+1)=2.650233e+002; foe(n+1)=4.501941e+002; krok(n+1)=6.197037e-005; ng(n+1)=1.399378e+003;
n=16; farx(n+1)=2.614912e+002; foe(n+1)=4.316417e+002; krok(n+1)=1.895397e-004; ng(n+1)=7.470928e+002;
n=17; farx(n+1)=2.335767e+002; foe(n+1)=4.140906e+002; krok(n+1)=5.618497e-005; ng(n+1)=1.324971e+003;
n=18; farx(n+1)=2.297049e+002; foe(n+1)=3.978060e+002; krok(n+1)=1.848518e-004; ng(n+1)=7.265108e+002;
n=19; farx(n+1)=2.069482e+002; foe(n+1)=3.827166e+002; krok(n+1)=5.111781e-005; ng(n+1)=1.311072e+003;
n=20; farx(n+1)=2.035621e+002; foe(n+1)=3.689511e+002; krok(n+1)=1.736749e-004; ng(n+1)=7.054082e+002;
n=21; farx(n+1)=1.851121e+002; foe(n+1)=3.562976e+002; krok(n+1)=4.687315e-005; ng(n+1)=1.269188e+003;
n=22; farx(n+1)=1.821153e+002; foe(n+1)=3.445076e+002; krok(n+1)=1.662532e-004; ng(n+1)=6.811739e+002;
n=23; farx(n+1)=1.668999e+002; foe(n+1)=3.336856e+002; krok(n+1)=4.279756e-005; ng(n+1)=1.241612e+003;
n=24; farx(n+1)=1.639297e+002; foe(n+1)=3.231709e+002; krok(n+1)=1.662532e-004; ng(n+1)=6.545789e+002;
n=25; farx(n+1)=1.509869e+002; foe(n+1)=3.135684e+002; krok(n+1)=3.909086e-005; ng(n+1)=1.244397e+003;
n=26; farx(n+1)=1.484162e+002; foe(n+1)=3.044262e+002; krok(n+1)=1.601067e-004; ng(n+1)=6.313723e+002;
n=27; farx(n+1)=1.374242e+002; foe(n+1)=2.960296e+002; krok(n+1)=3.594329e-005; ng(n+1)=1.235525e+003;
n=28; farx(n+1)=1.353181e+002; foe(n+1)=2.881068e+002; krok(n+1)=1.502858e-004; ng(n+1)=6.101284e+002;
n=29; farx(n+1)=1.259757e+002; foe(n+1)=2.808546e+002; krok(n+1)=3.331725e-005; ng(n+1)=1.211390e+003;
n=30; farx(n+1)=1.243331e+002; foe(n+1)=2.739974e+002; krok(n+1)=1.395032e-004; ng(n+1)=5.915029e+002;
n=31; farx(n+1)=1.163453e+002; foe(n+1)=2.676724e+002; krok(n+1)=3.079431e-005; ng(n+1)=1.188152e+003;
n=32; farx(n+1)=1.147915e+002; foe(n+1)=2.613767e+002; krok(n+1)=1.376014e-004; ng(n+1)=5.758036e+002;
n=33; farx(n+1)=1.076609e+002; foe(n+1)=2.555470e+002; krok(n+1)=2.860122e-005; ng(n+1)=1.199086e+003;
n=34; farx(n+1)=1.065866e+002; foe(n+1)=2.501705e+002; krok(n+1)=1.222111e-004; ng(n+1)=5.663073e+002;
n=35; farx(n+1)=1.004694e+002; foe(n+1)=2.451940e+002; krok(n+1)=2.697415e-005; ng(n+1)=1.153110e+003;
n=36; farx(n+1)=9.965609e+001; foe(n+1)=2.404372e+002; krok(n+1)=1.131625e-004; ng(n+1)=5.572819e+002;
n=37; farx(n+1)=9.425535e+001; foe(n+1)=2.359838e+002; krok(n+1)=2.533263e-005; ng(n+1)=1.133828e+003;
n=38; farx(n+1)=9.356990e+001; foe(n+1)=2.316536e+002; krok(n+1)=1.070782e-004; ng(n+1)=5.482377e+002;
n=39; farx(n+1)=8.874822e+001; foe(n+1)=2.276096e+002; krok(n+1)=2.383042e-005; ng(n+1)=1.120777e+003;
n=40; farx(n+1)=8.816401e+001; foe(n+1)=2.236440e+002; krok(n+1)=1.018184e-004; ng(n+1)=5.405382e+002;
n=41; farx(n+1)=8.385500e+001; foe(n+1)=2.199417e+002; krok(n+1)=2.230894e-005; ng(n+1)=1.110195e+003;
n=42; farx(n+1)=8.320311e+001; foe(n+1)=2.161073e+002; krok(n+1)=1.042068e-004; ng(n+1)=5.334354e+002;
n=43; farx(n+1)=7.919735e+001; foe(n+1)=2.124912e+002; krok(n+1)=2.078166e-005; ng(n+1)=1.137695e+003;
n=44; farx(n+1)=7.856479e+001; foe(n+1)=2.088769e+002; krok(n+1)=1.019540e-004; ng(n+1)=5.258315e+002;
n=45; farx(n+1)=7.492684e+001; foe(n+1)=2.055002e+002; krok(n+1)=1.954543e-005; ng(n+1)=1.136242e+003;
n=46; farx(n+1)=7.433438e+001; foe(n+1)=2.021139e+002; krok(n+1)=9.926586e-005; ng(n+1)=5.179931e+002;
n=47; farx(n+1)=7.101330e+001; foe(n+1)=1.989594e+002; krok(n+1)=1.850362e-005; ng(n+1)=1.133013e+003;
n=48; farx(n+1)=7.053056e+001; foe(n+1)=1.958788e+002; krok(n+1)=9.345155e-005; ng(n+1)=5.100405e+002;
n=49; farx(n+1)=6.752769e+001; foe(n+1)=1.929970e+002; krok(n+1)=1.763687e-005; ng(n+1)=1.114391e+003;
n=50; farx(n+1)=6.713363e+001; foe(n+1)=1.901794e+002; krok(n+1)=8.876552e-005; ng(n+1)=5.022094e+002;
n=51; farx(n+1)=6.440453e+001; foe(n+1)=1.875133e+002; krok(n+1)=1.675917e-005; ng(n+1)=1.100926e+003;
n=52; farx(n+1)=6.400417e+001; foe(n+1)=1.848106e+002; krok(n+1)=8.876552e-005; ng(n+1)=4.943946e+002;
n=53; farx(n+1)=6.145736e+001; foe(n+1)=1.822496e+002; krok(n+1)=1.596631e-005; ng(n+1)=1.107171e+003;
n=54; farx(n+1)=6.113419e+001; foe(n+1)=1.797712e+002; krok(n+1)=8.408970e-005; ng(n+1)=4.865137e+002;
n=55; farx(n+1)=5.880466e+001; foe(n+1)=1.774124e+002; krok(n+1)=1.535151e-005; ng(n+1)=1.088345e+003;
n=56; farx(n+1)=5.855675e+001; foe(n+1)=1.751456e+002; krok(n+1)=7.872898e-005; ng(n+1)=4.789026e+002;
n=57; farx(n+1)=5.643642e+001; foe(n+1)=1.730011e+002; krok(n+1)=1.482745e-005; ng(n+1)=1.060945e+003;
n=58; farx(n+1)=5.624082e+001; foe(n+1)=1.709014e+002; krok(n+1)=7.546424e-005; ng(n+1)=4.716322e+002;
n=59; farx(n+1)=5.429240e+001; foe(n+1)=1.688930e+002; krok(n+1)=1.416371e-005; ng(n+1)=1.048054e+003;
n=60; farx(n+1)=5.402549e+001; foe(n+1)=1.667571e+002; krok(n+1)=8.074766e-005; ng(n+1)=4.643290e+002;
n=61; farx(n+1)=5.212610e+001; foe(n+1)=1.647045e+002; krok(n+1)=1.356957e-005; ng(n+1)=1.080201e+003;
n=62; farx(n+1)=5.194360e+001; foe(n+1)=1.627747e+002; krok(n+1)=7.447980e-005; ng(n+1)=4.570018e+002;
n=63; farx(n+1)=5.021607e+001; foe(n+1)=1.609256e+002; krok(n+1)=1.318762e-005; ng(n+1)=1.045566e+003;
n=64; farx(n+1)=5.007418e+001; foe(n+1)=1.591240e+002; krok(n+1)=7.148615e-005; ng(n+1)=4.501607e+002;
n=65; farx(n+1)=4.846765e+001; foe(n+1)=1.573937e+002; krok(n+1)=1.278344e-005; ng(n+1)=1.029354e+003;
n=66; farx(n+1)=4.835065e+001; foe(n+1)=1.556881e+002; krok(n+1)=6.952670e-005; ng(n+1)=4.436154e+002;
n=67; farx(n+1)=4.683993e+001; foe(n+1)=1.540518e+002; krok(n+1)=1.245057e-005; ng(n+1)=1.017544e+003;
n=68; farx(n+1)=4.677644e+001; foe(n+1)=1.524878e+002; krok(n+1)=6.531497e-005; ng(n+1)=4.374713e+002;
n=69; farx(n+1)=4.537183e+001; foe(n+1)=1.509734e+002; krok(n+1)=1.216245e-005; ng(n+1)=9.942571e+002;
n=70; farx(n+1)=4.533769e+001; foe(n+1)=1.495050e+002; krok(n+1)=6.257477e-005; ng(n+1)=4.317097e+002;
n=71; farx(n+1)=4.401782e+001; foe(n+1)=1.480860e+002; krok(n+1)=1.191391e-005; ng(n+1)=9.762932e+002;
n=72; farx(n+1)=4.402131e+001; foe(n+1)=1.467257e+002; krok(n+1)=5.904967e-005; ng(n+1)=4.290214e+002;
n=73; farx(n+1)=4.278689e+001; foe(n+1)=1.454061e+002; krok(n+1)=1.167393e-005; ng(n+1)=9.541474e+002;
n=74; farx(n+1)=4.280744e+001; foe(n+1)=1.441109e+002; krok(n+1)=5.770216e-005; ng(n+1)=4.259299e+002;
n=75; farx(n+1)=4.163297e+001; foe(n+1)=1.428441e+002; krok(n+1)=1.138543e-005; ng(n+1)=9.463297e+002;
n=76; farx(n+1)=4.165708e+001; foe(n+1)=1.415859e+002; krok(n+1)=5.735983e-005; ng(n+1)=4.219589e+002;
n=77; farx(n+1)=4.053328e+001; foe(n+1)=1.403589e+002; krok(n+1)=1.109774e-005; ng(n+1)=9.422001e+002;
n=78; farx(n+1)=4.056002e+001; foe(n+1)=1.391308e+002; krok(n+1)=5.770216e-005; ng(n+1)=4.179246e+002;
n=79; farx(n+1)=3.948391e+001; foe(n+1)=1.379239e+002; krok(n+1)=1.069939e-005; ng(n+1)=9.447415e+002;
n=80; farx(n+1)=3.945901e+001; foe(n+1)=1.366098e+002; krok(n+1)=6.408173e-005; ng(n+1)=4.094775e+002;
n=81; farx(n+1)=3.837275e+001; foe(n+1)=1.353348e+002; krok(n+1)=1.038706e-005; ng(n+1)=9.796005e+002;
n=82; farx(n+1)=3.838228e+001; foe(n+1)=1.341146e+002; krok(n+1)=6.065098e-005; ng(n+1)=4.077966e+002;
n=83; farx(n+1)=3.736159e+001; foe(n+1)=1.329251e+002; krok(n+1)=1.018440e-005; ng(n+1)=9.576613e+002;
n=84; farx(n+1)=3.738119e+001; foe(n+1)=1.317444e+002; krok(n+1)=6.013826e-005; ng(n+1)=4.044952e+002;
n=85; farx(n+1)=3.639629e+001; foe(n+1)=1.305911e+002; krok(n+1)=9.999543e-006; ng(n+1)=9.526767e+002;
n=86; farx(n+1)=3.644603e+001; foe(n+1)=1.294939e+002; krok(n+1)=5.720244e-005; ng(n+1)=4.032492e+002;
n=87; farx(n+1)=3.551293e+001; foe(n+1)=1.284061e+002; krok(n+1)=9.808732e-006; ng(n+1)=9.351352e+002;
n=88; farx(n+1)=3.556473e+001; foe(n+1)=1.273320e+002; krok(n+1)=5.665485e-005; ng(n+1)=3.997851e+002;
n=89; farx(n+1)=3.466886e+001; foe(n+1)=1.262888e+002; krok(n+1)=9.645404e-006; ng(n+1)=9.246369e+002;
n=90; farx(n+1)=3.473986e+001; foe(n+1)=1.252705e+002; krok(n+1)=5.538211e-005; ng(n+1)=3.979329e+002;
n=91; farx(n+1)=3.387542e+001; foe(n+1)=1.242590e+002; krok(n+1)=9.456672e-006; ng(n+1)=9.186158e+002;
n=92; farx(n+1)=3.395101e+001; foe(n+1)=1.232662e+002; krok(n+1)=5.491616e-005; ng(n+1)=3.948892e+002;
n=93; farx(n+1)=3.311632e+001; foe(n+1)=1.222893e+002; krok(n+1)=9.301319e-006; ng(n+1)=9.107018e+002;
n=94; farx(n+1)=3.320414e+001; foe(n+1)=1.213420e+002; krok(n+1)=5.314479e-005; ng(n+1)=3.928994e+002;
n=95; farx(n+1)=3.240575e+001; foe(n+1)=1.204159e+002; krok(n+1)=9.173490e-006; ng(n+1)=8.947732e+002;
n=96; farx(n+1)=3.250649e+001; foe(n+1)=1.195110e+002; krok(n+1)=5.210339e-005; ng(n+1)=3.909720e+002;
n=97; farx(n+1)=3.173361e+001; foe(n+1)=1.186143e+002; krok(n+1)=9.024321e-006; ng(n+1)=8.875090e+002;
n=98; farx(n+1)=3.184080e+001; foe(n+1)=1.177395e+002; krok(n+1)=5.121654e-005; ng(n+1)=3.885486e+002;
n=99; farx(n+1)=3.109338e+001; foe(n+1)=1.168778e+002; krok(n+1)=8.913402e-006; ng(n+1)=8.769925e+002;
n=100; farx(n+1)=3.121335e+001; foe(n+1)=1.160530e+002; krok(n+1)=4.886973e-005; ng(n+1)=3.873274e+002;
n=101; farx(n+1)=3.050057e+001; foe(n+1)=1.152442e+002; krok(n+1)=8.818436e-006; ng(n+1)=8.568640e+002;
n=102; farx(n+1)=3.063076e+001; foe(n+1)=1.144530e+002; krok(n+1)=4.837788e-005; ng(n+1)=3.854050e+002;
n=103; farx(n+1)=2.993484e+001; foe(n+1)=1.136599e+002; krok(n+1)=8.667355e-006; ng(n+1)=8.541912e+002;
n=104; farx(n+1)=3.006358e+001; foe(n+1)=1.128826e+002; krok(n+1)=4.789661e-005; ng(n+1)=3.823776e+002;
n=105; farx(n+1)=2.939288e+001; foe(n+1)=1.121235e+002; krok(n+1)=8.546025e-006; ng(n+1)=8.418765e+002;
n=106; farx(n+1)=2.952352e+001; foe(n+1)=1.113624e+002; krok(n+1)=4.789661e-005; ng(n+1)=3.795003e+002;
n=107; farx(n+1)=2.886721e+001; foe(n+1)=1.106217e+002; krok(n+1)=8.445281e-006; ng(n+1)=8.370267e+002;
n=108; farx(n+1)=2.900963e+001; foe(n+1)=1.099034e+002; krok(n+1)=4.645269e-005; ng(n+1)=3.782340e+002;
n=109; farx(n+1)=2.837349e+001; foe(n+1)=1.091908e+002; krok(n+1)=8.342260e-006; ng(n+1)=8.266421e+002;
n=110; farx(n+1)=2.851892e+001; foe(n+1)=1.084970e+002; krok(n+1)=4.554170e-005; ng(n+1)=3.759500e+002;
n=111; farx(n+1)=2.790251e+001; foe(n+1)=1.078151e+002; krok(n+1)=8.269595e-006; ng(n+1)=8.141387e+002;
n=112; farx(n+1)=2.805398e+001; foe(n+1)=1.071598e+002; krok(n+1)=4.341873e-005; ng(n+1)=3.746770e+002;
n=113; farx(n+1)=2.746544e+001; foe(n+1)=1.065223e+002; krok(n+1)=8.206521e-006; ng(n+1)=7.929975e+002;
n=114; farx(n+1)=2.762330e+001; foe(n+1)=1.058903e+002; krok(n+1)=4.334515e-005; ng(n+1)=3.725328e+002;
n=115; farx(n+1)=2.704233e+001; foe(n+1)=1.052607e+002; krok(n+1)=8.106422e-006; ng(n+1)=7.922481e+002;
n=116; farx(n+1)=2.720368e+001; foe(n+1)=1.046569e+002; krok(n+1)=4.204485e-005; ng(n+1)=3.707514e+002;
n=117; farx(n+1)=2.664473e+001; foe(n+1)=1.040584e+002; krok(n+1)=8.010217e-006; ng(n+1)=7.773883e+002;
n=118; farx(n+1)=2.680491e+001; foe(n+1)=1.034595e+002; krok(n+1)=4.285449e-005; ng(n+1)=3.672772e+002;
n=119; farx(n+1)=2.625193e+001; foe(n+1)=1.028658e+002; krok(n+1)=7.902407e-006; ng(n+1)=7.779927e+002;
n=120; farx(n+1)=2.641328e+001; foe(n+1)=1.022876e+002; krok(n+1)=4.211190e-005; ng(n+1)=3.648949e+002;
n=121; farx(n+1)=2.587705e+001; foe(n+1)=1.017184e+002; krok(n+1)=7.821846e-006; ng(n+1)=7.663117e+002;
n=122; farx(n+1)=2.604042e+001; foe(n+1)=1.011569e+002; krok(n+1)=4.204485e-005; ng(n+1)=3.622445e+002;
n=123; farx(n+1)=2.551597e+001; foe(n+1)=1.006003e+002; krok(n+1)=7.698578e-006; ng(n+1)=7.616711e+002;
n=124; farx(n+1)=2.567767e+001; foe(n+1)=1.000364e+002; krok(n+1)=4.381503e-005; ng(n+1)=3.581115e+002;
n=125; farx(n+1)=2.515325e+001; foe(n+1)=9.947243e+001; krok(n+1)=7.559984e-006; ng(n+1)=7.693795e+002;
n=126; farx(n+1)=2.531089e+001; foe(n+1)=9.891345e+001; krok(n+1)=4.418089e-005; ng(n+1)=3.543886e+002;
n=127; farx(n+1)=2.479825e+001; foe(n+1)=9.836640e+001; krok(n+1)=7.475688e-006; ng(n+1)=7.618028e+002;
n=128; farx(n+1)=2.495802e+001; foe(n+1)=9.782660e+001; krok(n+1)=4.396667e-005; ng(n+1)=3.516162e+002;
n=129; farx(n+1)=2.445366e+001; foe(n+1)=9.729422e+001; krok(n+1)=7.425417e-006; ng(n+1)=7.552937e+002;
n=130; farx(n+1)=2.461785e+001; foe(n+1)=9.679561e+001; krok(n+1)=4.120586e-005; ng(n+1)=3.505439e+002;
n=131; farx(n+1)=2.413674e+001; foe(n+1)=9.630515e+001; krok(n+1)=7.425417e-006; ng(n+1)=7.301832e+002;
n=132; farx(n+1)=2.430520e+001; foe(n+1)=9.584464e+001; krok(n+1)=3.870095e-005; ng(n+1)=3.496932e+002;
n=133; farx(n+1)=2.384448e+001; foe(n+1)=9.538965e+001; krok(n+1)=7.425417e-006; ng(n+1)=7.077889e+002;
n=134; farx(n+1)=2.401511e+001; foe(n+1)=9.496053e+001; krok(n+1)=3.655818e-005; ng(n+1)=3.486803e+002;
n=135; farx(n+1)=2.357563e+001; foe(n+1)=9.453692e+001; krok(n+1)=7.381362e-006; ng(n+1)=6.868842e+002;
n=136; farx(n+1)=2.374506e+001; foe(n+1)=9.411403e+001; krok(n+1)=3.704903e-005; ng(n+1)=3.453765e+002;
n=137; farx(n+1)=2.331239e+001; foe(n+1)=9.369609e+001; krok(n+1)=7.265431e-006; ng(n+1)=6.848194e+002;
n=138; farx(n+1)=2.347879e+001; foe(n+1)=9.326766e+001; krok(n+1)=3.882721e-005; ng(n+1)=3.409054e+002;
n=139; farx(n+1)=2.304473e+001; foe(n+1)=9.284446e+001; krok(n+1)=7.156975e-006; ng(n+1)=6.908923e+002;
n=140; farx(n+1)=2.320955e+001; foe(n+1)=9.242402e+001; krok(n+1)=3.909086e-005; ng(n+1)=3.376056e+002;
n=141; farx(n+1)=2.278081e+001; foe(n+1)=9.201007e+001; krok(n+1)=7.114163e-006; ng(n+1)=6.858575e+002;
n=142; farx(n+1)=2.294721e+001; foe(n+1)=9.161619e+001; krok(n+1)=3.732546e-005; ng(n+1)=3.360824e+002;
n=143; farx(n+1)=2.253550e+001; foe(n+1)=9.122700e+001; krok(n+1)=7.065934e-006; ng(n+1)=6.684376e+002;
n=144; farx(n+1)=2.270048e+001; foe(n+1)=9.083871e+001; krok(n+1)=3.782543e-005; ng(n+1)=3.326892e+002;
n=145; farx(n+1)=2.229312e+001; foe(n+1)=9.045504e+001; krok(n+1)=6.997257e-006; ng(n+1)=6.657142e+002;
n=146; farx(n+1)=2.245740e+001; foe(n+1)=9.007897e+001; krok(n+1)=3.750283e-005; ng(n+1)=3.299409e+002;
n=147; farx(n+1)=2.205879e+001; foe(n+1)=8.970746e+001; krok(n+1)=6.939469e-006; ng(n+1)=6.574462e+002;
n=148; farx(n+1)=2.222205e+001; foe(n+1)=8.934081e+001; krok(n+1)=3.746225e-005; ng(n+1)=3.269681e+002;
n=149; farx(n+1)=2.183082e+001; foe(n+1)=8.897881e+001; krok(n+1)=6.872937e-006; ng(n+1)=6.510499e+002;
n=150; farx(n+1)=2.199288e+001; foe(n+1)=8.861836e+001; krok(n+1)=3.782669e-005; ng(n+1)=3.236992e+002;
n=151; farx(n+1)=2.160631e+001; foe(n+1)=8.826217e+001; krok(n+1)=6.807713e-006; ng(n+1)=6.474860e+002;
n=152; farx(n+1)=2.176726e+001; foe(n+1)=8.791089e+001; krok(n+1)=3.773938e-005; ng(n+1)=3.208025e+002;
n=153; farx(n+1)=2.138744e+001; foe(n+1)=8.756415e+001; krok(n+1)=6.755832e-006; ng(n+1)=6.407084e+002;
n=154; farx(n+1)=2.154770e+001; foe(n+1)=8.722352e+001; krok(n+1)=3.746225e-005; ng(n+1)=3.181470e+002;
n=155; farx(n+1)=2.117461e+001; foe(n+1)=8.688685e+001; krok(n+1)=6.715934e-006; ng(n+1)=6.331508e+002;
n=156; farx(n+1)=2.133448e+001; foe(n+1)=8.656007e+001; krok(n+1)=3.669396e-005; ng(n+1)=3.158422e+002;
n=157; farx(n+1)=2.097065e+001; foe(n+1)=8.623650e+001; krok(n+1)=6.676010e-006; ng(n+1)=6.226660e+002;
n=158; farx(n+1)=2.112954e+001; foe(n+1)=8.591861e+001; krok(n+1)=3.650218e-005; ng(n+1)=3.131715e+002;
n=159; farx(n+1)=2.077197e+001; foe(n+1)=8.560370e+001; krok(n+1)=6.631510e-006; ng(n+1)=6.157872e+002;
n=160; farx(n+1)=2.092971e+001; foe(n+1)=8.529530e+001; krok(n+1)=3.613934e-005; ng(n+1)=3.106012e+002;
n=161; farx(n+1)=2.057928e+001; foe(n+1)=8.499015e+001; krok(n+1)=6.589629e-006; ng(n+1)=6.076852e+002;
n=162; farx(n+1)=2.073611e+001; foe(n+1)=8.468989e+001; krok(n+1)=3.598915e-005; ng(n+1)=3.079535e+002;
n=163; farx(n+1)=2.039141e+001; foe(n+1)=8.439222e+001; krok(n+1)=6.542741e-006; ng(n+1)=6.014509e+002;
n=164; farx(n+1)=2.054658e+001; foe(n+1)=8.409896e+001; krok(n+1)=3.584012e-005; ng(n+1)=3.052271e+002;
n=165; farx(n+1)=2.020772e+001; foe(n+1)=8.380935e+001; krok(n+1)=6.505235e-006; ng(n+1)=5.944887e+002;
n=166; farx(n+1)=2.036232e+001; foe(n+1)=8.352525e+001; krok(n+1)=3.551969e-005; ng(n+1)=3.027551e+002;
n=167; farx(n+1)=2.002926e+001; foe(n+1)=8.324324e+001; krok(n+1)=6.464322e-006; ng(n+1)=5.877541e+002;
n=168; farx(n+1)=2.018231e+001; foe(n+1)=8.296615e+001; krok(n+1)=3.531947e-005; ng(n+1)=3.001121e+002;
n=169; farx(n+1)=1.985489e+001; foe(n+1)=8.269177e+001; krok(n+1)=6.429293e-006; ng(n+1)=5.808301e+002;
n=170; farx(n+1)=2.000693e+001; foe(n+1)=8.242291e+001; krok(n+1)=3.498787e-005; ng(n+1)=2.976377e+002;
n=171; farx(n+1)=1.968628e+001; foe(n+1)=8.215607e+001; krok(n+1)=6.371971e-006; ng(n+1)=5.741953e+002;
n=172; farx(n+1)=1.983629e+001; foe(n+1)=8.188586e+001; krok(n+1)=3.606196e-005; ng(n+1)=2.942021e+002;
n=173; farx(n+1)=1.951523e+001; foe(n+1)=8.161838e+001; krok(n+1)=6.325379e-006; ng(n+1)=5.768614e+002;
n=174; farx(n+1)=1.966430e+001; foe(n+1)=8.135853e+001; krok(n+1)=3.535372e-005; ng(n+1)=2.919251e+002;
n=175; farx(n+1)=1.935077e+001; foe(n+1)=8.110057e+001; krok(n+1)=6.289522e-006; ng(n+1)=5.680805e+002;
n=176; farx(n+1)=1.949851e+001; foe(n+1)=8.084459e+001; krok(n+1)=3.565361e-005; ng(n+1)=2.890464e+002;
n=177; farx(n+1)=1.918809e+001; foe(n+1)=8.058996e+001; krok(n+1)=6.240680e-006; ng(n+1)=5.660983e+002;
n=178; farx(n+1)=1.933406e+001; foe(n+1)=8.033799e+001; krok(n+1)=3.583315e-005; ng(n+1)=2.861817e+002;
n=179; farx(n+1)=1.902821e+001; foe(n+1)=8.008784e+001; krok(n+1)=6.182960e-006; ng(n+1)=5.626863e+002;
n=180; farx(n+1)=1.917216e+001; foe(n+1)=7.983458e+001; krok(n+1)=3.693653e-005; ng(n+1)=2.827805e+002;
n=181; farx(n+1)=1.886597e+001; foe(n+1)=7.958360e+001; krok(n+1)=6.137476e-006; ng(n+1)=5.652774e+002;
n=182; farx(n+1)=1.900853e+001; foe(n+1)=7.933874e+001; krok(n+1)=3.639457e-005; ng(n+1)=2.803123e+002;
n=183; farx(n+1)=1.870838e+001; foe(n+1)=7.909596e+001; krok(n+1)=6.116772e-006; ng(n+1)=5.573946e+002;
n=184; farx(n+1)=1.884928e+001; foe(n+1)=7.885987e+001; krok(n+1)=3.565785e-005; ng(n+1)=2.779669e+002;
n=185; farx(n+1)=1.855612e+001; foe(n+1)=7.862636e+001; krok(n+1)=6.098994e-006; ng(n+1)=5.479512e+002;
n=186; farx(n+1)=1.869583e+001; foe(n+1)=7.839770e+001; krok(n+1)=3.521996e-005; ng(n+1)=2.755754e+002;
n=187; farx(n+1)=1.840762e+001; foe(n+1)=7.817085e+001; krok(n+1)=6.082025e-006; ng(n+1)=5.413357e+002;
n=188; farx(n+1)=1.854495e+001; foe(n+1)=7.795106e+001; krok(n+1)=3.418410e-005; ng(n+1)=2.734217e+002;
n=189; farx(n+1)=1.826487e+001; foe(n+1)=7.773515e+001; krok(n+1)=6.082025e-006; ng(n+1)=5.291425e+002;
n=190; farx(n+1)=1.840200e+001; foe(n+1)=7.752439e+001; krok(n+1)=3.351433e-005; ng(n+1)=2.714289e+002;
n=191; farx(n+1)=1.812771e+001; foe(n+1)=7.731464e+001; krok(n+1)=6.043919e-006; ng(n+1)=5.226046e+002;
n=192; farx(n+1)=1.826261e+001; foe(n+1)=7.710377e+001; krok(n+1)=3.418410e-005; ng(n+1)=2.684065e+002;
n=193; farx(n+1)=1.798976e+001; foe(n+1)=7.689523e+001; krok(n+1)=6.005855e-006; ng(n+1)=5.222741e+002;
n=194; farx(n+1)=1.812435e+001; foe(n+1)=7.668867e+001; krok(n+1)=3.440034e-005; ng(n+1)=2.658998e+002;
n=195; farx(n+1)=1.785299e+001; foe(n+1)=7.648155e+001; krok(n+1)=5.972097e-006; ng(n+1)=5.216114e+002;
n=196; farx(n+1)=1.798517e+001; foe(n+1)=7.628004e+001; krok(n+1)=3.395536e-005; ng(n+1)=2.634996e+002;
n=197; farx(n+1)=1.771978e+001; foe(n+1)=7.607992e+001; krok(n+1)=5.952427e-006; ng(n+1)=5.136195e+002;
n=198; farx(n+1)=1.785044e+001; foe(n+1)=7.588213e+001; krok(n+1)=3.395536e-005; ng(n+1)=2.610033e+002;
n=199; farx(n+1)=1.758780e+001; foe(n+1)=7.568564e+001; krok(n+1)=5.939356e-006; ng(n+1)=5.098860e+002;
n=200; farx(n+1)=1.771680e+001; foe(n+1)=7.549610e+001; krok(n+1)=3.291660e-005; ng(n+1)=2.591179e+002;
n=201; farx(n+1)=1.746185e+001; foe(n+1)=7.530824e+001; krok(n+1)=5.921156e-006; ng(n+1)=4.992608e+002;
n=202; farx(n+1)=1.758955e+001; foe(n+1)=7.511978e+001; krok(n+1)=3.345851e-005; ng(n+1)=2.564189e+002;
n=203; farx(n+1)=1.733598e+001; foe(n+1)=7.493246e+001; krok(n+1)=5.876381e-006; ng(n+1)=5.000889e+002;
n=204; farx(n+1)=1.746237e+001; foe(n+1)=7.474441e+001; krok(n+1)=3.415827e-005; ng(n+1)=2.536346e+002;
n=205; farx(n+1)=1.720948e+001; foe(n+1)=7.455692e+001; krok(n+1)=5.834398e-006; ng(n+1)=5.026479e+002;
n=206; farx(n+1)=1.733390e+001; foe(n+1)=7.437048e+001; krok(n+1)=3.444097e-005; ng(n+1)=2.509895e+002;
n=207; farx(n+1)=1.708308e+001; foe(n+1)=7.418549e+001; krok(n+1)=5.819064e-006; ng(n+1)=5.011224e+002;
n=208; farx(n+1)=1.720650e+001; foe(n+1)=7.400609e+001; krok(n+1)=3.369455e-005; ng(n+1)=2.490443e+002;
n=209; farx(n+1)=1.696096e+001; foe(n+1)=7.382712e+001; krok(n+1)=5.807097e-006; ng(n+1)=4.940217e+002;
n=210; farx(n+1)=1.708232e+001; foe(n+1)=7.365172e+001; krok(n+1)=3.333298e-005; ng(n+1)=2.468045e+002;
n=211; farx(n+1)=1.684153e+001; foe(n+1)=7.347800e+001; krok(n+1)=5.793243e-006; ng(n+1)=4.878793e+002;
n=212; farx(n+1)=1.696192e+001; foe(n+1)=7.330610e+001; krok(n+1)=3.331725e-005; ng(n+1)=2.445625e+002;
n=213; farx(n+1)=1.672372e+001; foe(n+1)=7.313486e+001; krok(n+1)=5.769722e-006; ng(n+1)=4.859530e+002;
n=214; farx(n+1)=1.684245e+001; foe(n+1)=7.296588e+001; krok(n+1)=3.324441e-005; ng(n+1)=2.422779e+002;
n=215; farx(n+1)=1.660754e+001; foe(n+1)=7.279800e+001; krok(n+1)=5.751510e-006; ng(n+1)=4.824916e+002;
n=216; farx(n+1)=1.672491e+001; foe(n+1)=7.263224e+001; krok(n+1)=3.315087e-005; ng(n+1)=2.400714e+002;
n=217; farx(n+1)=1.649270e+001; foe(n+1)=7.246743e+001; krok(n+1)=5.742854e-006; ng(n+1)=4.793828e+002;
n=218; farx(n+1)=1.660802e+001; foe(n+1)=7.230734e+001; krok(n+1)=3.223481e-005; ng(n+1)=2.382538e+002;
n=219; farx(n+1)=1.638179e+001; foe(n+1)=7.214959e+001; krok(n+1)=5.751510e-006; ng(n+1)=4.695125e+002;
n=220; farx(n+1)=1.649668e+001; foe(n+1)=7.199514e+001; krok(n+1)=3.167303e-005; ng(n+1)=2.365499e+002;
n=221; farx(n+1)=1.627445e+001; foe(n+1)=7.184118e+001; krok(n+1)=5.731843e-006; ng(n+1)=4.648090e+002;
n=222; farx(n+1)=1.638805e+001; foe(n+1)=7.168767e+001; krok(n+1)=3.204809e-005; ng(n+1)=2.342210e+002;
n=223; farx(n+1)=1.616685e+001; foe(n+1)=7.153454e+001; krok(n+1)=5.714341e-006; ng(n+1)=4.650548e+002;
n=224; farx(n+1)=1.627908e+001; foe(n+1)=7.138510e+001; krok(n+1)=3.158272e-005; ng(n+1)=2.323709e+002;
n=225; farx(n+1)=1.606197e+001; foe(n+1)=7.123615e+001; krok(n+1)=5.705400e-006; ng(n+1)=4.594479e+002;
n=226; farx(n+1)=1.617312e+001; foe(n+1)=7.108932e+001; krok(n+1)=3.153604e-005; ng(n+1)=2.303630e+002;
n=227; farx(n+1)=1.595856e+001; foe(n+1)=7.094257e+001; krok(n+1)=5.685733e-006; ng(n+1)=4.571514e+002;
n=228; farx(n+1)=1.606819e+001; foe(n+1)=7.079737e+001; krok(n+1)=3.161211e-005; ng(n+1)=2.282406e+002;
n=229; farx(n+1)=1.585609e+001; foe(n+1)=7.065275e+001; krok(n+1)=5.667521e-006; ng(n+1)=4.549147e+002;
n=230; farx(n+1)=1.596441e+001; foe(n+1)=7.050939e+001; krok(n+1)=3.167304e-005; ng(n+1)=2.261762e+002;
n=231; farx(n+1)=1.575463e+001; foe(n+1)=7.036671e+001; krok(n+1)=5.649012e-006; ng(n+1)=4.529078e+002;
n=232; farx(n+1)=1.586203e+001; foe(n+1)=7.022508e+001; krok(n+1)=3.185663e-005; ng(n+1)=2.241209e+002;
n=233; farx(n+1)=1.565375e+001; foe(n+1)=7.008349e+001; krok(n+1)=5.630525e-006; ng(n+1)=4.523652e+002;
n=234; farx(n+1)=1.575986e+001; foe(n+1)=6.994422e+001; krok(n+1)=3.175337e-005; ng(n+1)=2.221842e+002;
n=235; farx(n+1)=1.555426e+001; foe(n+1)=6.980509e+001; krok(n+1)=5.618744e-006; ng(n+1)=4.492907e+002;
n=236; farx(n+1)=1.565905e+001; foe(n+1)=6.966821e+001; krok(n+1)=3.160963e-005; ng(n+1)=2.202889e+002;
n=237; farx(n+1)=1.545629e+001; foe(n+1)=6.953169e+001; krok(n+1)=5.607521e-006; ng(n+1)=4.458723e+002;
n=238; farx(n+1)=1.555988e+001; foe(n+1)=6.939694e+001; krok(n+1)=3.154696e-005; ng(n+1)=2.183931e+002;
n=239; farx(n+1)=1.535951e+001; foe(n+1)=6.926259e+001; krok(n+1)=5.597181e-006; ng(n+1)=4.431754e+002;
n=240; farx(n+1)=1.546185e+001; foe(n+1)=6.913035e+001; krok(n+1)=3.134037e-005; ng(n+1)=2.165939e+002;
n=241; farx(n+1)=1.526432e+001; foe(n+1)=6.899866e+001; krok(n+1)=5.589743e-006; ng(n+1)=4.394739e+002;
n=242; farx(n+1)=1.536575e+001; foe(n+1)=6.886887e+001; krok(n+1)=3.122446e-005; ng(n+1)=2.148248e+002;
n=243; farx(n+1)=1.517052e+001; foe(n+1)=6.873916e+001; krok(n+1)=5.578968e-006; ng(n+1)=4.369362e+002;
n=244; farx(n+1)=1.527055e+001; foe(n+1)=6.861149e+001; krok(n+1)=3.104372e-005; ng(n+1)=2.130414e+002;
n=245; farx(n+1)=1.507823e+001; foe(n+1)=6.848441e+001; krok(n+1)=5.570525e-006; ng(n+1)=4.331269e+002;
n=246; farx(n+1)=1.517725e+001; foe(n+1)=6.835855e+001; krok(n+1)=3.104372e-005; ng(n+1)=2.112443e+002;
n=247; farx(n+1)=1.498686e+001; foe(n+1)=6.823308e+001; krok(n+1)=5.560895e-006; ng(n+1)=4.311379e+002;
n=248; farx(n+1)=1.508527e+001; foe(n+1)=6.810957e+001; krok(n+1)=3.098519e-005; ng(n+1)=2.095627e+002;
n=249; farx(n+1)=1.489698e+001; foe(n+1)=6.798550e+001; krok(n+1)=5.540176e-006; ng(n+1)=4.294911e+002;
n=250; farx(n+1)=1.499369e+001; foe(n+1)=6.786206e+001; krok(n+1)=3.128738e-005; ng(n+1)=2.075304e+002;
n=251; farx(n+1)=1.480723e+001; foe(n+1)=6.773925e+001; krok(n+1)=5.523363e-006; ng(n+1)=4.280952e+002;
n=252; farx(n+1)=1.490300e+001; foe(n+1)=6.761675e+001; krok(n+1)=3.154696e-005; ng(n+1)=2.056527e+002;
n=253; farx(n+1)=1.471758e+001; foe(n+1)=6.749461e+001; krok(n+1)=5.513365e-006; ng(n+1)=4.277602e+002;
n=254; farx(n+1)=1.481244e+001; foe(n+1)=6.737490e+001; krok(n+1)=3.122446e-005; ng(n+1)=2.040982e+002;
n=255; farx(n+1)=1.462962e+001; foe(n+1)=6.725525e+001; krok(n+1)=5.513365e-006; ng(n+1)=4.238838e+002;
n=256; farx(n+1)=1.472316e+001; foe(n+1)=6.713845e+001; krok(n+1)=3.070447e-005; ng(n+1)=2.026150e+002;
n=257; farx(n+1)=1.454400e+001; foe(n+1)=6.702226e+001; krok(n+1)=5.513365e-006; ng(n+1)=4.180021e+002;
n=258; farx(n+1)=1.463646e+001; foe(n+1)=6.690717e+001; krok(n+1)=3.061523e-005; ng(n+1)=2.009792e+002;
n=259; farx(n+1)=1.445929e+001; foe(n+1)=6.679280e+001; krok(n+1)=5.513365e-006; ng(n+1)=4.152424e+002;
n=260; farx(n+1)=1.455119e+001; foe(n+1)=6.668074e+001; krok(n+1)=3.023762e-005; ng(n+1)=1.996097e+002;
n=261; farx(n+1)=1.437642e+001; foe(n+1)=6.656855e+001; krok(n+1)=5.513365e-006; ng(n+1)=4.116818e+002;
n=262; farx(n+1)=1.446702e+001; foe(n+1)=6.645894e+001; krok(n+1)=2.979489e-005; ng(n+1)=1.981856e+002;
n=263; farx(n+1)=1.429557e+001; foe(n+1)=6.634977e+001; krok(n+1)=5.513365e-006; ng(n+1)=4.063265e+002;
n=264; farx(n+1)=1.438544e+001; foe(n+1)=6.624159e+001; krok(n+1)=2.983290e-005; ng(n+1)=1.966202e+002;
n=265; farx(n+1)=1.421547e+001; foe(n+1)=6.613346e+001; krok(n+1)=5.504150e-006; ng(n+1)=4.049635e+002;
n=266; farx(n+1)=1.430421e+001; foe(n+1)=6.602686e+001; krok(n+1)=2.969568e-005; ng(n+1)=1.950991e+002;
n=267; farx(n+1)=1.413640e+001; foe(n+1)=6.592062e+001; krok(n+1)=5.504150e-006; ng(n+1)=4.017722e+002;
n=268; farx(n+1)=1.422427e+001; foe(n+1)=6.581623e+001; krok(n+1)=2.940339e-005; ng(n+1)=1.937327e+002;
n=269; farx(n+1)=1.405885e+001; foe(n+1)=6.571207e+001; krok(n+1)=5.504150e-006; ng(n+1)=3.981387e+002;
n=270; farx(n+1)=1.414543e+001; foe(n+1)=6.560954e+001; krok(n+1)=2.907853e-005; ng(n+1)=1.923386e+002;
n=271; farx(n+1)=1.398262e+001; foe(n+1)=6.550788e+001; krok(n+1)=5.513365e-006; ng(n+1)=3.935176e+002;
n=272; farx(n+1)=1.406840e+001; foe(n+1)=6.540809e+001; krok(n+1)=2.857275e-005; ng(n+1)=1.911641e+002;
n=273; farx(n+1)=1.390836e+001; foe(n+1)=6.530896e+001; krok(n+1)=5.519987e-006; ng(n+1)=3.887069e+002;
n=274; farx(n+1)=1.399437e+001; foe(n+1)=6.521123e+001; krok(n+1)=2.860122e-005; ng(n+1)=1.898977e+002;
n=275; farx(n+1)=1.383514e+001; foe(n+1)=6.511223e+001; krok(n+1)=5.499859e-006; ng(n+1)=3.891112e+002;
n=276; farx(n+1)=1.391946e+001; foe(n+1)=6.501470e+001; krok(n+1)=2.867992e-005; ng(n+1)=1.882609e+002;
n=277; farx(n+1)=1.376224e+001; foe(n+1)=6.491737e+001; krok(n+1)=5.499859e-006; ng(n+1)=3.861008e+002;
n=278; farx(n+1)=1.384525e+001; foe(n+1)=6.482156e+001; krok(n+1)=2.832743e-005; ng(n+1)=1.869567e+002;
n=279; farx(n+1)=1.369061e+001; foe(n+1)=6.472670e+001; krok(n+1)=5.513365e-006; ng(n+1)=3.812202e+002;
n=280; farx(n+1)=1.377301e+001; foe(n+1)=6.463361e+001; krok(n+1)=2.780862e-005; ng(n+1)=1.859147e+002;
n=281; farx(n+1)=1.362092e+001; foe(n+1)=6.454101e+001; krok(n+1)=5.519987e-006; ng(n+1)=3.767293e+002;
n=282; farx(n+1)=1.370222e+001; foe(n+1)=6.444966e+001; krok(n+1)=2.746916e-005; ng(n+1)=1.847016e+002;
n=283; farx(n+1)=1.355271e+001; foe(n+1)=6.435923e+001; krok(n+1)=5.523363e-006; ng(n+1)=3.723271e+002;
n=284; farx(n+1)=1.363412e+001; foe(n+1)=6.426909e+001; krok(n+1)=2.769106e-005; ng(n+1)=1.834041e+002;
n=285; farx(n+1)=1.348487e+001; foe(n+1)=6.417833e+001; krok(n+1)=5.504150e-006; ng(n+1)=3.736771e+002;
n=286; farx(n+1)=1.356427e+001; foe(n+1)=6.408878e+001; krok(n+1)=2.746916e-005; ng(n+1)=1.819656e+002;
n=287; farx(n+1)=1.341777e+001; foe(n+1)=6.400063e+001; krok(n+1)=5.521422e-006; ng(n+1)=3.681665e+002;
n=288; farx(n+1)=1.349701e+001; foe(n+1)=6.391352e+001; krok(n+1)=2.711505e-005; ng(n+1)=1.809924e+002;
n=289; farx(n+1)=1.335225e+001; foe(n+1)=6.382672e+001; krok(n+1)=5.521422e-006; ng(n+1)=3.655184e+002;
n=290; farx(n+1)=1.343085e+001; foe(n+1)=6.374090e+001; krok(n+1)=2.703354e-005; ng(n+1)=1.797801e+002;
n=291; farx(n+1)=1.328752e+001; foe(n+1)=6.365516e+001; krok(n+1)=5.519987e-006; ng(n+1)=3.635598e+002;
n=292; farx(n+1)=1.336503e+001; foe(n+1)=6.357067e+001; krok(n+1)=2.677366e-005; ng(n+1)=1.785986e+002;
n=293; farx(n+1)=1.322381e+001; foe(n+1)=6.348679e+001; krok(n+1)=5.530601e-006; ng(n+1)=3.596123e+002;
n=294; farx(n+1)=1.330086e+001; foe(n+1)=6.340419e+001; krok(n+1)=2.648762e-005; ng(n+1)=1.775872e+002;
n=295; farx(n+1)=1.316161e+001; foe(n+1)=6.332167e+001; krok(n+1)=5.523363e-006; ng(n+1)=3.568065e+002;
n=296; farx(n+1)=1.323770e+001; foe(n+1)=6.323928e+001; krok(n+1)=2.669060e-005; ng(n+1)=1.761740e+002;
n=297; farx(n+1)=1.309940e+001; foe(n+1)=6.315736e+001; krok(n+1)=5.519987e-006; ng(n+1)=3.558327e+002;
n=298; farx(n+1)=1.317514e+001; foe(n+1)=6.307631e+001; krok(n+1)=2.663711e-005; ng(n+1)=1.750558e+002;
n=299; farx(n+1)=1.303821e+001; foe(n+1)=6.299510e+001; krok(n+1)=5.505138e-006; ng(n+1)=3.545804e+002;
n=300; farx(n+1)=1.311310e+001; foe(n+1)=6.291394e+001; krok(n+1)=2.700051e-005; ng(n+1)=1.735837e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)