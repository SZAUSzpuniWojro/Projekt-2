%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.574860e+003; foe(n+1)=4.458767e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
%odnowa zmiennej metryki
n=1; farx(n+1)=3.771272e+003; foe(n+1)=3.770560e+003; krok(n+1)=5.823131e-004; ng(n+1)=3.739537e+003;
n=2; farx(n+1)=1.133623e+003; foe(n+1)=1.152335e+003; krok(n+1)=1.736904e-003; ng(n+1)=3.835953e+003;
n=3; farx(n+1)=1.351233e+003; foe(n+1)=9.083730e+002; krok(n+1)=1.131625e-004; ng(n+1)=8.986781e+003;
n=4; farx(n+1)=1.675164e+003; foe(n+1)=8.831728e+002; krok(n+1)=3.790794e-004; ng(n+1)=5.253772e+003;
n=5; farx(n+1)=1.428839e+003; foe(n+1)=7.776080e+002; krok(n+1)=3.050294e-003; ng(n+1)=1.607261e+003;
n=6; farx(n+1)=1.305437e+003; foe(n+1)=7.570272e+002; krok(n+1)=4.987505e-004; ng(n+1)=1.705972e+003;
n=7; farx(n+1)=9.730380e+002; foe(n+1)=6.895398e+002; krok(n+1)=6.647718e-004; ng(n+1)=5.777121e+003;
n=8; farx(n+1)=3.754978e+002; foe(n+1)=5.124832e+002; krok(n+1)=2.778799e-003; ng(n+1)=3.065752e+003;
n=9; farx(n+1)=3.411544e+002; foe(n+1)=4.928298e+002; krok(n+1)=2.199768e-005; ng(n+1)=6.582950e+003;
n=10; farx(n+1)=1.900567e+002; foe(n+1)=3.304836e+002; krok(n+1)=4.223924e-003; ng(n+1)=6.628305e+003;
n=11; farx(n+1)=1.886796e+002; foe(n+1)=3.271537e+002; krok(n+1)=1.954543e-005; ng(n+1)=4.896563e+003;
n=12; farx(n+1)=2.136403e+002; foe(n+1)=3.062542e+002; krok(n+1)=2.444561e-003; ng(n+1)=4.743640e+003;
n=13; farx(n+1)=1.157206e+002; foe(n+1)=2.775002e+002; krok(n+1)=3.159755e-003; ng(n+1)=3.793606e+003;
n=14; farx(n+1)=9.191980e+001; foe(n+1)=2.493027e+002; krok(n+1)=1.266922e-004; ng(n+1)=9.156168e+003;
n=15; farx(n+1)=8.097673e+001; foe(n+1)=2.382239e+002; krok(n+1)=2.601359e-004; ng(n+1)=2.955358e+003;
n=16; farx(n+1)=6.983290e+001; foe(n+1)=2.275583e+002; krok(n+1)=1.682049e-003; ng(n+1)=2.488558e+003;
n=17; farx(n+1)=6.798655e+001; foe(n+1)=2.188838e+002; krok(n+1)=5.418302e-004; ng(n+1)=2.719513e+003;
n=18; farx(n+1)=5.524298e+001; foe(n+1)=2.027306e+002; krok(n+1)=1.448033e-003; ng(n+1)=2.862623e+003;
n=19; farx(n+1)=4.495227e+001; foe(n+1)=1.868161e+002; krok(n+1)=5.499825e-004; ng(n+1)=7.048050e+003;
n=20; farx(n+1)=3.722624e+001; foe(n+1)=1.760519e+002; krok(n+1)=3.517334e-004; ng(n+1)=3.184996e+003;
n=21; farx(n+1)=3.292342e+001; foe(n+1)=1.638306e+002; krok(n+1)=1.358387e-003; ng(n+1)=5.068643e+003;
n=22; farx(n+1)=3.427956e+001; foe(n+1)=1.311165e+002; krok(n+1)=3.023106e-003; ng(n+1)=4.469994e+003;
n=23; farx(n+1)=3.764899e+001; foe(n+1)=1.178636e+002; krok(n+1)=3.192302e-004; ng(n+1)=1.079717e+004;
n=24; farx(n+1)=3.923885e+001; foe(n+1)=1.136259e+002; krok(n+1)=5.524716e-004; ng(n+1)=7.616195e+003;
n=25; farx(n+1)=3.871834e+001; foe(n+1)=1.047694e+002; krok(n+1)=3.032635e-003; ng(n+1)=3.297704e+003;
%odnowa zmiennej metryki
n=26; farx(n+1)=3.824733e+001; foe(n+1)=1.029575e+002; krok(n+1)=2.026389e-006; ng(n+1)=4.108286e+003;
n=27; farx(n+1)=2.969955e+001; foe(n+1)=8.089891e+001; krok(n+1)=8.876552e-005; ng(n+1)=2.623559e+003;
n=28; farx(n+1)=3.006514e+001; foe(n+1)=7.969208e+001; krok(n+1)=1.961746e-005; ng(n+1)=1.357569e+003;
n=29; farx(n+1)=3.168538e+001; foe(n+1)=7.235420e+001; krok(n+1)=1.864752e-004; ng(n+1)=1.127303e+003;
n=30; farx(n+1)=2.849759e+001; foe(n+1)=6.585818e+001; krok(n+1)=1.915637e-003; ng(n+1)=2.289421e+003;
n=31; farx(n+1)=2.535255e+001; foe(n+1)=6.363566e+001; krok(n+1)=2.054008e-004; ng(n+1)=2.380673e+003;
n=32; farx(n+1)=1.950722e+001; foe(n+1)=5.837784e+001; krok(n+1)=1.220985e-003; ng(n+1)=1.355778e+003;
n=33; farx(n+1)=1.621229e+001; foe(n+1)=4.511654e+001; krok(n+1)=4.312903e-003; ng(n+1)=3.243252e+003;
n=34; farx(n+1)=1.466156e+001; foe(n+1)=3.976332e+001; krok(n+1)=2.045232e-003; ng(n+1)=3.513267e+003;
n=35; farx(n+1)=1.282506e+001; foe(n+1)=3.646691e+001; krok(n+1)=1.965086e-003; ng(n+1)=1.627139e+003;
n=36; farx(n+1)=1.259274e+001; foe(n+1)=3.379459e+001; krok(n+1)=1.532691e-003; ng(n+1)=1.677309e+003;
n=37; farx(n+1)=8.782099e+000; foe(n+1)=2.879483e+001; krok(n+1)=2.426108e-002; ng(n+1)=7.579701e+002;
n=38; farx(n+1)=8.473928e+000; foe(n+1)=2.834303e+001; krok(n+1)=4.724072e-004; ng(n+1)=1.837870e+003;
n=39; farx(n+1)=7.829665e+000; foe(n+1)=2.742253e+001; krok(n+1)=3.683227e-003; ng(n+1)=1.242809e+003;
n=40; farx(n+1)=7.432424e+000; foe(n+1)=2.687224e+001; krok(n+1)=5.733303e-004; ng(n+1)=2.701399e+003;
n=41; farx(n+1)=6.624284e+000; foe(n+1)=2.536459e+001; krok(n+1)=2.804195e-003; ng(n+1)=2.052266e+003;
n=42; farx(n+1)=5.900196e+000; foe(n+1)=2.368164e+001; krok(n+1)=1.923079e-003; ng(n+1)=2.569214e+003;
n=43; farx(n+1)=5.092017e+000; foe(n+1)=2.050954e+001; krok(n+1)=8.017859e-003; ng(n+1)=1.219904e+003;
n=44; farx(n+1)=4.884715e+000; foe(n+1)=1.910007e+001; krok(n+1)=6.319511e-003; ng(n+1)=2.088875e+003;
n=45; farx(n+1)=4.738350e+000; foe(n+1)=1.767783e+001; krok(n+1)=4.743099e-003; ng(n+1)=8.422840e+002;
n=46; farx(n+1)=4.794305e+000; foe(n+1)=1.696210e+001; krok(n+1)=4.266621e-003; ng(n+1)=8.154813e+002;
n=47; farx(n+1)=4.976884e+000; foe(n+1)=1.611811e+001; krok(n+1)=4.352152e-003; ng(n+1)=6.601885e+002;
n=48; farx(n+1)=4.771086e+000; foe(n+1)=1.570907e+001; krok(n+1)=2.804195e-003; ng(n+1)=4.943820e+002;
n=49; farx(n+1)=4.610979e+000; foe(n+1)=1.524436e+001; krok(n+1)=1.129495e-002; ng(n+1)=8.483317e+002;
n=50; farx(n+1)=4.606488e+000; foe(n+1)=1.496057e+001; krok(n+1)=7.592571e-003; ng(n+1)=1.123411e+003;
%odnowa zmiennej metryki
n=51; farx(n+1)=4.527102e+000; foe(n+1)=1.473489e+001; krok(n+1)=9.006227e-006; ng(n+1)=6.889997e+002;
n=52; farx(n+1)=4.505238e+000; foe(n+1)=1.467060e+001; krok(n+1)=3.602457e-006; ng(n+1)=5.907367e+002;
n=53; farx(n+1)=4.551085e+000; foe(n+1)=1.453390e+001; krok(n+1)=8.066838e-005; ng(n+1)=1.937475e+002;
n=54; farx(n+1)=4.576907e+000; foe(n+1)=1.446435e+001; krok(n+1)=2.494331e-004; ng(n+1)=8.993675e+001;
n=55; farx(n+1)=4.511277e+000; foe(n+1)=1.431446e+001; krok(n+1)=3.505243e-004; ng(n+1)=1.144593e+002;
n=56; farx(n+1)=4.090842e+000; foe(n+1)=1.337705e+001; krok(n+1)=2.488900e-003; ng(n+1)=1.024056e+002;
n=57; farx(n+1)=3.842966e+000; foe(n+1)=1.260539e+001; krok(n+1)=3.223493e-003; ng(n+1)=3.594097e+002;
n=58; farx(n+1)=3.551238e+000; foe(n+1)=1.163651e+001; krok(n+1)=1.675019e-003; ng(n+1)=1.019501e+003;
n=59; farx(n+1)=3.417010e+000; foe(n+1)=1.057351e+001; krok(n+1)=4.904166e-003; ng(n+1)=4.696055e+002;
n=60; farx(n+1)=3.207650e+000; foe(n+1)=9.588267e+000; krok(n+1)=5.852382e-003; ng(n+1)=6.869585e+002;
n=61; farx(n+1)=3.065587e+000; foe(n+1)=8.916427e+000; krok(n+1)=1.700633e-003; ng(n+1)=1.420836e+003;
n=62; farx(n+1)=2.992436e+000; foe(n+1)=8.574229e+000; krok(n+1)=6.254538e-004; ng(n+1)=2.245142e+003;
n=63; farx(n+1)=3.070123e+000; foe(n+1)=8.165975e+000; krok(n+1)=4.679447e-003; ng(n+1)=1.227644e+003;
n=64; farx(n+1)=3.174452e+000; foe(n+1)=7.814534e+000; krok(n+1)=1.970836e-003; ng(n+1)=1.040963e+003;
n=65; farx(n+1)=3.004105e+000; foe(n+1)=7.377222e+000; krok(n+1)=1.321526e-002; ng(n+1)=9.367455e+002;
n=66; farx(n+1)=2.907786e+000; foe(n+1)=6.991292e+000; krok(n+1)=3.187345e-003; ng(n+1)=1.169425e+003;
n=67; farx(n+1)=2.898675e+000; foe(n+1)=6.858741e+000; krok(n+1)=1.935356e-003; ng(n+1)=7.324660e+002;
n=68; farx(n+1)=2.936499e+000; foe(n+1)=6.727576e+000; krok(n+1)=7.032104e-003; ng(n+1)=3.871213e+002;
n=69; farx(n+1)=2.904136e+000; foe(n+1)=6.368367e+000; krok(n+1)=1.525498e-002; ng(n+1)=7.147772e+002;
n=70; farx(n+1)=3.008810e+000; foe(n+1)=5.978755e+000; krok(n+1)=1.473291e-002; ng(n+1)=6.577605e+002;
n=71; farx(n+1)=2.928380e+000; foe(n+1)=5.865167e+000; krok(n+1)=1.818608e-002; ng(n+1)=2.855299e+002;
n=72; farx(n+1)=2.808968e+000; foe(n+1)=5.605577e+000; krok(n+1)=1.333847e-002; ng(n+1)=4.652920e+002;
n=73; farx(n+1)=2.643226e+000; foe(n+1)=5.264698e+000; krok(n+1)=8.899418e-003; ng(n+1)=1.085088e+003;
n=74; farx(n+1)=2.463701e+000; foe(n+1)=4.905440e+000; krok(n+1)=7.761346e-003; ng(n+1)=7.062582e+002;
n=75; farx(n+1)=2.217002e+000; foe(n+1)=4.559113e+000; krok(n+1)=1.417782e-002; ng(n+1)=1.070427e+003;
%odnowa zmiennej metryki
n=76; farx(n+1)=2.215147e+000; foe(n+1)=4.518901e+000; krok(n+1)=3.679859e-006; ng(n+1)=5.126630e+002;
n=77; farx(n+1)=2.203306e+000; foe(n+1)=4.437040e+000; krok(n+1)=5.751510e-006; ng(n+1)=4.515706e+002;
n=78; farx(n+1)=2.195603e+000; foe(n+1)=4.403283e+000; krok(n+1)=9.671965e-006; ng(n+1)=2.758970e+002;
n=79; farx(n+1)=2.182026e+000; foe(n+1)=4.377991e+000; krok(n+1)=1.553088e-004; ng(n+1)=6.763852e+001;
n=80; farx(n+1)=2.183365e+000; foe(n+1)=4.304236e+000; krok(n+1)=4.526501e-004; ng(n+1)=7.554600e+001;
n=81; farx(n+1)=2.182491e+000; foe(n+1)=4.255371e+000; krok(n+1)=6.935224e-004; ng(n+1)=5.501572e+001;
n=82; farx(n+1)=2.138935e+000; foe(n+1)=4.146910e+000; krok(n+1)=2.187783e-003; ng(n+1)=6.506330e+001;
n=83; farx(n+1)=2.080097e+000; foe(n+1)=4.076571e+000; krok(n+1)=9.532898e-004; ng(n+1)=1.125911e+002;
n=84; farx(n+1)=2.006301e+000; foe(n+1)=3.993655e+000; krok(n+1)=3.671029e-003; ng(n+1)=2.851180e+002;
n=85; farx(n+1)=1.874747e+000; foe(n+1)=3.835882e+000; krok(n+1)=4.323984e-003; ng(n+1)=2.921759e+002;
n=86; farx(n+1)=1.753601e+000; foe(n+1)=3.730489e+000; krok(n+1)=2.841023e-003; ng(n+1)=4.230005e+002;
n=87; farx(n+1)=1.490608e+000; foe(n+1)=3.531393e+000; krok(n+1)=6.390530e-003; ng(n+1)=6.887806e+002;
n=88; farx(n+1)=1.350740e+000; foe(n+1)=3.396880e+000; krok(n+1)=6.402952e-003; ng(n+1)=4.589675e+002;
n=89; farx(n+1)=1.264184e+000; foe(n+1)=3.293590e+000; krok(n+1)=7.637224e-003; ng(n+1)=2.682281e+002;
n=90; farx(n+1)=1.211791e+000; foe(n+1)=3.198314e+000; krok(n+1)=8.417037e-003; ng(n+1)=5.795911e+002;
n=91; farx(n+1)=1.153195e+000; foe(n+1)=3.120479e+000; krok(n+1)=5.999763e-003; ng(n+1)=1.133542e+003;
n=92; farx(n+1)=1.014732e+000; foe(n+1)=2.925882e+000; krok(n+1)=3.126305e-002; ng(n+1)=1.085218e+003;
n=93; farx(n+1)=9.562202e-001; foe(n+1)=2.760890e+000; krok(n+1)=1.479563e-002; ng(n+1)=8.637486e+002;
n=94; farx(n+1)=9.164076e-001; foe(n+1)=2.640480e+000; krok(n+1)=1.226153e-002; ng(n+1)=1.434238e+003;
n=95; farx(n+1)=8.917433e-001; foe(n+1)=2.549255e+000; krok(n+1)=9.045235e-003; ng(n+1)=7.306137e+002;
n=96; farx(n+1)=8.497453e-001; foe(n+1)=2.363418e+000; krok(n+1)=9.690117e-003; ng(n+1)=1.546105e+003;
n=97; farx(n+1)=7.871101e-001; foe(n+1)=2.132687e+000; krok(n+1)=1.370811e-002; ng(n+1)=1.735367e+003;
n=98; farx(n+1)=7.498098e-001; foe(n+1)=1.978315e+000; krok(n+1)=4.795168e-003; ng(n+1)=7.783200e+002;
n=99; farx(n+1)=7.392505e-001; foe(n+1)=1.895734e+000; krok(n+1)=8.996953e-003; ng(n+1)=1.146593e+003;
n=100; farx(n+1)=7.343523e-001; foe(n+1)=1.849325e+000; krok(n+1)=6.192369e-003; ng(n+1)=1.388626e+003;
%odnowa zmiennej metryki
n=101; farx(n+1)=7.342053e-001; foe(n+1)=1.841961e+000; krok(n+1)=3.511715e-007; ng(n+1)=4.890537e+002;
n=102; farx(n+1)=7.343480e-001; foe(n+1)=1.834782e+000; krok(n+1)=2.828453e-006; ng(n+1)=2.473903e+002;
n=103; farx(n+1)=7.349065e-001; foe(n+1)=1.828686e+000; krok(n+1)=5.032633e-006; ng(n+1)=1.705176e+002;
n=104; farx(n+1)=7.335968e-001; foe(n+1)=1.769844e+000; krok(n+1)=2.431917e-004; ng(n+1)=9.334399e+001;
n=105; farx(n+1)=7.316129e-001; foe(n+1)=1.760606e+000; krok(n+1)=8.904993e-005; ng(n+1)=5.814595e+001;
n=106; farx(n+1)=7.295730e-001; foe(n+1)=1.743716e+000; krok(n+1)=9.776439e-005; ng(n+1)=7.118185e+001;
n=107; farx(n+1)=7.249642e-001; foe(n+1)=1.722083e+000; krok(n+1)=4.168272e-004; ng(n+1)=4.803461e+001;
n=108; farx(n+1)=7.226121e-001; foe(n+1)=1.693835e+000; krok(n+1)=7.010405e-004; ng(n+1)=9.261363e+001;
n=109; farx(n+1)=7.054623e-001; foe(n+1)=1.657586e+000; krok(n+1)=1.387045e-003; ng(n+1)=1.467417e+002;
n=110; farx(n+1)=6.778375e-001; foe(n+1)=1.596571e+000; krok(n+1)=4.743099e-003; ng(n+1)=1.646402e+002;
n=111; farx(n+1)=6.510403e-001; foe(n+1)=1.522451e+000; krok(n+1)=8.073226e-003; ng(n+1)=1.641541e+002;
n=112; farx(n+1)=6.316557e-001; foe(n+1)=1.472123e+000; krok(n+1)=6.970444e-003; ng(n+1)=7.889574e+002;
n=113; farx(n+1)=6.257555e-001; foe(n+1)=1.448789e+000; krok(n+1)=4.008930e-003; ng(n+1)=5.020816e+002;
n=114; farx(n+1)=6.271316e-001; foe(n+1)=1.424968e+000; krok(n+1)=7.366454e-003; ng(n+1)=4.086794e+002;
n=115; farx(n+1)=6.263708e-001; foe(n+1)=1.371948e+000; krok(n+1)=9.928643e-003; ng(n+1)=5.186282e+002;
n=116; farx(n+1)=6.134648e-001; foe(n+1)=1.349075e+000; krok(n+1)=1.285062e-002; ng(n+1)=6.217602e+002;
n=117; farx(n+1)=6.019419e-001; foe(n+1)=1.331529e+000; krok(n+1)=1.778267e-002; ng(n+1)=3.338327e+002;
n=118; farx(n+1)=5.929668e-001; foe(n+1)=1.300739e+000; krok(n+1)=1.468412e-002; ng(n+1)=5.818082e+002;
n=119; farx(n+1)=5.919131e-001; foe(n+1)=1.291372e+000; krok(n+1)=5.027940e-003; ng(n+1)=3.014542e+002;
n=120; farx(n+1)=5.861760e-001; foe(n+1)=1.272236e+000; krok(n+1)=1.670157e-002; ng(n+1)=7.402753e+002;
n=121; farx(n+1)=5.866264e-001; foe(n+1)=1.254811e+000; krok(n+1)=1.664870e-002; ng(n+1)=4.422908e+002;
n=122; farx(n+1)=5.872442e-001; foe(n+1)=1.236082e+000; krok(n+1)=1.809047e-002; ng(n+1)=1.468714e+002;
n=123; farx(n+1)=5.929927e-001; foe(n+1)=1.206671e+000; krok(n+1)=1.274938e-002; ng(n+1)=2.930400e+002;
n=124; farx(n+1)=5.921709e-001; foe(n+1)=1.190373e+000; krok(n+1)=1.301252e-002; ng(n+1)=5.260610e+002;
n=125; farx(n+1)=5.973699e-001; foe(n+1)=1.161173e+000; krok(n+1)=1.378623e-002; ng(n+1)=5.991688e+002;
%odnowa zmiennej metryki
n=126; farx(n+1)=5.973463e-001; foe(n+1)=1.160246e+000; krok(n+1)=2.851733e-007; ng(n+1)=2.991536e+002;
n=127; farx(n+1)=5.971336e-001; foe(n+1)=1.155406e+000; krok(n+1)=1.957485e-006; ng(n+1)=2.007336e+002;
n=128; farx(n+1)=5.967227e-001; foe(n+1)=1.152140e+000; krok(n+1)=3.073770e-006; ng(n+1)=1.543004e+002;
n=129; farx(n+1)=5.958908e-001; foe(n+1)=1.147618e+000; krok(n+1)=3.440034e-005; ng(n+1)=6.320556e+001;
n=130; farx(n+1)=5.948393e-001; foe(n+1)=1.141830e+000; krok(n+1)=6.667236e-005; ng(n+1)=4.922052e+001;
n=131; farx(n+1)=5.938697e-001; foe(n+1)=1.138199e+000; krok(n+1)=1.733806e-004; ng(n+1)=2.728994e+001;
n=132; farx(n+1)=5.940087e-001; foe(n+1)=1.131789e+000; krok(n+1)=5.157570e-004; ng(n+1)=1.739518e+001;
n=133; farx(n+1)=5.972951e-001; foe(n+1)=1.118618e+000; krok(n+1)=1.499367e-003; ng(n+1)=2.485603e+001;
n=134; farx(n+1)=6.000063e-001; foe(n+1)=1.110852e+000; krok(n+1)=6.710424e-004; ng(n+1)=1.549179e+002;
n=135; farx(n+1)=6.053890e-001; foe(n+1)=1.103829e+000; krok(n+1)=6.130766e-003; ng(n+1)=2.737944e+002;
n=136; farx(n+1)=6.113522e-001; foe(n+1)=1.084797e+000; krok(n+1)=9.454744e-003; ng(n+1)=7.138748e+002;
n=137; farx(n+1)=6.114638e-001; foe(n+1)=1.080499e+000; krok(n+1)=1.358387e-003; ng(n+1)=4.166607e+002;
n=138; farx(n+1)=6.129522e-001; foe(n+1)=1.076491e+000; krok(n+1)=4.036613e-003; ng(n+1)=1.787716e+002;
n=139; farx(n+1)=6.139300e-001; foe(n+1)=1.057603e+000; krok(n+1)=1.336240e-002; ng(n+1)=2.473134e+002;
n=140; farx(n+1)=6.114042e-001; foe(n+1)=1.045353e+000; krok(n+1)=1.797739e-002; ng(n+1)=1.453264e+002;
n=141; farx(n+1)=6.164664e-001; foe(n+1)=1.018858e+000; krok(n+1)=2.643052e-002; ng(n+1)=1.132006e+002;
n=142; farx(n+1)=6.148961e-001; foe(n+1)=1.011036e+000; krok(n+1)=6.131265e-003; ng(n+1)=8.648041e+002;
n=143; farx(n+1)=6.143881e-001; foe(n+1)=9.946757e-001; krok(n+1)=7.821513e-003; ng(n+1)=3.927226e+002;
n=144; farx(n+1)=6.156648e-001; foe(n+1)=9.837306e-001; krok(n+1)=4.624718e-003; ng(n+1)=4.505312e+002;
n=145; farx(n+1)=6.153545e-001; foe(n+1)=9.625250e-001; krok(n+1)=2.386882e-002; ng(n+1)=4.606364e+002;
n=146; farx(n+1)=6.110307e-001; foe(n+1)=9.226800e-001; krok(n+1)=4.852216e-002; ng(n+1)=1.978921e+002;
n=147; farx(n+1)=6.082921e-001; foe(n+1)=8.967431e-001; krok(n+1)=3.196395e-002; ng(n+1)=7.140616e+002;
n=148; farx(n+1)=6.068970e-001; foe(n+1)=8.850565e-001; krok(n+1)=7.763326e-003; ng(n+1)=2.950757e+002;
n=149; farx(n+1)=6.020522e-001; foe(n+1)=8.586019e-001; krok(n+1)=9.636858e-003; ng(n+1)=7.794374e+002;
n=150; farx(n+1)=5.974744e-001; foe(n+1)=8.354276e-001; krok(n+1)=1.495979e-002; ng(n+1)=4.884728e+002;
%odnowa zmiennej metryki
n=151; farx(n+1)=5.974934e-001; foe(n+1)=8.277291e-001; krok(n+1)=8.508023e-008; ng(n+1)=1.221335e+003;
n=152; farx(n+1)=5.975251e-001; foe(n+1)=8.253345e-001; krok(n+1)=1.948696e-006; ng(n+1)=1.741605e+002;
n=153; farx(n+1)=5.974928e-001; foe(n+1)=8.202756e-001; krok(n+1)=5.918863e-006; ng(n+1)=1.442939e+002;
n=154; farx(n+1)=5.965024e-001; foe(n+1)=8.081989e-001; krok(n+1)=2.428432e-005; ng(n+1)=1.169384e+002;
n=155; farx(n+1)=5.944470e-001; foe(n+1)=7.890893e-001; krok(n+1)=3.782669e-005; ng(n+1)=1.303300e+002;
n=156; farx(n+1)=5.940817e-001; foe(n+1)=7.697939e-001; krok(n+1)=4.187132e-004; ng(n+1)=4.552862e+001;
n=157; farx(n+1)=5.940153e-001; foe(n+1)=7.657113e-001; krok(n+1)=8.254679e-005; ng(n+1)=5.165062e+001;
n=158; farx(n+1)=5.935456e-001; foe(n+1)=7.470655e-001; krok(n+1)=1.548346e-003; ng(n+1)=3.871634e+001;
n=159; farx(n+1)=5.927659e-001; foe(n+1)=7.456882e-001; krok(n+1)=6.696507e-004; ng(n+1)=2.407515e+002;
n=160; farx(n+1)=5.897402e-001; foe(n+1)=7.403114e-001; krok(n+1)=2.883948e-003; ng(n+1)=2.819826e+002;
n=161; farx(n+1)=5.903417e-001; foe(n+1)=7.303395e-001; krok(n+1)=7.837498e-003; ng(n+1)=6.025664e+002;
n=162; farx(n+1)=5.905532e-001; foe(n+1)=7.234350e-001; krok(n+1)=2.850058e-003; ng(n+1)=8.696949e+001;
n=163; farx(n+1)=5.905580e-001; foe(n+1)=7.163793e-001; krok(n+1)=6.485153e-003; ng(n+1)=3.933034e+002;
n=164; farx(n+1)=5.882176e-001; foe(n+1)=7.104689e-001; krok(n+1)=1.703305e-002; ng(n+1)=5.625516e+002;
n=165; farx(n+1)=5.821896e-001; foe(n+1)=6.943826e-001; krok(n+1)=2.876644e-002; ng(n+1)=2.710854e+002;
n=166; farx(n+1)=5.830134e-001; foe(n+1)=6.830787e-001; krok(n+1)=8.092065e-003; ng(n+1)=1.125023e+002;
n=167; farx(n+1)=5.816619e-001; foe(n+1)=6.735546e-001; krok(n+1)=1.552039e-002; ng(n+1)=2.491068e+002;
n=168; farx(n+1)=5.784053e-001; foe(n+1)=6.657260e-001; krok(n+1)=1.080437e-002; ng(n+1)=3.464123e+002;
n=169; farx(n+1)=5.671535e-001; foe(n+1)=6.441369e-001; krok(n+1)=1.376894e-002; ng(n+1)=5.863931e+002;
n=170; farx(n+1)=5.634633e-001; foe(n+1)=6.355256e-001; krok(n+1)=3.096185e-003; ng(n+1)=1.010346e+003;
n=171; farx(n+1)=5.616401e-001; foe(n+1)=6.194871e-001; krok(n+1)=6.853004e-003; ng(n+1)=1.015215e+003;
n=172; farx(n+1)=5.584676e-001; foe(n+1)=6.064101e-001; krok(n+1)=2.753897e-002; ng(n+1)=3.151423e+002;
n=173; farx(n+1)=5.575944e-001; foe(n+1)=6.031082e-001; krok(n+1)=9.778245e-003; ng(n+1)=3.017760e+002;
n=174; farx(n+1)=5.461983e-001; foe(n+1)=5.843298e-001; krok(n+1)=2.513448e-002; ng(n+1)=7.176853e+002;
n=175; farx(n+1)=5.371787e-001; foe(n+1)=5.659557e-001; krok(n+1)=2.991958e-002; ng(n+1)=2.082296e+002;
%odnowa zmiennej metryki
n=176; farx(n+1)=5.371984e-001; foe(n+1)=5.628194e-001; krok(n+1)=8.395551e-008; ng(n+1)=7.516347e+002;
n=177; farx(n+1)=5.372229e-001; foe(n+1)=5.602883e-001; krok(n+1)=2.471149e-006; ng(n+1)=1.754622e+002;
n=178; farx(n+1)=5.371880e-001; foe(n+1)=5.477061e-001; krok(n+1)=7.622483e-006; ng(n+1)=2.106639e+002;
n=179; farx(n+1)=5.368368e-001; foe(n+1)=5.322417e-001; krok(n+1)=1.085468e-005; ng(n+1)=1.852172e+002;
n=180; farx(n+1)=5.368296e-001; foe(n+1)=5.314508e-001; krok(n+1)=1.839286e-005; ng(n+1)=3.864847e+001;
n=181; farx(n+1)=5.367869e-001; foe(n+1)=5.268146e-001; krok(n+1)=2.325873e-004; ng(n+1)=2.654892e+001;
n=182; farx(n+1)=5.370750e-001; foe(n+1)=5.252252e-001; krok(n+1)=1.888305e-004; ng(n+1)=1.824275e+001;
n=183; farx(n+1)=5.366493e-001; foe(n+1)=5.237877e-001; krok(n+1)=6.856719e-004; ng(n+1)=2.102338e+001;
n=184; farx(n+1)=5.351822e-001; foe(n+1)=5.202223e-001; krok(n+1)=2.426783e-003; ng(n+1)=2.144152e+001;
n=185; farx(n+1)=5.323886e-001; foe(n+1)=5.158285e-001; krok(n+1)=2.329278e-003; ng(n+1)=2.207341e+001;
n=186; farx(n+1)=5.293651e-001; foe(n+1)=5.059299e-001; krok(n+1)=9.612071e-003; ng(n+1)=1.941163e+002;
n=187; farx(n+1)=5.307213e-001; foe(n+1)=5.024268e-001; krok(n+1)=4.150488e-003; ng(n+1)=8.698686e+001;
n=188; farx(n+1)=5.302917e-001; foe(n+1)=4.997468e-001; krok(n+1)=4.508154e-003; ng(n+1)=1.177367e+002;
n=189; farx(n+1)=5.309246e-001; foe(n+1)=4.962175e-001; krok(n+1)=4.738412e-003; ng(n+1)=2.944044e+002;
n=190; farx(n+1)=5.332775e-001; foe(n+1)=4.802187e-001; krok(n+1)=2.946581e-002; ng(n+1)=3.722972e+002;
n=191; farx(n+1)=5.307696e-001; foe(n+1)=4.746399e-001; krok(n+1)=1.263902e-002; ng(n+1)=2.694239e+002;
n=192; farx(n+1)=5.269074e-001; foe(n+1)=4.629671e-001; krok(n+1)=2.768352e-002; ng(n+1)=2.690501e+002;
n=193; farx(n+1)=5.239751e-001; foe(n+1)=4.560789e-001; krok(n+1)=7.251821e-003; ng(n+1)=9.219417e+002;
n=194; farx(n+1)=5.199285e-001; foe(n+1)=4.461423e-001; krok(n+1)=3.065020e-002; ng(n+1)=5.651834e+002;
n=195; farx(n+1)=5.209183e-001; foe(n+1)=4.415581e-001; krok(n+1)=1.193571e-002; ng(n+1)=6.054322e+002;
n=196; farx(n+1)=5.211418e-001; foe(n+1)=4.372614e-001; krok(n+1)=1.204622e-002; ng(n+1)=4.453672e+002;
n=197; farx(n+1)=5.203044e-001; foe(n+1)=4.318658e-001; krok(n+1)=1.901220e-002; ng(n+1)=1.805283e+002;
n=198; farx(n+1)=5.160648e-001; foe(n+1)=4.190826e-001; krok(n+1)=7.581460e-002; ng(n+1)=4.813583e+002;
n=199; farx(n+1)=5.125080e-001; foe(n+1)=4.117344e-001; krok(n+1)=1.650423e-002; ng(n+1)=4.267451e+002;
n=200; farx(n+1)=5.100690e-001; foe(n+1)=4.047069e-001; krok(n+1)=1.302239e-002; ng(n+1)=3.353512e+002;
%odnowa zmiennej metryki
n=201; farx(n+1)=5.100750e-001; foe(n+1)=4.042941e-001; krok(n+1)=1.484034e-007; ng(n+1)=2.748318e+002;
n=202; farx(n+1)=5.100626e-001; foe(n+1)=4.023686e-001; krok(n+1)=1.762952e-006; ng(n+1)=1.527773e+002;
n=203; farx(n+1)=5.100643e-001; foe(n+1)=4.003033e-001; krok(n+1)=2.272267e-006; ng(n+1)=1.531507e+002;
n=204; farx(n+1)=5.100144e-001; foe(n+1)=3.992446e-001; krok(n+1)=1.185899e-005; ng(n+1)=4.891368e+001;
n=205; farx(n+1)=5.099541e-001; foe(n+1)=3.980641e-001; krok(n+1)=7.250092e-005; ng(n+1)=2.207814e+001;
n=206; farx(n+1)=5.098586e-001; foe(n+1)=3.955336e-001; krok(n+1)=4.129769e-004; ng(n+1)=1.513292e+001;
n=207; farx(n+1)=5.099100e-001; foe(n+1)=3.941505e-001; krok(n+1)=2.266194e-004; ng(n+1)=1.582598e+001;
n=208; farx(n+1)=5.085712e-001; foe(n+1)=3.911289e-001; krok(n+1)=1.708199e-003; ng(n+1)=1.490348e+001;
n=209; farx(n+1)=5.080311e-001; foe(n+1)=3.892572e-001; krok(n+1)=6.767654e-004; ng(n+1)=4.907351e+001;
n=210; farx(n+1)=5.083833e-001; foe(n+1)=3.877646e-001; krok(n+1)=4.375565e-003; ng(n+1)=1.337830e+002;
n=211; farx(n+1)=5.091151e-001; foe(n+1)=3.851585e-001; krok(n+1)=3.657312e-003; ng(n+1)=2.458814e+002;
n=212; farx(n+1)=5.079073e-001; foe(n+1)=3.826937e-001; krok(n+1)=3.255599e-003; ng(n+1)=3.362542e+002;
n=213; farx(n+1)=5.060631e-001; foe(n+1)=3.780236e-001; krok(n+1)=1.462925e-002; ng(n+1)=5.077074e+002;
n=214; farx(n+1)=5.066861e-001; foe(n+1)=3.744457e-001; krok(n+1)=2.223039e-002; ng(n+1)=6.576622e+001;
n=215; farx(n+1)=5.029273e-001; foe(n+1)=3.701896e-001; krok(n+1)=1.694661e-002; ng(n+1)=2.979899e+002;
n=216; farx(n+1)=4.986563e-001; foe(n+1)=3.637448e-001; krok(n+1)=1.321323e-002; ng(n+1)=2.069562e+002;
n=217; farx(n+1)=4.971982e-001; foe(n+1)=3.611134e-001; krok(n+1)=8.471180e-003; ng(n+1)=1.389457e+002;
n=218; farx(n+1)=4.958748e-001; foe(n+1)=3.591790e-001; krok(n+1)=1.450364e-002; ng(n+1)=1.863941e+002;
n=219; farx(n+1)=4.934629e-001; foe(n+1)=3.562785e-001; krok(n+1)=1.290032e-002; ng(n+1)=1.786308e+002;
n=220; farx(n+1)=4.912725e-001; foe(n+1)=3.539464e-001; krok(n+1)=1.086709e-002; ng(n+1)=2.151480e+002;
n=221; farx(n+1)=4.896667e-001; foe(n+1)=3.499493e-001; krok(n+1)=6.601690e-002; ng(n+1)=2.849451e+002;
n=222; farx(n+1)=4.892070e-001; foe(n+1)=3.468929e-001; krok(n+1)=1.825465e-002; ng(n+1)=1.462289e+002;
n=223; farx(n+1)=4.893280e-001; foe(n+1)=3.430934e-001; krok(n+1)=2.835564e-002; ng(n+1)=1.761440e+002;
n=224; farx(n+1)=4.906681e-001; foe(n+1)=3.383049e-001; krok(n+1)=5.625683e-002; ng(n+1)=2.045754e+002;
n=225; farx(n+1)=4.917228e-001; foe(n+1)=3.351364e-001; krok(n+1)=1.729712e-002; ng(n+1)=3.072099e+002;
%odnowa zmiennej metryki
n=226; farx(n+1)=4.917184e-001; foe(n+1)=3.336805e-001; krok(n+1)=2.224567e-007; ng(n+1)=3.543833e+002;
n=227; farx(n+1)=4.917052e-001; foe(n+1)=3.318209e-001; krok(n+1)=4.536859e-007; ng(n+1)=2.798934e+002;
n=228; farx(n+1)=4.917449e-001; foe(n+1)=3.306061e-001; krok(n+1)=1.099167e-005; ng(n+1)=5.421406e+001;
n=229; farx(n+1)=4.917349e-001; foe(n+1)=3.303259e-001; krok(n+1)=8.318815e-006; ng(n+1)=3.248776e+001;
n=230; farx(n+1)=4.916910e-001; foe(n+1)=3.285941e-001; krok(n+1)=1.222501e-004; ng(n+1)=1.861285e+001;
n=231; farx(n+1)=4.917930e-001; foe(n+1)=3.278509e-001; krok(n+1)=1.418485e-004; ng(n+1)=1.117536e+001;
n=232; farx(n+1)=4.920125e-001; foe(n+1)=3.271895e-001; krok(n+1)=2.008920e-004; ng(n+1)=1.016742e+001;
n=233; farx(n+1)=4.925374e-001; foe(n+1)=3.265556e-001; krok(n+1)=8.160960e-004; ng(n+1)=8.207493e+000;
n=234; farx(n+1)=4.926339e-001; foe(n+1)=3.261732e-001; krok(n+1)=1.369522e-003; ng(n+1)=1.375596e+001;
n=235; farx(n+1)=4.923402e-001; foe(n+1)=3.256723e-001; krok(n+1)=4.539152e-003; ng(n+1)=1.769912e+001;
n=236; farx(n+1)=4.919845e-001; foe(n+1)=3.253030e-001; krok(n+1)=2.293768e-003; ng(n+1)=2.189820e+001;
n=237; farx(n+1)=4.922760e-001; foe(n+1)=3.236634e-001; krok(n+1)=1.389524e-002; ng(n+1)=1.640899e+001;
n=238; farx(n+1)=4.901651e-001; foe(n+1)=3.216679e-001; krok(n+1)=1.767909e-002; ng(n+1)=1.686039e+002;
n=239; farx(n+1)=4.890335e-001; foe(n+1)=3.196899e-001; krok(n+1)=2.334578e-002; ng(n+1)=5.322454e+001;
n=240; farx(n+1)=4.882463e-001; foe(n+1)=3.183937e-001; krok(n+1)=8.161186e-003; ng(n+1)=3.132280e+002;
n=241; farx(n+1)=4.890175e-001; foe(n+1)=3.176797e-001; krok(n+1)=1.660195e-002; ng(n+1)=1.700645e+002;
n=242; farx(n+1)=4.887651e-001; foe(n+1)=3.167912e-001; krok(n+1)=1.391254e-002; ng(n+1)=1.895507e+002;
n=243; farx(n+1)=4.883446e-001; foe(n+1)=3.156665e-001; krok(n+1)=1.776017e-002; ng(n+1)=1.973834e+002;
n=244; farx(n+1)=4.891977e-001; foe(n+1)=3.139092e-001; krok(n+1)=2.243330e-002; ng(n+1)=3.913770e+002;
n=245; farx(n+1)=4.892809e-001; foe(n+1)=3.112593e-001; krok(n+1)=4.674985e-002; ng(n+1)=3.035160e+002;
n=246; farx(n+1)=4.902363e-001; foe(n+1)=3.077902e-001; krok(n+1)=6.918848e-002; ng(n+1)=2.661132e+002;
n=247; farx(n+1)=4.896897e-001; foe(n+1)=3.051458e-001; krok(n+1)=3.771383e-002; ng(n+1)=3.243662e+002;
n=248; farx(n+1)=4.890614e-001; foe(n+1)=3.010394e-001; krok(n+1)=4.134280e-002; ng(n+1)=1.428729e+002;
n=249; farx(n+1)=4.893698e-001; foe(n+1)=2.999999e-001; krok(n+1)=2.378798e-002; ng(n+1)=1.696482e+002;
n=250; farx(n+1)=4.907791e-001; foe(n+1)=2.963626e-001; krok(n+1)=3.353295e-002; ng(n+1)=2.991752e+002;
%odnowa zmiennej metryki
n=251; farx(n+1)=4.907798e-001; foe(n+1)=2.962949e-001; krok(n+1)=2.668615e-007; ng(n+1)=6.844160e+001;
n=252; farx(n+1)=4.907762e-001; foe(n+1)=2.961543e-001; krok(n+1)=5.950526e-007; ng(n+1)=7.790724e+001;
n=253; farx(n+1)=4.907626e-001; foe(n+1)=2.958086e-001; krok(n+1)=4.182313e-006; ng(n+1)=4.755824e+001;
n=254; farx(n+1)=4.907621e-001; foe(n+1)=2.956480e-001; krok(n+1)=9.471564e-006; ng(n+1)=2.311712e+001;
n=255; farx(n+1)=4.907721e-001; foe(n+1)=2.954148e-001; krok(n+1)=6.836820e-005; ng(n+1)=9.880082e+000;
n=256; farx(n+1)=4.908319e-001; foe(n+1)=2.946518e-001; krok(n+1)=5.823194e-004; ng(n+1)=6.139051e+000;
n=257; farx(n+1)=4.908608e-001; foe(n+1)=2.940902e-001; krok(n+1)=2.967298e-004; ng(n+1)=8.297344e+000;
n=258; farx(n+1)=4.902881e-001; foe(n+1)=2.926746e-001; krok(n+1)=1.604967e-003; ng(n+1)=6.894042e+000;
n=259; farx(n+1)=4.899797e-001; foe(n+1)=2.923851e-001; krok(n+1)=4.616173e-004; ng(n+1)=1.605874e+001;
n=260; farx(n+1)=4.891775e-001; foe(n+1)=2.920261e-001; krok(n+1)=3.427028e-003; ng(n+1)=1.881734e+001;
n=261; farx(n+1)=4.883442e-001; foe(n+1)=2.913659e-001; krok(n+1)=5.557598e-003; ng(n+1)=2.462578e+001;
n=262; farx(n+1)=4.882208e-001; foe(n+1)=2.908068e-001; krok(n+1)=2.436016e-003; ng(n+1)=3.644462e+001;
n=263; farx(n+1)=4.889610e-001; foe(n+1)=2.890378e-001; krok(n+1)=1.987953e-002; ng(n+1)=3.954862e+001;
n=264; farx(n+1)=4.885300e-001; foe(n+1)=2.878796e-001; krok(n+1)=1.007731e-002; ng(n+1)=1.051607e+002;
n=265; farx(n+1)=4.884774e-001; foe(n+1)=2.863161e-001; krok(n+1)=1.540867e-002; ng(n+1)=1.491318e+002;
n=266; farx(n+1)=4.877603e-001; foe(n+1)=2.850515e-001; krok(n+1)=2.378377e-002; ng(n+1)=2.482017e+002;
n=267; farx(n+1)=4.866161e-001; foe(n+1)=2.828960e-001; krok(n+1)=2.592196e-002; ng(n+1)=1.785853e+002;
n=268; farx(n+1)=4.862113e-001; foe(n+1)=2.813712e-001; krok(n+1)=1.352965e-002; ng(n+1)=1.349943e+002;
n=269; farx(n+1)=4.865847e-001; foe(n+1)=2.805984e-001; krok(n+1)=1.381076e-002; ng(n+1)=1.938203e+002;
n=270; farx(n+1)=4.865930e-001; foe(n+1)=2.797970e-001; krok(n+1)=2.121282e-002; ng(n+1)=9.414985e+001;
n=271; farx(n+1)=4.868149e-001; foe(n+1)=2.790376e-001; krok(n+1)=2.779047e-002; ng(n+1)=1.276160e+002;
n=272; farx(n+1)=4.871251e-001; foe(n+1)=2.778174e-001; krok(n+1)=5.857529e-002; ng(n+1)=9.957151e+001;
n=273; farx(n+1)=4.874537e-001; foe(n+1)=2.754717e-001; krok(n+1)=1.247927e-001; ng(n+1)=7.412333e+001;
n=274; farx(n+1)=4.875621e-001; foe(n+1)=2.744874e-001; krok(n+1)=1.706813e-002; ng(n+1)=2.272639e+002;
n=275; farx(n+1)=4.880790e-001; foe(n+1)=2.706637e-001; krok(n+1)=6.905442e-002; ng(n+1)=3.284782e+002;
%odnowa zmiennej metryki
n=276; farx(n+1)=4.880804e-001; foe(n+1)=2.702510e-001; krok(n+1)=1.365385e-007; ng(n+1)=1.882012e+002;
n=277; farx(n+1)=4.880800e-001; foe(n+1)=2.700188e-001; krok(n+1)=8.143365e-007; ng(n+1)=8.470797e+001;
n=278; farx(n+1)=4.880796e-001; foe(n+1)=2.697732e-001; krok(n+1)=3.333783e-006; ng(n+1)=4.938470e+001;
n=279; farx(n+1)=4.881052e-001; foe(n+1)=2.696612e-001; krok(n+1)=3.686149e-005; ng(n+1)=1.028154e+001;
n=280; farx(n+1)=4.881608e-001; foe(n+1)=2.695220e-001; krok(n+1)=3.400479e-005; ng(n+1)=1.118793e+001;
n=281; farx(n+1)=4.886897e-001; foe(n+1)=2.688367e-001; krok(n+1)=4.532388e-004; ng(n+1)=6.940065e+000;
n=282; farx(n+1)=4.888820e-001; foe(n+1)=2.686269e-001; krok(n+1)=3.149159e-004; ng(n+1)=5.094659e+000;
n=283; farx(n+1)=4.891120e-001; foe(n+1)=2.683389e-001; krok(n+1)=5.593002e-004; ng(n+1)=4.412839e+000;
n=284; farx(n+1)=4.893034e-001; foe(n+1)=2.680793e-001; krok(n+1)=7.432430e-004; ng(n+1)=5.772027e+000;
n=285; farx(n+1)=4.894096e-001; foe(n+1)=2.676158e-001; krok(n+1)=2.728971e-003; ng(n+1)=4.423672e+000;
n=286; farx(n+1)=4.892797e-001; foe(n+1)=2.673287e-001; krok(n+1)=2.707062e-003; ng(n+1)=6.278970e+000;
n=287; farx(n+1)=4.895825e-001; foe(n+1)=2.661009e-001; krok(n+1)=9.683310e-003; ng(n+1)=6.332981e+000;
n=288; farx(n+1)=4.897274e-001; foe(n+1)=2.652244e-001; krok(n+1)=9.030079e-003; ng(n+1)=2.664918e+001;
n=289; farx(n+1)=4.886411e-001; foe(n+1)=2.645888e-001; krok(n+1)=1.238677e-002; ng(n+1)=7.586333e+001;
n=290; farx(n+1)=4.890058e-001; foe(n+1)=2.639703e-001; krok(n+1)=1.932256e-002; ng(n+1)=1.049305e+002;
n=291; farx(n+1)=4.897429e-001; foe(n+1)=2.632723e-001; krok(n+1)=5.793922e-002; ng(n+1)=1.329410e+002;
n=292; farx(n+1)=4.899431e-001; foe(n+1)=2.620517e-001; krok(n+1)=2.024365e-002; ng(n+1)=1.369320e+002;
n=293; farx(n+1)=4.902028e-001; foe(n+1)=2.600123e-001; krok(n+1)=3.209254e-002; ng(n+1)=1.977262e+002;
n=294; farx(n+1)=4.902277e-001; foe(n+1)=2.594354e-001; krok(n+1)=1.169933e-002; ng(n+1)=1.777123e+002;
n=295; farx(n+1)=4.903469e-001; foe(n+1)=2.583860e-001; krok(n+1)=5.857529e-002; ng(n+1)=3.911479e+001;
n=296; farx(n+1)=4.905824e-001; foe(n+1)=2.564205e-001; krok(n+1)=7.487116e-002; ng(n+1)=1.756537e+002;
n=297; farx(n+1)=4.907609e-001; foe(n+1)=2.557115e-001; krok(n+1)=1.556427e-002; ng(n+1)=6.730169e+001;
n=298; farx(n+1)=4.909812e-001; foe(n+1)=2.551933e-001; krok(n+1)=1.820180e-002; ng(n+1)=9.083472e+001;
n=299; farx(n+1)=4.906474e-001; foe(n+1)=2.543457e-001; krok(n+1)=2.741202e-002; ng(n+1)=1.959345e+002;
n=300; farx(n+1)=4.913230e-001; foe(n+1)=2.525340e-001; krok(n+1)=5.234463e-002; ng(n+1)=2.721112e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
