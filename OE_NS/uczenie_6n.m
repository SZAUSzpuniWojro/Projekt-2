%uczenie predyktora oe
clear all;
n=0; farx(n+1)=4.406945e+003; foe(n+1)=4.557609e+003; wspucz(n+1)=0.000000e+000; ng(n+1)=0.000000e+000;
n=1; farx(n+1)=3.758812e+003; foe(n+1)=3.785570e+003; krok(n+1)=5.534736e-004; ng(n+1)=3.361259e+003;
n=2; farx(n+1)=1.128715e+003; foe(n+1)=1.291581e+003; krok(n+1)=2.041615e-003; ng(n+1)=2.655039e+003;
n=3; farx(n+1)=5.350142e+002; foe(n+1)=8.520528e+002; krok(n+1)=1.319295e-004; ng(n+1)=7.274537e+003;
n=4; farx(n+1)=5.597741e+002; foe(n+1)=8.200504e+002; krok(n+1)=2.630324e-004; ng(n+1)=1.617871e+003;
n=5; farx(n+1)=4.813698e+002; foe(n+1)=7.973255e+002; krok(n+1)=1.154043e-004; ng(n+1)=2.148220e+003;
n=6; farx(n+1)=5.023098e+002; foe(n+1)=7.799785e+002; krok(n+1)=1.614953e-004; ng(n+1)=1.272051e+003;
n=7; farx(n+1)=4.415656e+002; foe(n+1)=7.622684e+002; krok(n+1)=1.107642e-004; ng(n+1)=1.698639e+003;
n=8; farx(n+1)=4.603699e+002; foe(n+1)=7.461360e+002; krok(n+1)=1.520709e-004; ng(n+1)=1.381319e+003;
n=9; farx(n+1)=4.052721e+002; foe(n+1)=7.290599e+002; krok(n+1)=1.046887e-004; ng(n+1)=1.705678e+003;
n=10; farx(n+1)=4.209976e+002; foe(n+1)=7.131755e+002; krok(n+1)=1.457243e-004; ng(n+1)=1.405360e+003;
n=11; farx(n+1)=3.704763e+002; foe(n+1)=6.962540e+002; krok(n+1)=9.833971e-005; ng(n+1)=1.752209e+003;
n=12; farx(n+1)=3.837529e+002; foe(n+1)=6.804281e+002; krok(n+1)=1.389822e-004; ng(n+1)=1.440104e+003;
n=13; farx(n+1)=3.375172e+002; foe(n+1)=6.635554e+002; krok(n+1)=9.243013e-005; ng(n+1)=1.807617e+003;
n=14; farx(n+1)=3.488621e+002; foe(n+1)=6.477823e+002; krok(n+1)=1.316270e-004; ng(n+1)=1.484291e+003;
n=15; farx(n+1)=3.072772e+002; foe(n+1)=6.309996e+002; krok(n+1)=8.559513e-005; ng(n+1)=1.863888e+003;
n=16; farx(n+1)=3.166118e+002; foe(n+1)=6.150004e+002; krok(n+1)=1.307587e-004; ng(n+1)=1.486600e+003;
n=17; farx(n+1)=2.782027e+002; foe(n+1)=5.980661e+002; krok(n+1)=7.961977e-005; ng(n+1)=1.960146e+003;
n=18; farx(n+1)=2.863369e+002; foe(n+1)=5.824019e+002; krok(n+1)=1.205665e-004; ng(n+1)=1.537962e+003;
n=19; farx(n+1)=2.520762e+002; foe(n+1)=5.659297e+002; krok(n+1)=7.492449e-005; ng(n+1)=1.995252e+003;
n=20; farx(n+1)=2.591119e+002; foe(n+1)=5.505803e+002; krok(n+1)=1.133097e-004; ng(n+1)=1.570048e+003;
n=21; farx(n+1)=2.284269e+002; foe(n+1)=5.345429e+002; krok(n+1)=7.005353e-005; ng(n+1)=2.046369e+003;
n=22; farx(n+1)=2.345731e+002; foe(n+1)=5.196060e+002; krok(n+1)=1.067624e-004; ng(n+1)=1.596709e+003;
n=23; farx(n+1)=2.070617e+002; foe(n+1)=5.040649e+002; krok(n+1)=6.543012e-005; ng(n+1)=2.099848e+003;
n=24; farx(n+1)=2.125073e+002; foe(n+1)=4.896612e+002; krok(n+1)=1.001340e-004; ng(n+1)=1.623498e+003;
n=25; farx(n+1)=1.879360e+002; foe(n+1)=4.747115e+002; krok(n+1)=6.097986e-005; ng(n+1)=2.151272e+003;
n=26; farx(n+1)=1.927219e+002; foe(n+1)=4.608555e+002; krok(n+1)=9.441524e-005; ng(n+1)=1.641214e+003;
n=27; farx(n+1)=1.706318e+002; foe(n+1)=4.465401e+002; krok(n+1)=5.720244e-005; ng(n+1)=2.202738e+003;
n=28; farx(n+1)=1.749002e+002; foe(n+1)=4.334332e+002; krok(n+1)=8.559513e-005; ng(n+1)=1.679949e+003;
n=29; farx(n+1)=1.553703e+002; foe(n+1)=4.200603e+002; krok(n+1)=5.440277e-005; ng(n+1)=2.208226e+003;
n=30; farx(n+1)=1.592892e+002; foe(n+1)=4.076885e+002; krok(n+1)=7.910691e-005; ng(n+1)=1.718447e+003;
n=31; farx(n+1)=1.418493e+002; foe(n+1)=3.950696e+002; krok(n+1)=5.113377e-005; ng(n+1)=2.233215e+003;
n=32; farx(n+1)=1.454076e+002; foe(n+1)=3.833549e+002; krok(n+1)=7.447980e-005; ng(n+1)=1.742077e+003;
n=33; farx(n+1)=1.297699e+002; foe(n+1)=3.714046e+002; krok(n+1)=4.766085e-005; ng(n+1)=2.266888e+003;
n=34; farx(n+1)=1.329163e+002; foe(n+1)=3.602837e+002; krok(n+1)=7.066590e-005; ng(n+1)=1.747028e+003;
n=35; farx(n+1)=1.187358e+002; foe(n+1)=3.490107e+002; krok(n+1)=4.506169e-005; ng(n+1)=2.298077e+003;
n=36; farx(n+1)=1.217071e+002; foe(n+1)=3.386642e+002; krok(n+1)=6.485137e-005; ng(n+1)=1.790513e+003;
n=37; farx(n+1)=1.090955e+002; foe(n+1)=3.281519e+002; krok(n+1)=4.229784e-005; ng(n+1)=2.303982e+003;
n=38; farx(n+1)=1.117102e+002; foe(n+1)=3.183414e+002; krok(n+1)=6.244891e-005; ng(n+1)=1.778663e+003;
n=39; farx(n+1)=1.002307e+002; foe(n+1)=3.084331e+002; krok(n+1)=3.980988e-005; ng(n+1)=2.331161e+003;
n=40; farx(n+1)=1.026476e+002; foe(n+1)=2.992976e+002; krok(n+1)=5.805985e-005; ng(n+1)=1.800949e+003;
n=41; farx(n+1)=9.235867e+001; foe(n+1)=2.900951e+002; krok(n+1)=3.762088e-005; ng(n+1)=2.326179e+003;
n=42; farx(n+1)=9.450646e+001; foe(n+1)=2.815246e+002; krok(n+1)=5.493831e-005; ng(n+1)=1.795172e+003;
n=43; farx(n+1)=8.519271e+001; foe(n+1)=2.729790e+002; krok(n+1)=3.583315e-005; ng(n+1)=2.318884e+003;
n=44; farx(n+1)=8.723645e+001; foe(n+1)=2.650315e+002; krok(n+1)=5.210339e-005; ng(n+1)=1.807519e+003;
n=45; farx(n+1)=7.877322e+001; foe(n+1)=2.570377e+002; krok(n+1)=3.369455e-005; ng(n+1)=2.324860e+003;
n=46; farx(n+1)=8.060640e+001; foe(n+1)=2.496154e+002; krok(n+1)=4.960813e-005; ng(n+1)=1.792179e+003;
n=47; farx(n+1)=7.295638e+001; foe(n+1)=2.422127e+002; krok(n+1)=3.195675e-005; ng(n+1)=2.308558e+003;
n=48; farx(n+1)=7.462987e+001; foe(n+1)=2.353004e+002; krok(n+1)=4.765565e-005; ng(n+1)=1.776093e+003;
n=49; farx(n+1)=6.766811e+001; foe(n+1)=2.284186e+002; krok(n+1)=3.023762e-005; ng(n+1)=2.298452e+003;
n=50; farx(n+1)=6.919394e+001; foe(n+1)=2.219978e+002; krok(n+1)=4.565415e-005; ng(n+1)=1.758404e+003;
n=51; farx(n+1)=6.284463e+001; foe(n+1)=2.156314e+002; krok(n+1)=2.885108e-005; ng(n+1)=2.277010e+003;
n=52; farx(n+1)=6.424120e+001; foe(n+1)=2.097103e+002; krok(n+1)=4.279756e-005; ng(n+1)=1.748099e+003;
n=53; farx(n+1)=5.849284e+001; foe(n+1)=2.038987e+002; krok(n+1)=2.791016e-005; ng(n+1)=2.226719e+003;
n=54; farx(n+1)=5.981865e+001; foe(n+1)=1.984536e+002; krok(n+1)=4.073849e-005; ng(n+1)=1.738595e+003;
n=55; farx(n+1)=5.456977e+001; foe(n+1)=1.930901e+002; krok(n+1)=2.668806e-005; ng(n+1)=2.196218e+003;
n=56; farx(n+1)=5.580706e+001; foe(n+1)=1.880633e+002; krok(n+1)=3.923630e-005; ng(n+1)=1.714381e+003;
n=57; farx(n+1)=5.100911e+001; foe(n+1)=1.831076e+002; krok(n+1)=2.545461e-005; ng(n+1)=2.166517e+003;
n=58; farx(n+1)=5.214391e+001; foe(n+1)=1.784588e+002; krok(n+1)=3.782669e-005; ng(n+1)=1.681028e+003;
n=59; farx(n+1)=4.774835e+001; foe(n+1)=1.738992e+002; krok(n+1)=2.449700e-005; ng(n+1)=2.127200e+003;
n=60; farx(n+1)=4.881330e+001; foe(n+1)=1.696269e+002; krok(n+1)=3.613934e-005; ng(n+1)=1.656597e+003;
n=61; farx(n+1)=4.480206e+001; foe(n+1)=1.654428e+002; krok(n+1)=2.357491e-005; ng(n+1)=2.081116e+003;
n=62; farx(n+1)=4.579207e+001; foe(n+1)=1.615004e+002; krok(n+1)=3.501276e-005; ng(n+1)=1.620719e+003;
n=63; farx(n+1)=4.209274e+001; foe(n+1)=1.576453e+002; krok(n+1)=2.277170e-005; ng(n+1)=2.041576e+003;
n=64; farx(n+1)=4.303208e+001; foe(n+1)=1.540341e+002; krok(n+1)=3.333618e-005; ng(n+1)=1.597174e+003;
n=65; farx(n+1)=3.965535e+001; foe(n+1)=1.505047e+002; krok(n+1)=2.201176e-005; ng(n+1)=1.988307e+003;
n=66; farx(n+1)=4.052166e+001; foe(n+1)=1.471708e+002; krok(n+1)=3.223481e-005; ng(n+1)=1.556432e+003;
n=67; farx(n+1)=3.737642e+001; foe(n+1)=1.439326e+002; krok(n+1)=2.167257e-005; ng(n+1)=1.938273e+003;
n=68; farx(n+1)=3.823611e+001; foe(n+1)=1.409127e+002; krok(n+1)=2.992557e-005; ng(n+1)=1.556593e+003;
n=69; farx(n+1)=3.536486e+001; foe(n+1)=1.379617e+002; krok(n+1)=2.122479e-005; ng(n+1)=1.871727e+003;
n=70; farx(n+1)=3.617041e+001; foe(n+1)=1.351781e+002; krok(n+1)=2.890345e-005; ng(n+1)=1.519659e+003;
n=71; farx(n+1)=3.352829e+001; foe(n+1)=1.324647e+002; krok(n+1)=2.060293e-005; ng(n+1)=1.821701e+003;
n=72; farx(n+1)=3.428808e+001; foe(n+1)=1.298940e+002; krok(n+1)=2.829063e-005; ng(n+1)=1.478744e+003;
n=73; farx(n+1)=3.185279e+001; foe(n+1)=1.273840e+002; krok(n+1)=1.977673e-005; ng(n+1)=1.784310e+003;
n=74; farx(n+1)=3.255357e+001; foe(n+1)=1.249964e+002; krok(n+1)=2.800232e-005; ng(n+1)=1.425736e+003;
n=75; farx(n+1)=3.029027e+001; foe(n+1)=1.226718e+002; krok(n+1)=1.912864e-005; ng(n+1)=1.749124e+003;
n=76; farx(n+1)=3.095198e+001; foe(n+1)=1.204690e+002; krok(n+1)=2.726938e-005; ng(n+1)=1.388992e+003;
n=77; farx(n+1)=2.885491e+001; foe(n+1)=1.183239e+002; krok(n+1)=1.860264e-005; ng(n+1)=1.705025e+003;
n=78; farx(n+1)=2.947620e+001; foe(n+1)=1.162900e+002; krok(n+1)=2.647223e-005; ng(n+1)=1.352223e+003;
n=79; farx(n+1)=2.753343e+001; foe(n+1)=1.143150e+002; krok(n+1)=1.819728e-005; ng(n+1)=1.655382e+003;
n=80; farx(n+1)=2.812497e+001; foe(n+1)=1.124383e+002; krok(n+1)=2.580710e-005; ng(n+1)=1.318449e+003;
n=81; farx(n+1)=2.631922e+001; foe(n+1)=1.106120e+002; krok(n+1)=1.773106e-005; ng(n+1)=1.612972e+003;
n=82; farx(n+1)=2.687536e+001; foe(n+1)=1.088764e+002; krok(n+1)=2.509886e-005; ng(n+1)=1.283157e+003;
n=83; farx(n+1)=2.519431e+001; foe(n+1)=1.071922e+002; krok(n+1)=1.743790e-005; ng(n+1)=1.564168e+003;
n=84; farx(n+1)=2.573047e+001; foe(n+1)=1.055918e+002; krok(n+1)=2.444110e-005; ng(n+1)=1.255508e+003;
n=85; farx(n+1)=2.416690e+001; foe(n+1)=1.040307e+002; krok(n+1)=1.697768e-005; ng(n+1)=1.523784e+003;
n=86; farx(n+1)=2.467369e+001; foe(n+1)=1.025442e+002; krok(n+1)=2.418894e-005; ng(n+1)=1.218409e+003;
n=87; farx(n+1)=2.321760e+001; foe(n+1)=1.010914e+002; krok(n+1)=1.643712e-005; ng(n+1)=1.489618e+003;
n=88; farx(n+1)=2.369356e+001; foe(n+1)=9.970508e+001; krok(n+1)=2.418894e-005; ng(n+1)=1.177477e+003;
n=89; farx(n+1)=2.233568e+001; foe(n+1)=9.834765e+001; krok(n+1)=1.584198e-005; ng(n+1)=1.464017e+003;
n=90; farx(n+1)=2.277008e+001; foe(n+1)=9.704942e+001; krok(n+1)=2.394830e-005; ng(n+1)=1.135248e+003;
n=91; farx(n+1)=2.149814e+001; foe(n+1)=9.579365e+001; krok(n+1)=1.565988e-005; ng(n+1)=1.419441e+003;
n=92; farx(n+1)=2.191644e+001; foe(n+1)=9.459607e+001; krok(n+1)=2.310648e-005; ng(n+1)=1.110804e+003;
n=93; farx(n+1)=2.073599e+001; foe(n+1)=9.343562e+001; krok(n+1)=1.539716e-005; ng(n+1)=1.371524e+003;
n=94; farx(n+1)=2.113431e+001; foe(n+1)=9.231818e+001; krok(n+1)=2.307378e-005; ng(n+1)=1.078809e+003;
n=95; farx(n+1)=2.001956e+001; foe(n+1)=9.123038e+001; krok(n+1)=1.502405e-005; ng(n+1)=1.347366e+003;
n=96; farx(n+1)=2.039759e+001; foe(n+1)=9.019150e+001; krok(n+1)=2.254192e-005; ng(n+1)=1.051567e+003;
n=97; farx(n+1)=1.935401e+001; foe(n+1)=8.918216e+001; krok(n+1)=1.484032e-005; ng(n+1)=1.304940e+003;
n=98; farx(n+1)=1.971751e+001; foe(n+1)=8.821580e+001; krok(n+1)=2.219138e-005; ng(n+1)=1.026654e+003;
n=99; farx(n+1)=1.873787e+001; foe(n+1)=8.727277e+001; krok(n+1)=1.453671e-005; ng(n+1)=1.272661e+003;
n=100; farx(n+1)=1.907718e+001; foe(n+1)=8.636888e+001; krok(n+1)=2.170937e-005; ng(n+1)=9.989802e+002;
n=101; farx(n+1)=1.815772e+001; foe(n+1)=8.549476e+001; krok(n+1)=1.448879e-005; ng(n+1)=1.226953e+003;
n=102; farx(n+1)=1.848656e+001; foe(n+1)=8.465501e+001; krok(n+1)=2.114892e-005; ng(n+1)=9.782538e+002;
n=103; farx(n+1)=1.761605e+001; foe(n+1)=8.383945e+001; krok(n+1)=1.442554e-005; ng(n+1)=1.190454e+003;
n=104; farx(n+1)=1.793637e+001; foe(n+1)=8.305971e+001; krok(n+1)=2.047517e-005; ng(n+1)=9.600883e+002;
n=105; farx(n+1)=1.712560e+001; foe(n+1)=8.229907e+001; krok(n+1)=1.416371e-005; ng(n+1)=1.152080e+003;
n=106; farx(n+1)=1.742795e+001; foe(n+1)=8.156210e+001; krok(n+1)=2.077412e-005; ng(n+1)=9.298256e+002;
n=107; farx(n+1)=1.666110e+001; foe(n+1)=8.084233e+001; krok(n+1)=1.373458e-005; ng(n+1)=1.135502e+003;
n=108; farx(n+1)=1.694673e+001; foe(n+1)=8.014532e+001; krok(n+1)=2.102242e-005; ng(n+1)=9.010366e+002;
n=109; farx(n+1)=1.621988e+001; foe(n+1)=7.946414e+001; krok(n+1)=1.338683e-005; ng(n+1)=1.116967e+003;
n=110; farx(n+1)=1.648907e+001; foe(n+1)=7.880609e+001; krok(n+1)=2.097007e-005; ng(n+1)=8.749696e+002;
n=111; farx(n+1)=1.579787e+001; foe(n+1)=7.816577e+001; krok(n+1)=1.328171e-005; ng(n+1)=1.088565e+003;
n=112; farx(n+1)=1.605891e+001; foe(n+1)=7.755130e+001; krok(n+1)=2.044242e-005; ng(n+1)=8.568136e+002;
n=113; farx(n+1)=1.540811e+001; foe(n+1)=7.695159e+001; krok(n+1)=1.315377e-005; ng(n+1)=1.054720e+003;
n=114; farx(n+1)=1.565773e+001; foe(n+1)=7.637223e+001; krok(n+1)=2.031935e-005; ng(n+1)=8.345683e+002;
n=115; farx(n+1)=1.504394e+001; foe(n+1)=7.580666e+001; krok(n+1)=1.290553e-005; ng(n+1)=1.029103e+003;
n=116; farx(n+1)=1.528162e+001; foe(n+1)=7.525641e+001; krok(n+1)=2.058283e-005; ng(n+1)=8.097059e+002;
n=117; farx(n+1)=1.469553e+001; foe(n+1)=7.471828e+001; krok(n+1)=1.264803e-005; ng(n+1)=1.013231e+003;
n=118; farx(n+1)=1.492280e+001; foe(n+1)=7.419820e+001; krok(n+1)=2.047517e-005; ng(n+1)=7.885795e+002;
n=119; farx(n+1)=1.436559e+001; foe(n+1)=7.368990e+001; krok(n+1)=1.251675e-005; ng(n+1)=9.883335e+002;
n=120; farx(n+1)=1.458450e+001; foe(n+1)=7.319931e+001; krok(n+1)=2.026093e-005; ng(n+1)=7.704364e+002;
n=121; farx(n+1)=1.405733e+001; foe(n+1)=7.271923e+001; krok(n+1)=1.234279e-005; ng(n+1)=9.638161e+002;
n=122; farx(n+1)=1.426643e+001; foe(n+1)=7.225275e+001; krok(n+1)=2.036925e-005; ng(n+1)=7.502442e+002;
n=123; farx(n+1)=1.376113e+001; foe(n+1)=7.179632e+001; krok(n+1)=1.222055e-005; ng(n+1)=9.468050e+002;
n=124; farx(n+1)=1.396373e+001; foe(n+1)=7.135693e+001; krok(n+1)=1.999909e-005; ng(n+1)=7.351172e+002;
n=125; farx(n+1)=1.348424e+001; foe(n+1)=7.092628e+001; krok(n+1)=1.214216e-005; ng(n+1)=9.203208e+002;
n=126; farx(n+1)=1.367996e+001; foe(n+1)=7.050984e+001; krok(n+1)=1.990494e-005; ng(n+1)=7.185021e+002;
n=127; farx(n+1)=1.321942e+001; foe(n+1)=7.010069e+001; krok(n+1)=1.209447e-005; ng(n+1)=9.012617e+002;
n=128; farx(n+1)=1.340892e+001; foe(n+1)=6.970873e+001; krok(n+1)=1.928098e-005; ng(n+1)=7.059773e+002;
n=129; farx(n+1)=1.297814e+001; foe(n+1)=6.932459e+001; krok(n+1)=1.197415e-005; ng(n+1)=8.683327e+002;
n=130; farx(n+1)=1.315935e+001; foe(n+1)=6.894707e+001; krok(n+1)=1.990494e-005; ng(n+1)=6.848589e+002;
n=131; farx(n+1)=1.274231e+001; foe(n+1)=6.857576e+001; krok(n+1)=1.171380e-005; ng(n+1)=8.666832e+002;
n=132; farx(n+1)=1.291725e+001; foe(n+1)=6.821617e+001; krok(n+1)=1.990494e-005; ng(n+1)=6.690302e+002;
n=133; farx(n+1)=1.251852e+001; foe(n+1)=6.786187e+001; krok(n+1)=1.157184e-005; ng(n+1)=8.505424e+002;
n=134; farx(n+1)=1.268671e+001; foe(n+1)=6.751857e+001; krok(n+1)=1.990494e-005; ng(n+1)=6.532265e+002;
n=135; farx(n+1)=1.230460e+001; foe(n+1)=6.718073e+001; krok(n+1)=1.148266e-005; ng(n+1)=8.333990e+002;
n=136; farx(n+1)=1.246650e+001; foe(n+1)=6.685411e+001; krok(n+1)=1.964039e-005; ng(n+1)=6.394904e+002;
n=137; farx(n+1)=1.210288e+001; foe(n+1)=6.653391e+001; krok(n+1)=1.143043e-005; ng(n+1)=8.102151e+002;
n=138; farx(n+1)=1.225873e+001; foe(n+1)=6.622193e+001; krok(n+1)=1.954543e-005; ng(n+1)=6.252723e+002;
n=139; farx(n+1)=1.190935e+001; foe(n+1)=6.591687e+001; krok(n+1)=1.141354e-005; ng(n+1)=7.917128e+002;
n=140; farx(n+1)=1.206151e+001; foe(n+1)=6.562095e+001; krok(n+1)=1.920477e-005; ng(n+1)=6.142635e+002;
n=141; farx(n+1)=1.172657e+001; foe(n+1)=6.533079e+001; krok(n+1)=1.141354e-005; ng(n+1)=7.713451e+002;
n=142; farx(n+1)=1.187485e+001; foe(n+1)=6.504937e+001; krok(n+1)=1.891271e-005; ng(n+1)=6.032635e+002;
n=143; farx(n+1)=1.155466e+001; foe(n+1)=6.477296e+001; krok(n+1)=1.136531e-005; ng(n+1)=7.521089e+002;
n=144; farx(n+1)=1.169852e+001; foe(n+1)=6.450344e+001; krok(n+1)=1.890200e-005; ng(n+1)=5.906135e+002;
n=145; farx(n+1)=1.139364e+001; foe(n+1)=6.423839e+001; krok(n+1)=1.115447e-005; ng(n+1)=7.388100e+002;
n=146; farx(n+1)=1.153192e+001; foe(n+1)=6.397604e+001; krok(n+1)=1.961815e-005; ng(n+1)=5.749169e+002;
n=147; farx(n+1)=1.123669e+001; foe(n+1)=6.371745e+001; krok(n+1)=1.085468e-005; ng(n+1)=7.397348e+002;
n=148; farx(n+1)=1.136942e+001; foe(n+1)=6.346235e+001; krok(n+1)=2.018691e-005; ng(n+1)=5.606807e+002;
n=149; farx(n+1)=1.108199e+001; foe(n+1)=6.321108e+001; krok(n+1)=1.071362e-005; ng(n+1)=7.361319e+002;
n=150; farx(n+1)=1.121040e+001; foe(n+1)=6.296712e+001; krok(n+1)=1.995189e-005; ng(n+1)=5.500637e+002;
n=151; farx(n+1)=1.093173e+001; foe(n+1)=6.272776e+001; krok(n+1)=1.083629e-005; ng(n+1)=7.172721e+002;
n=152; farx(n+1)=1.105834e+001; foe(n+1)=6.249901e+001; krok(n+1)=1.892601e-005; ng(n+1)=5.439746e+002;
n=153; farx(n+1)=1.079651e+001; foe(n+1)=6.227420e+001; krok(n+1)=1.085468e-005; ng(n+1)=6.869772e+002;
n=154; farx(n+1)=1.091883e+001; foe(n+1)=6.205222e+001; krok(n+1)=1.931405e-005; ng(n+1)=5.317570e+002;
n=155; farx(n+1)=1.066190e+001; foe(n+1)=6.183408e+001; krok(n+1)=1.083629e-005; ng(n+1)=6.818758e+002;
n=156; farx(n+1)=1.078204e+001; foe(n+1)=6.162402e+001; krok(n+1)=1.859953e-005; ng(n+1)=5.250093e+002;
n=157; farx(n+1)=1.053814e+001; foe(n+1)=6.141737e+001; krok(n+1)=1.085468e-005; ng(n+1)=6.581430e+002;
n=158; farx(n+1)=1.065488e+001; foe(n+1)=6.121433e+001; krok(n+1)=1.877358e-005; ng(n+1)=5.145363e+002;
n=159; farx(n+1)=1.041681e+001; foe(n+1)=6.101440e+001; krok(n+1)=1.083629e-005; ng(n+1)=6.504522e+002;
n=160; farx(n+1)=1.053103e+001; foe(n+1)=6.082088e+001; krok(n+1)=1.827909e-005; ng(n+1)=5.073456e+002;
n=161; farx(n+1)=1.030147e+001; foe(n+1)=6.063035e+001; krok(n+1)=1.095376e-005; ng(n+1)=6.311249e+002;
n=162; farx(n+1)=1.041342e+001; foe(n+1)=6.044588e+001; krok(n+1)=1.773994e-005; ng(n+1)=5.008479e+002;
n=163; farx(n+1)=1.019320e+001; foe(n+1)=6.026423e+001; krok(n+1)=1.103022e-005; ng(n+1)=6.116354e+002;
n=164; farx(n+1)=1.030239e+001; foe(n+1)=6.008695e+001; krok(n+1)=1.751338e-005; ng(n+1)=4.930730e+002;
n=165; farx(n+1)=1.008963e+001; foe(n+1)=5.991250e+001; krok(n+1)=1.106017e-005; ng(n+1)=5.978486e+002;
n=166; farx(n+1)=1.019660e+001; foe(n+1)=5.974205e+001; krok(n+1)=1.733471e-005; ng(n+1)=4.856488e+002;
n=167; farx(n+1)=9.990940e+000; foe(n+1)=5.957402e+001; krok(n+1)=1.104492e-005; ng(n+1)=5.863997e+002;
n=168; farx(n+1)=1.009560e+001; foe(n+1)=5.940956e+001; krok(n+1)=1.729605e-005; ng(n+1)=4.777678e+002;
n=169; farx(n+1)=9.896301e+000; foe(n+1)=5.924713e+001; krok(n+1)=1.099167e-005; ng(n+1)=5.777383e+002;
n=170; farx(n+1)=9.998526e+000; foe(n+1)=5.908817e+001; krok(n+1)=1.729605e-005; ng(n+1)=4.698101e+002;
n=171; farx(n+1)=9.806612e+000; foe(n+1)=5.893103e+001; krok(n+1)=1.085468e-005; ng(n+1)=5.695954e+002;
n=172; farx(n+1)=9.905457e+000; foe(n+1)=5.877588e+001; krok(n+1)=1.763713e-005; ng(n+1)=4.599899e+002;
n=173; farx(n+1)=9.716893e+000; foe(n+1)=5.862296e+001; krok(n+1)=1.083629e-005; ng(n+1)=5.658699e+002;
n=174; farx(n+1)=9.813844e+000; foe(n+1)=5.847429e+001; krok(n+1)=1.718695e-005; ng(n+1)=4.544924e+002;
n=175; farx(n+1)=9.633015e+000; foe(n+1)=5.832773e+001; krok(n+1)=1.085468e-005; ng(n+1)=5.501980e+002;
n=176; farx(n+1)=9.727465e+000; foe(n+1)=5.818347e+001; krok(n+1)=1.726474e-005; ng(n+1)=4.466236e+002;
n=177; farx(n+1)=9.553227e+000; foe(n+1)=5.804140e+001; krok(n+1)=1.069939e-005; ng(n+1)=5.434077e+002;
n=178; farx(n+1)=9.644968e+000; foe(n+1)=5.790002e+001; krok(n+1)=1.777520e-005; ng(n+1)=4.370797e+002;
n=179; farx(n+1)=9.474098e+000; foe(n+1)=5.776084e+001; krok(n+1)=1.055114e-005; ng(n+1)=5.442776e+002;
n=180; farx(n+1)=9.564111e+000; foe(n+1)=5.762368e+001; krok(n+1)=1.791657e-005; ng(n+1)=4.297629e+002;
n=181; farx(n+1)=9.397498e+000; foe(n+1)=5.748815e+001; krok(n+1)=1.046329e-005; ng(n+1)=5.403633e+002;
n=182; farx(n+1)=9.485238e+000; foe(n+1)=5.735510e+001; krok(n+1)=1.792006e-005; ng(n+1)=4.227103e+002;
n=183; farx(n+1)=9.324051e+000; foe(n+1)=5.722386e+001; krok(n+1)=1.039083e-005; ng(n+1)=5.325710e+002;
n=184; farx(n+1)=9.409764e+000; foe(n+1)=5.709419e+001; krok(n+1)=1.816250e-005; ng(n+1)=4.151076e+002;
n=185; farx(n+1)=9.252345e+000; foe(n+1)=5.696618e+001; krok(n+1)=1.030147e-005; ng(n+1)=5.298265e+002;
n=186; farx(n+1)=9.336170e+000; foe(n+1)=5.684031e+001; krok(n+1)=1.821554e-005; ng(n+1)=4.083821e+002;
n=187; farx(n+1)=9.182421e+000; foe(n+1)=5.671600e+001; krok(n+1)=1.029142e-005; ng(n+1)=5.239480e+002;
n=188; farx(n+1)=9.264642e+000; foe(n+1)=5.659458e+001; krok(n+1)=1.799457e-005; ng(n+1)=4.029131e+002;
n=189; farx(n+1)=9.116455e+000; foe(n+1)=5.647459e+001; krok(n+1)=1.023989e-005; ng(n+1)=5.141136e+002;
n=190; farx(n+1)=9.196643e+000; foe(n+1)=5.635590e+001; krok(n+1)=1.827909e-005; ng(n+1)=3.955588e+002;
n=191; farx(n+1)=9.051510e+000; foe(n+1)=5.623871e+001; krok(n+1)=1.017255e-005; ng(n+1)=5.116630e+002;
n=192; farx(n+1)=9.130045e+000; foe(n+1)=5.612368e+001; krok(n+1)=1.821997e-005; ng(n+1)=3.897210e+002;
n=193; farx(n+1)=8.988846e+000; foe(n+1)=5.601008e+001; krok(n+1)=1.015759e-005; ng(n+1)=5.044266e+002;
n=194; farx(n+1)=9.065805e+000; foe(n+1)=5.589858e+001; krok(n+1)=1.815977e-005; ng(n+1)=3.840551e+002;
n=195; farx(n+1)=8.928259e+000; foe(n+1)=5.578841e+001; krok(n+1)=1.015216e-005; ng(n+1)=4.973975e+002;
n=196; farx(n+1)=9.003744e+000; foe(n+1)=5.568043e+001; krok(n+1)=1.804864e-005; ng(n+1)=3.787441e+002;
n=197; farx(n+1)=8.869873e+000; foe(n+1)=5.557367e+001; krok(n+1)=1.015216e-005; ng(n+1)=4.898702e+002;
n=198; farx(n+1)=8.943760e+000; foe(n+1)=5.546901e+001; krok(n+1)=1.792006e-005; ng(n+1)=3.735112e+002;
n=199; farx(n+1)=8.813599e+000; foe(n+1)=5.536566e+001; krok(n+1)=1.015968e-005; ng(n+1)=4.817411e+002;
n=200; farx(n+1)=8.886031e+000; foe(n+1)=5.526414e+001; krok(n+1)=1.782892e-005; ng(n+1)=3.683686e+002;
n=201; farx(n+1)=8.759251e+000; foe(n+1)=5.516389e+001; krok(n+1)=1.015968e-005; ng(n+1)=4.746517e+002;
n=202; farx(n+1)=8.830303e+000; foe(n+1)=5.506540e+001; krok(n+1)=1.773994e-005; ng(n+1)=3.633766e+002;
n=203; farx(n+1)=8.706769e+000; foe(n+1)=5.496811e+001; krok(n+1)=1.015759e-005; ng(n+1)=4.678418e+002;
n=204; farx(n+1)=8.776692e+000; foe(n+1)=5.487249e+001; krok(n+1)=1.771309e-005; ng(n+1)=3.584721e+002;
n=205; farx(n+1)=8.656349e+000; foe(n+1)=5.477779e+001; krok(n+1)=1.010485e-005; ng(n+1)=4.627570e+002;
n=206; farx(n+1)=8.724822e+000; foe(n+1)=5.468446e+001; krok(n+1)=1.782680e-005; ng(n+1)=3.529480e+002;
n=207; farx(n+1)=8.606776e+000; foe(n+1)=5.459209e+001; krok(n+1)=1.009346e-005; ng(n+1)=4.587928e+002;
n=208; farx(n+1)=8.673942e+000; foe(n+1)=5.450170e+001; krok(n+1)=1.760998e-005; ng(n+1)=3.486838e+002;
n=209; farx(n+1)=8.559057e+000; foe(n+1)=5.441230e+001; krok(n+1)=1.013046e-005; ng(n+1)=4.503453e+002;
n=210; farx(n+1)=8.624957e+000; foe(n+1)=5.432463e+001; krok(n+1)=1.744002e-005; ng(n+1)=3.444002e+002;
n=211; farx(n+1)=8.513073e+000; foe(n+1)=5.423797e+001; krok(n+1)=1.015216e-005; ng(n+1)=4.427646e+002;
n=212; farx(n+1)=8.577891e+000; foe(n+1)=5.415279e+001; krok(n+1)=1.736589e-005; ng(n+1)=3.400425e+002;
n=213; farx(n+1)=8.468696e+000; foe(n+1)=5.406847e+001; krok(n+1)=1.013760e-005; ng(n+1)=4.371912e+002;
n=214; farx(n+1)=8.532269e+000; foe(n+1)=5.398552e+001; krok(n+1)=1.733471e-005; ng(n+1)=3.354869e+002;
n=215; farx(n+1)=8.425854e+000; foe(n+1)=5.390346e+001; krok(n+1)=1.010485e-005; ng(n+1)=4.317422e+002;
n=216; farx(n+1)=8.488058e+000; foe(n+1)=5.382244e+001; krok(n+1)=1.738167e-005; ng(n+1)=3.306428e+002;
n=217; farx(n+1)=8.383577e+000; foe(n+1)=5.374247e+001; krok(n+1)=1.013046e-005; ng(n+1)=4.269882e+002;
n=218; farx(n+1)=8.444745e+000; foe(n+1)=5.366405e+001; krok(n+1)=1.709205e-005; ng(n+1)=3.272403e+002;
n=219; farx(n+1)=8.343159e+000; foe(n+1)=5.358663e+001; krok(n+1)=1.017076e-005; ng(n+1)=4.185749e+002;
n=220; farx(n+1)=8.403213e+000; foe(n+1)=5.351031e+001; krok(n+1)=1.700060e-005; ng(n+1)=3.232096e+002;
n=221; farx(n+1)=8.303946e+000; foe(n+1)=5.343500e+001; krok(n+1)=1.018440e-005; ng(n+1)=4.127277e+002;
n=222; farx(n+1)=8.363261e+000; foe(n+1)=5.336073e+001; krok(n+1)=1.697768e-005; ng(n+1)=3.193534e+002;
n=223; farx(n+1)=8.265980e+000; foe(n+1)=5.328718e+001; krok(n+1)=1.015759e-005; ng(n+1)=4.090558e+002;
n=224; farx(n+1)=8.324306e+000; foe(n+1)=5.321477e+001; krok(n+1)=1.695813e-005; ng(n+1)=3.153271e+002;
n=225; farx(n+1)=8.229142e+000; foe(n+1)=5.314302e+001; krok(n+1)=1.013760e-005; ng(n+1)=4.046199e+002;
n=226; farx(n+1)=8.286460e+000; foe(n+1)=5.307236e+001; krok(n+1)=1.695813e-005; ng(n+1)=3.112799e+002;
n=227; farx(n+1)=8.193428e+000; foe(n+1)=5.300235e+001; krok(n+1)=1.010485e-005; ng(n+1)=4.003817e+002;
n=228; farx(n+1)=8.249943e+000; foe(n+1)=5.293327e+001; krok(n+1)=1.707914e-005; ng(n+1)=3.070867e+002;
n=229; farx(n+1)=8.158149e+000; foe(n+1)=5.286465e+001; krok(n+1)=1.009346e-005; ng(n+1)=3.983985e+002;
n=230; farx(n+1)=8.213729e+000; foe(n+1)=5.279757e+001; krok(n+1)=1.684728e-005; ng(n+1)=3.039583e+002;
n=231; farx(n+1)=8.124628e+000; foe(n+1)=5.273096e+001; krok(n+1)=1.008355e-005; ng(n+1)=3.913492e+002;
n=232; farx(n+1)=8.179074e+000; foe(n+1)=5.266517e+001; krok(n+1)=1.697768e-005; ng(n+1)=2.995455e+002;
n=233; farx(n+1)=8.091276e+000; foe(n+1)=5.260000e+001; krok(n+1)=1.009346e-005; ng(n+1)=3.883402e+002;
n=234; farx(n+1)=8.144857e+000; foe(n+1)=5.253612e+001; krok(n+1)=1.672925e-005; ng(n+1)=2.966351e+002;
n=235; farx(n+1)=8.059480e+000; foe(n+1)=5.247284e+001; krok(n+1)=1.010485e-005; ng(n+1)=3.813747e+002;
n=236; farx(n+1)=8.112300e+000; foe(n+1)=5.241032e+001; krok(n+1)=1.684728e-005; ng(n+1)=2.927167e+002;
n=237; farx(n+1)=8.028040e+000; foe(n+1)=5.234827e+001; krok(n+1)=1.009346e-005; ng(n+1)=3.794364e+002;
n=238; farx(n+1)=8.079976e+000; foe(n+1)=5.228745e+001; krok(n+1)=1.663384e-005; ng(n+1)=2.897809e+002;
n=239; farx(n+1)=7.997948e+000; foe(n+1)=5.222715e+001; krok(n+1)=1.010485e-005; ng(n+1)=3.729155e+002;
n=240; farx(n+1)=8.048932e+000; foe(n+1)=5.216761e+001; krok(n+1)=1.666649e-005; ng(n+1)=2.860481e+002;
n=241; farx(n+1)=7.968546e+000; foe(n+1)=5.210868e+001; krok(n+1)=1.008355e-005; ng(n+1)=3.692351e+002;
n=242; farx(n+1)=8.018653e+000; foe(n+1)=5.205044e+001; krok(n+1)=1.666649e-005; ng(n+1)=2.825360e+002;
n=243; farx(n+1)=7.939624e+000; foe(n+1)=5.199286e+001; krok(n+1)=1.010255e-005; ng(n+1)=3.653948e+002;
n=244; farx(n+1)=7.988978e+000; foe(n+1)=5.193615e+001; krok(n+1)=1.649119e-005; ng(n+1)=2.797402e+002;
n=245; farx(n+1)=7.911647e+000; foe(n+1)=5.188007e+001; krok(n+1)=1.013760e-005; ng(n+1)=3.598670e+002;
n=246; farx(n+1)=7.960344e+000; foe(n+1)=5.182475e+001; krok(n+1)=1.641304e-005; ng(n+1)=2.768295e+002;
n=247; farx(n+1)=7.884411e+000; foe(n+1)=5.176996e+001; krok(n+1)=1.015216e-005; ng(n+1)=3.558684e+002;
n=248; farx(n+1)=7.932474e+000; foe(n+1)=5.171599e+001; krok(n+1)=1.632874e-005; ng(n+1)=2.741866e+002;
n=249; farx(n+1)=7.858015e+000; foe(n+1)=5.166246e+001; krok(n+1)=1.015216e-005; ng(n+1)=3.519013e+002;
n=250; farx(n+1)=7.905239e+000; foe(n+1)=5.160969e+001; krok(n+1)=1.625849e-005; ng(n+1)=2.712173e+002;
n=251; farx(n+1)=7.832246e+000; foe(n+1)=5.155743e+001; krok(n+1)=1.017255e-005; ng(n+1)=3.473829e+002;
n=252; farx(n+1)=7.878588e+000; foe(n+1)=5.150588e+001; krok(n+1)=1.611741e-005; ng(n+1)=2.685291e+002;
n=253; farx(n+1)=7.806657e+000; foe(n+1)=5.145497e+001; krok(n+1)=1.029142e-005; ng(n+1)=3.418868e+002;
n=254; farx(n+1)=7.852725e+000; foe(n+1)=5.140504e+001; krok(n+1)=1.576802e-005; ng(n+1)=2.677444e+002;
n=255; farx(n+1)=7.782656e+000; foe(n+1)=5.135546e+001; krok(n+1)=1.030147e-005; ng(n+1)=3.361065e+002;
n=256; farx(n+1)=7.827913e+000; foe(n+1)=5.130644e+001; krok(n+1)=1.583652e-005; ng(n+1)=2.642667e+002;
n=257; farx(n+1)=7.759298e+000; foe(n+1)=5.125786e+001; krok(n+1)=1.023989e-005; ng(n+1)=3.333930e+002;
n=258; farx(n+1)=7.804003e+000; foe(n+1)=5.120964e+001; krok(n+1)=1.607373e-005; ng(n+1)=2.605150e+002;
n=259; farx(n+1)=7.736166e+000; foe(n+1)=5.116176e+001; krok(n+1)=1.015968e-005; ng(n+1)=3.335055e+002;
n=260; farx(n+1)=7.780139e+000; foe(n+1)=5.111448e+001; krok(n+1)=1.608814e-005; ng(n+1)=2.575374e+002;
n=261; farx(n+1)=7.713455e+000; foe(n+1)=5.106760e+001; krok(n+1)=1.015759e-005; ng(n+1)=3.303595e+002;
n=262; farx(n+1)=7.756761e+000; foe(n+1)=5.102134e+001; krok(n+1)=1.602405e-005; ng(n+1)=2.549774e+002;
n=263; farx(n+1)=7.691284e+000; foe(n+1)=5.097550e+001; krok(n+1)=1.017255e-005; ng(n+1)=3.265601e+002;
n=264; farx(n+1)=7.733964e+000; foe(n+1)=5.093025e+001; krok(n+1)=1.595012e-005; ng(n+1)=2.525393e+002;
n=265; farx(n+1)=7.669626e+000; foe(n+1)=5.088541e+001; krok(n+1)=1.019428e-005; ng(n+1)=3.228130e+002;
n=266; farx(n+1)=7.711683e+000; foe(n+1)=5.084116e+001; krok(n+1)=1.584198e-005; ng(n+1)=2.502645e+002;
n=267; farx(n+1)=7.648537e+000; foe(n+1)=5.079734e+001; krok(n+1)=1.022121e-005; ng(n+1)=3.187062e+002;
n=268; farx(n+1)=7.690132e+000; foe(n+1)=5.075404e+001; krok(n+1)=1.580481e-005; ng(n+1)=2.479649e+002;
n=269; farx(n+1)=7.628171e+000; foe(n+1)=5.071108e+001; krok(n+1)=1.018440e-005; ng(n+1)=3.159898e+002;
n=270; farx(n+1)=7.669095e+000; foe(n+1)=5.066853e+001; krok(n+1)=1.588017e-005; ng(n+1)=2.449614e+002;
n=271; farx(n+1)=7.607971e+000; foe(n+1)=5.062638e+001; krok(n+1)=1.018462e-005; ng(n+1)=3.137144e+002;
n=272; farx(n+1)=7.648397e+000; foe(n+1)=5.058476e+001; krok(n+1)=1.579136e-005; ng(n+1)=2.428560e+002;
n=273; farx(n+1)=7.588239e+000; foe(n+1)=5.054350e+001; krok(n+1)=1.021265e-005; ng(n+1)=3.103084e+002;
n=274; farx(n+1)=7.628147e+000; foe(n+1)=5.050280e+001; krok(n+1)=1.567018e-005; ng(n+1)=2.408752e+002;
n=275; farx(n+1)=7.569069e+000; foe(n+1)=5.046243e+001; krok(n+1)=1.023758e-005; ng(n+1)=3.064961e+002;
n=276; farx(n+1)=7.608455e+000; foe(n+1)=5.042257e+001; krok(n+1)=1.560740e-005; ng(n+1)=2.386879e+002;
n=277; farx(n+1)=7.550132e+000; foe(n+1)=5.038303e+001; krok(n+1)=1.029142e-005; ng(n+1)=3.033028e+002;
n=278; farx(n+1)=7.588989e+000; foe(n+1)=5.034414e+001; krok(n+1)=1.535224e-005; ng(n+1)=2.372743e+002;
n=279; farx(n+1)=7.531888e+000; foe(n+1)=5.030556e+001; krok(n+1)=1.035782e-005; ng(n+1)=2.980525e+002;
n=280; farx(n+1)=7.570283e+000; foe(n+1)=5.026747e+001; krok(n+1)=1.527639e-005; ng(n+1)=2.353223e+002;
n=281; farx(n+1)=7.514030e+000; foe(n+1)=5.022967e+001; krok(n+1)=1.038706e-005; ng(n+1)=2.949788e+002;
n=282; farx(n+1)=7.551938e+000; foe(n+1)=5.019238e+001; krok(n+1)=1.515663e-005; ng(n+1)=2.336024e+002;
n=283; farx(n+1)=7.496442e+000; foe(n+1)=5.015537e+001; krok(n+1)=1.046074e-005; ng(n+1)=2.913558e+002;
n=284; farx(n+1)=7.533846e+000; foe(n+1)=5.011895e+001; krok(n+1)=1.489745e-005; ng(n+1)=2.325065e+002;
n=285; farx(n+1)=7.479543e+000; foe(n+1)=5.008282e+001; krok(n+1)=1.052881e-005; ng(n+1)=2.861831e+002;
n=286; farx(n+1)=7.516492e+000; foe(n+1)=5.004710e+001; krok(n+1)=1.484784e-005; ng(n+1)=2.305699e+002;
n=287; farx(n+1)=7.462980e+000; foe(n+1)=5.001165e+001; krok(n+1)=1.054828e-005; ng(n+1)=2.834570e+002;
n=288; farx(n+1)=7.499506e+000; foe(n+1)=4.997663e+001; krok(n+1)=1.477879e-005; ng(n+1)=2.288253e+002;
n=289; farx(n+1)=7.446587e+000; foe(n+1)=4.994186e+001; krok(n+1)=1.061240e-005; ng(n+1)=2.806490e+002;
n=290; farx(n+1)=7.482730e+000; foe(n+1)=4.990762e+001; krok(n+1)=1.455793e-005; ng(n+1)=2.279972e+002;
n=291; farx(n+1)=7.431260e+000; foe(n+1)=4.987360e+001; krok(n+1)=1.057446e-005; ng(n+1)=2.763843e+002;
n=292; farx(n+1)=7.466880e+000; foe(n+1)=4.983977e+001; krok(n+1)=1.484241e-005; ng(n+1)=2.242961e+002;
n=293; farx(n+1)=7.415844e+000; foe(n+1)=4.980618e+001; krok(n+1)=1.048504e-005; ng(n+1)=2.769332e+002;
n=294; farx(n+1)=7.450988e+000; foe(n+1)=4.977291e+001; krok(n+1)=1.486954e-005; ng(n+1)=2.220562e+002;
n=295; farx(n+1)=7.400497e+000; foe(n+1)=4.973991e+001; krok(n+1)=1.051121e-005; ng(n+1)=2.749197e+002;
n=296; farx(n+1)=7.435357e+000; foe(n+1)=4.970733e+001; krok(n+1)=1.476011e-005; ng(n+1)=2.209510e+002;
n=297; farx(n+1)=7.385812e+000; foe(n+1)=4.967494e+001; krok(n+1)=1.048504e-005; ng(n+1)=2.722557e+002;
n=298; farx(n+1)=7.420161e+000; foe(n+1)=4.964285e+001; krok(n+1)=1.484241e-005; ng(n+1)=2.184943e+002;
n=299; farx(n+1)=7.371261e+000; foe(n+1)=4.961099e+001; krok(n+1)=1.046329e-005; ng(n+1)=2.706573e+002;
n=300; farx(n+1)=7.405206e+000; foe(n+1)=4.957945e+001; krok(n+1)=1.484032e-005; ng(n+1)=2.166629e+002;

figure; semilogy(farx,'b'); hold on; semilogy(foe,'r'); xlabel('Iteracje'); ylabel('Earx, Eoe'); legend('Earx','Eoe'); title('Uczenie predyktora OE');
figure; subplot(2,1,1); semilogy(krok); xlabel('Iteracje'); ylabel('Krok');
subplot(2,1,2); semilogy(ng); xlabel('Iteracje'); ylabel('Norma gradientu');
Earx=farx(n+1)
Eoe=foe(n+1)
