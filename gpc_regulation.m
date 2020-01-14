%symulacja obiektu -> pomiar y(k)
%linearyzacja modelu -> obliczyæ a,b
%obliczyæ odp skokow¹ -> sj(k), j=1, N - za pomoc¹ modelu zlinearyzowanego
%macierz dynamiczna M(k)
%K(k) = (M'(k)M(k) + Lam*I)^(-1)*M'(k)
%odpowiedŸ swobodn¹ na podstawie SN
%obliczyæ sterowanie -> du(k) = K(k)(Y_zad(k) - Y0(k))
%u(k) = u(k-1) + du(k|k)
%load('OE_BFGS\OE_BFGS_7neurony4proba.mat')
%main;
%model;
clearvars -except w1 w10 w2 w20 proba UCZ

tsim=500;
nB = 4;
tau = 3;
nA = 2;
lbd = 10;
N = 20;
Nu = 20;
Lbd = lbd * eye(Nu);

alfa1 = -1.535262;
alfa2 = 0.586646;
beta1 = 0.027970;
beta2 = 0.023414;

x = zeros(2, tsim);
y = zeros(1, tsim);

delta = 10^(-5);

y_zad1 = ones(1,tsim/4)*(0.7);
y_zad2 = ones(1,tsim/4)*(0.0);
y_zad3 = ones(1,tsim/4)*(-3.2);
y_zad4 = ones(1,tsim/4)*(0.3);
y_zad = [y_zad1, y_zad2, y_zad3, y_zad4];
u = zeros(1,tsim);

%linearyzacja modelu

dane_ucz = dlmread('dane.txt');
x_ucz = UCZ(:, 1)';
y_ucz = UCZ(:, 2)';
q =  [x_ucz(2), x_ucz(1), y_ucz(4), y_ucz(3)]';
q1 = [x_ucz(2)+delta, x_ucz(1), y_ucz(4), y_ucz(3)]';
q2 = [x_ucz(2), x_ucz(1)+delta, y_ucz(4), y_ucz(3)]';
q3 = [x_ucz(2), x_ucz(1), y_ucz(4)+delta, y_ucz(3)]';
q4 = [x_ucz(2), x_ucz(1), y_ucz(4), y_ucz(3)+delta]';

y_mod = w20 + w2*tanh(w10+ w1*q);
y_mod1 = w20 + w2*tanh(w10+ w1*q1);
y_mod2 = w20 + w2*tanh(w10+ w1*q2);
y_mod3 = w20 + w2*tanh(w10+ w1*q3);
y_mod4 = w20 + w2*tanh(w10+ w1*q4);

B3 = (y_mod1 - y_mod)/delta;
B4 = (y_mod2 - y_mod)/delta;
A1 = -(y_mod3 - y_mod)/delta;
A2 = -(y_mod4 - y_mod)/delta;

B = [0,0,B3,B4];
A = [A1,A2];


%odpowiedŸ skokowa
M = zeros(N,Nu);
s = zeros(N,1);

for j=1:N
    A_prim = 0;
    if(j-1)>0
        for i=1:min([j-1,nA])
            A_prim = A_prim + A(i)*s(j-i);
        end
    end
    
    s(j) = sum(B(1:min([j,4]))) - A_prim;
end

%macierz dynamiczna
for i=1:1:Nu
    M(i:N,i)=s(1:N-i+1)';
end

%wspó³czynnik K
K = (M'*M + Lbd)^(-1)*M';

for k=tau+2:tsim
    x(1, k) = -alfa1*x(1,k-1) + x(2,k-1) + beta1*g1(u(k-3));
    x(2, k) = -alfa2*x(1,k-1) + beta2*g1(u(k-3));
    y(k) = g2(x(1,k));% + 0.03*rand(1);
  
    y_mod = B(3)*u(k-3)+B(4)*u(k-4)-A(1)*y(k-1)-A(2)*y(k-2);
    d = y(k) - y_mod;
    
    Y0 = zeros(N,1);
    for i = 1:N
        if i == 1
            q0 = [u(k-2), u(k-3), y(k), y(k-1)]';
        elseif i == 2
            q0 = [u(k-1), u(k-2), Y0(1), y(k)]';
        else
            q0 = [u(k-1), u(k-1), Y0(i-1), Y0(i-2)]';
        end
        Y0(i) = B(3)*q0(1) + B(4)*q0(2) - A(1)*q0(3) - A(2)*q0(4) + d;
    end
    %%prawo regulacji
    du = K*(y_zad(k)*ones(N,1)-Y0);
    u(k) = u(k-1) + du(1);
    u(k) = min(u(k), 1);
    u(k) = max(u(k), -1);
end


figure(1)
plot(y)
hold on;
plot(y_zad);

figure(2)
plot(u)