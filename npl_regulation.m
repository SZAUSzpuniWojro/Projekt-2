load('OE_BFGS\OE_BFGS_7neurony4proba.mat')
%main;
%model;
clearvars -except w1 w10 w2 w20 proba

tsim=500;
nB = 4;
tau = 3;
nA = 2;
lbd = 10;
N = 20;
Nu = 1;
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

for k=tau+2:tsim
    x(1, k) = -alfa1*x(1,k-1) + x(2,k-1) + beta1*g1(u(k-3));
    x(2, k) = -alfa2*x(1,k-1) + beta2*g1(u(k-3));
    y(k) = g2(x(1,k));% + 0.03*rand(1);
    
    %linearyzacja modelu
    q =  [u(k-3), u(k-4), y(k-1), y(k-2)]';
    q1 = [u(k-3)+delta, u(k-4), y(k-1), y(k-2)]';
    q2 = [u(k-3), u(k-4)+delta, y(k-1), y(k-2)]';
    q3 = [u(k-3), u(k-4), y(k-1)+delta, y(k-2)]';
    q4 = [u(k-3), u(k-4), y(k-1), y(k-2)+delta]';
    
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
    
    %odpowied� skokowa
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
    
    %wsp�czynnik K
    K = (M'*M + Lbd)^(-1)*M';
    
    
    %odpowied� swobodna
    Y0 = zeros(N,1);
    d = y(k) - y_mod;
    for i = 1:N
        if i == 1
            q0 = [u(k-2), u(k-3), y(k), y(k-1)]';
        elseif i == 2
            q0 = [u(k-1), u(k-2), Y0(1), y(k)]';
        else
            q0 = [u(k-1), u(k-1), Y0(i-1), Y0(i-2)]';
        end
        Y0(i) = w20 + w2*tanh(w10+ w1*q0) + d;
    end
    
    %obliczenie sterowania
    Y_zad = ones(N,1)*y_zad(k);
    du = K*(Y_zad - Y0);
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