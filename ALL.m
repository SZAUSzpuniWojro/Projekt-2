clearvars -except proba
proba = proba+1;

%ZMIENNE OBIEKTU
a1 = -1.535262;
a2 = 0.586646;
b1 = 0.027970;
b2 = 0.023414;


%CHARAKTERYSTYKA STATYCZNA%
%-------------------------%
u = -1:0.01:1;
y = zeros(1,length(u));
for i = 1:length(u)
    y(i) = g2((g1(u(i))*(b1+b2)/(1+a1+a2)));
end

figure_number = 1;

figure(figure_number)
plot(u,y);
grid on;
title("Charakterystyka statyczna obiektu")
xlabel("u")
ylabel("y")


%SYMULACJA OBIEKTU%
%-----------------%
clearvars u y
tsim = 2000;
save_mode = 0; %informuje, czy nadpisywaæ dane.txt

for set = 1:2
    if set == 1
        rng(3); %uczenie
    end
    if set == 2
        rng(2);  %weryfikacja
    end
    
    x = zeros(2, tsim);
    y = zeros(tsim, 1);
    
    
    u = zeros(tsim, 1);
    itim = 50;
    for k = 1:itim:tsim
        u(k:k+itim-1) = 1 - 2*rand(1);
    end
    
    for k = 4:tsim
        x(1, k) = -a1*x(1,k-1) + x(2,k-1) + b1*g1(u(k-3));
        x(2, k) = -a2*x(1,k-1) + b2*g1(u(k-3));
        y(k) = g2(x(1,k)) + 0.03*rand(1);
    end
    
    data_set = [u , y];
    d_size = size(data_set);
    if set == 1
        UCZ = data_set;
        if save_mode == 1
            %zapisz dane w pliku "dane.txt"
            fid = fopen('dane.txt','w');
            for i = 1:d_size(1)
                fprintf(fid,'%d %d\n', data_set(i,:));
            end
            fclose(fid);
        end
    end
    
    if set == 2
        WER = data_set;
    end
    
end

figure_number = figure_number + 1;
figure(figure_number)
stairs(1:tsim, UCZ(:,1))
hold on;
grid on;
stairs(1:tsim, UCZ(:,2))
title("Symulacja obiektu: dane ucz¹ce")
xlabel("chwila k")
legend("u","y")

figure_number = figure_number + 1;
figure(figure_number)
stairs(1:tsim, WER(:,1))
hold on;
grid on;
stairs(1:tsim, WER(:,2))
title("Symulacja obiektu: dane weryfikuj¹ce")
xlabel("chwila k")
legend("u","y")


%TWORZENIE MODELU%
%----------------%
clearvars -except proba tsim UCZ WER figure_number

tau = 3;
nB = 4;
nA = 2;

new_model = 0;

if new_model == 1
    system('sieci.exe');
    model;
else
    load('OE_BFGS/OE_BFGS_7neurony4proba.mat'); %TODO
end


u_ucz = UCZ(:,1)'; 
u_wer = WER(:,1)';
y_true_ucz = UCZ(:,2)'; 
y_true_wer = WER(:,2)';

y_mod_ucz = zeros(1,tsim);
y_mod_wer = zeros(1,tsim);

if nB>nA
    i = nB + 1;
else
    i = nA + 1;
end

for k = i:tsim
    %dane ucz¹ce
    q_u = [flip(u_ucz(k-nB:k-tau)),flip(y_mod_ucz(k-nA:k-1))]';
    y_mod_ucz(k) = w20 + w2* tanh(w10 + w1*q_u);
    
    %dane weryfikuj¹ce
    q_w = [flip(u_wer(k-nB:k-tau)),flip(y_mod_wer(k-nA:k-1))]';
    y_mod_wer(k) = w20 + w2* tanh(w10 + w1*q_w);
end

figure_number = figure_number + 1;
figure(figure_number);
stairs(y_true_ucz);
hold on;
grid on;
stairs(y_mod_ucz);
title("Aproksymacja danych: zbiór ucz¹cy");
legend("y", "y_{mod}",'Location','southwest');
xlabel("Chwila k");
hold off;

figure_number = figure_number + 1;
figure(figure_number);
stairs(y_true_wer);
hold on;
grid on;
stairs(y_mod_wer);
title("Aproksymacja danych: zbiór weryfikuj¹cy");
legend("y", "y_{mod}",'Location','southwest');
xlabel("Chwila k");
hold off;


%liczenie b³êdu
err_ucz = (y_mod_ucz - y_true_ucz)*(y_mod_ucz - y_true_ucz)'
err_wer = (y_mod_wer - y_true_wer)*(y_mod_wer - y_true_wer)'
err = [err_ucz, err_wer];
%zapisywanie b³êdu
fid = fopen('bledy_OE_BFGS_2.txt', 'at');
fprintf(fid, '%d: %d %d\n',proba, err);
fclose(fid);

%saving part
if new_model == 1
    message= ['OE_BFGS_2/OE_BFGS_TESTneurony' num2str(proba) 'proba'];
    save(message); % PRZY PRZYWRACANIU WORKSPACE WYKOMENTOWAÆ
end


%REGULACJA%
%---------%
clearvars -except w1 w10 w2 w20 proba UCZ figure_number

%PARAMETRY
tsim=500;
nB = 4;
tau = 3;
nA = 2;
lbd = 10;
N = 20;
Nu = 1;
Lbd = lbd * eye(Nu);
delta = 10^(-5);

alfa1 = -1.535262;
alfa2 = 0.586646;
beta1 = 0.027970;
beta2 = 0.023414;

%budowanie sygna³u wejœciowego
y_zad1 = ones(1,tsim/4)*(0.7);
y_zad2 = ones(1,tsim/4)*(0.0);
y_zad3 = ones(1,tsim/4)*(-3.2);
y_zad4 = ones(1,tsim/4)*(0.3);
y_zad = [y_zad1, y_zad2, y_zad3, y_zad4];

for regulator = 1:2 %1 -GPC, 2 - NPL
    x = zeros(2, tsim);
    y = zeros(1, tsim);
    
    u = zeros(1,tsim);
    
    %LICZENIE OFF-LINE : GPC
    if regulator == 1
        
        u_ucz = UCZ(:, 1)';
        y_ucz = UCZ(:, 2)';
        q =  [u_ucz(2), u_ucz(1), y_ucz(4), y_ucz(3)]';
        q1 = [u_ucz(2)+delta, u_ucz(1), y_ucz(4), y_ucz(3)]';
        q2 = [u_ucz(2), u_ucz(1)+delta, y_ucz(4), y_ucz(3)]';
        q3 = [u_ucz(2), u_ucz(1), y_ucz(4)+delta, y_ucz(3)]';
        q4 = [u_ucz(2), u_ucz(1), y_ucz(4), y_ucz(3)+delta]';
        
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
    end
    
    %------g³ówna pêtla regulacji------
    for k=tau+2:tsim
        
        %symulacja obiektu
        x(1, k) = -alfa1*x(1,k-1) + x(2,k-1) + beta1*g1(u(k-3));
        x(2, k) = -alfa2*x(1,k-1) + beta2*g1(u(k-3));
        y(k) = g2(x(1,k));% + 0.03*rand(1);
        
        if regulator == 1
            y_mod = B(3)*u(k-3)+B(4)*u(k-4)-A(1)*y(k-1)-A(2)*y(k-2);
        end
        
        %LICZENIE ON-LINE : NPL
        if regulator == 2
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
        end
        
        d = y(k) - y_mod;
        
        %liczenie odpowiedzi skokowej
        Y0 = zeros(N,1);
        for i = 1:N
            if i == 1
                q0 = [u(k-2), u(k-3), y(k), y(k-1)]';
            elseif i == 2
                q0 = [u(k-1), u(k-2), Y0(1), y(k)]';
            else
                q0 = [u(k-1), u(k-1), Y0(i-1), Y0(i-2)]';
            end
            if regulator == 1
                Y0(i) = B(3)*q0(1) + B(4)*q0(2) - A(1)*q0(3) - A(2)*q0(4) + d;
            elseif regulator == 2
                Y0(i) = w20 + w2*tanh(w10+ w1*q0) + d;
            end
        end
        
        %------prawo regulacji------
        du = K*(y_zad(k)*ones(N,1)-Y0);
        u(k) = u(k-1) + du(1);
        u(k) = min(u(k), 1);
        u(k) = max(u(k), -1);
    end
    
    figure_number = figure_number + 1;
    figure(figure_number)
    plot(y)
    hold on;
    grid on;
    plot(y_zad);
    xlabel("Chwila k");
    ylabel("y");
    legend("y","y^{zad}",'Location','southwest')
    if regulator == 1
        title("Wyjœcie obiektu : regulator GPC")
    else
        title("Wyjœcie obiektu : regulator NPL")
    end
    
    figure_number = figure_number + 1;
    figure(figure_number)
    plot(u)
    grid on;
    xlabel("Chwila k");
    ylabel("u");
    if regulator == 1
        title("Sterowanie : GPC")
    else
        title("Sterowanie : NPL")
    end
end
