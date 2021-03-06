%execute sieci.exe
clearvars -except proba

system('sieci.exe'); % PRZY PRZYWRACANIU WORKSPACE WYKOMENTOWAĆ
save_mode = 0;
przebieg;

tau = 3;
nB = 4;
nA = 2;
model; % PRZY PRZYWRACANIU WORKSPACE WYKOMENTOWAĆ
%load('OE_BFGS/OE_BFGS_7neurony4proba.mat'); %PRZY UCZENIU WYKOMENTOWAĆ
%wyznaczanie wyjścia sieci neuronowej

u_ucz = UCZ(:,1)'; 
u_wer = WER(:,1)';
y_true_ucz = UCZ(:,2)'; 
y_true_wer = WER(:,2)';

y_mod_ucz = zeros(1,2000);
y_mod_wer = zeros(1,2000);

if nB>nA
    i = nB + 1;
else
    i = nA + 1;
end

for k = i:tsim
    %dane uczące
    q_u = [flip(u_ucz(k-nB:k-tau)),flip(y_mod_ucz(k-nA:k-1))]'; 
    y_mod_ucz(k) = w20 + w2* tanh(w10 + w1*q_u);
    
    %dane weryfikujące
    q_w = [flip(u_wer(k-nB:k-tau)),flip(y_mod_wer(k-nA:k-1))]';
    y_mod_wer(k) = w20 + w2* tanh(w10 + w1*q_w);
end

figure(1);
plot(y_mod_ucz);
hold on;
plot(y_true_ucz);
hold off;

figure(2);
plot(y_mod_wer);
hold on;
plot(y_true_wer);
hold off;
%liczenie błędu

err_ucz = (y_mod_ucz - y_true_ucz)*(y_mod_ucz - y_true_ucz)'
err_wer = (y_mod_wer - y_true_wer)*(y_mod_wer - y_true_wer)'
err = [err_ucz, err_wer];

fid = fopen('bledy_OE_BFGS_2.txt', 'at');
fprintf(fid, 'TEST N: %d %d\n', err);
fclose(fid);

%saving part
%proba must be declared in the command line
proba = proba+1; 
% if proba == 6
%    proba = 1;
% end
message= ['OE_BFGS_2/OE_BFGS_TESTneurony' num2str(proba) 'proba'];
save(message); % PRZY PRZYWRACANIU WORKSPACE WYKOMENTOWAĆ

