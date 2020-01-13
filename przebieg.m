%clear all;

a1 = -1.535262;
a2 = 0.586646;
b1 = 0.027970;
b2 = 0.023414;


for set = 1:2
    if set == 1
        rng(16); %UCZENIE
    end
    if set == 2
        rng(2);  %WARYFIKACJA
    end
    
    tsim = 4000;  % GLOBALS
    
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
    
    %T = table(u , y);
    data_set = [u , y];
    d_size = size(data_set);
    if set == 1
        UCZ = data_set;
        if save_mode == 1
            %write data set to file "dane.txt"
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

% figure(1)
% plot(1:tsim, UCZ(:,1))
% hold on;
% plot(1:tsim, UCZ(:,2))
% 
% figure(2)
% plot(1:tsim, WER(:,1))
% hold on;
% plot(1:tsim, WER(:,2))
