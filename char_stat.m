clearvars -except proba
a1 = -1.535262;
a2 = 0.586646;
b1 = 0.027970;
b2 = 0.023414;

u = -1:0.1:1;
y = zeros(1,length(u));
for i = 1:length(u)
    y(i) = g2((g1(u(i))*(b1+b2)/(1+a1+a2)));
end

plot(u,y);
