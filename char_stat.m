a1 = -1.535262;
a2 = 0.586646;
b1 = 0.027970;
b2 = 0.023414;

x = -1:0.1:1;
y = zeros(1,length(x));
for i = 1:length(x)
    y(i) = g2((g1(x(i))*(b1+b2)/(1+a1+a2)));
end

plot(x,y);