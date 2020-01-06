%clear all;

rng(2);

a1 = -1.535262;
a2 = 0.586646;
b1 = 0.027970;
b2 = 0.023414;

k_sim = 2000;
x = zeros(2, k_sim);
y = zeros(1, k_sim);

%TODO: Generate amplitude signal
u = zeros(1, k_sim);
itim = 50;
for k = 1:itim:k_sim
    u(k:k+itim-1) = 1 - 2*rand(1);
end

for k = 4:k_sim
    x(1, k) = -a1*x(1,k-1) + x(2,k-1) + b1*g1(u(k-3));
    x(2, k) = -a2*x(1,k-1) + b2*g1(u(k-3));
    y(k) = g2(x(1,k)) + 0.03*rand(1);
end

plot(1:k_sim, y)