close all; clc;

[xDim yDim] = size(cSpace);
cSpaceAxes = [0 xDim 0 yDim];

figure; hold on;
axis (cSpaceAxes); 
image(100*(1-cSpace)');
colormap(gray);

for i=1:10
    [st en] = generatePoints(cSpace);
    plot(st(1), st(2), 'g*', 'MarkerSize', 10);
    plot(en(1), en(2), 'r*', 'MarkerSize', 10);
    plot([st(1) en(1)], [st(2) en(2)], 'b');
end