clear; close all;
xBinSize = 0.24/100; % Size of each bin along the x-axis
yBinSize = 0.24/100; % Size of each bin along the y-axis

csvFile = 'E_rates_L180_d173_C100.csv'; 
data1 = (csvread(csvFile)+1e-8)';

% Compute the x and y axis ranges
[nRows, nCols] = size(data1); % Get the dimensions of the image
x = 2.13+(0:nCols-1) * xBinSize; % x-axis values
y = 2.13+(0:nCols-1) * yBinSize; % x-axis values

% Display the image with specified axes
figure();
data1(data1<0)=-.01;
imagesc(x,y,data1);
myColorMap = parula(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar; axis on; hold on;
set(gca,'YDir','normal');
clim([-.01,.4])
title("Average Transfer Rates from Simulation")
ylabel("E_A [eV]");xlabel("E_D [mV]");
%{
xlim(2.25+.12*[-1 1])
ylim(2.25+.12*[-1 1])
%}
