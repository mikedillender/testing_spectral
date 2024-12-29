clear; close all;
xBinSize = 0.25/200; % Size of each bin along the x-axis
yBinSize = 8/300; % Size of each bin along the y-axis

csvFile = 'output_3L_d100_C23.csv'; % Your CSV file
data1 = csvread(csvFile);

% Compute the x and y axis ranges
[nRows, nCols] = size(data1); % Get the dimensions of the image
x = 2.15+(0:nCols-1) * xBinSize; % x-axis values
y = (0:nRows-1) * yBinSize; % y-axis values

% Display the image with specified axes
subplot(1,2,1);
imshow(data1, [], 'XData', x, 'YData', y);
colorbar; axis on; hold on;

%% medians 1

% Calculate the median x position for each row
[nRows, nCols] = size(data1);
medianXIndices = zeros(nRows, 1); % Preallocate array for median x indices

for i = 1:nRows
    medianXIndices(i) = find(cumsum(data1(i, :)) >= sum(data1(i, :)) / 2, 1);
end

% Prepare the x and y coordinates for the line
xCoords = 2.15+(medianXIndices - 1) * xBinSize; % Convert indices to x-axis values

% Plot the line on top of the image
plot(xCoords, y, 'r', 'LineWidth', 2); % Red line with thickness 2
title("2 Full Layers");


ylabel('time (ns)');
xlabel('energy (eV)');
daspect([4 100 1]);
axis tight;
xlim([2.18 2.33]);
set(gcf, 'Position', [100, 100, 1000, 800]); % [x, y, width, height]

%%  part 2
csvFile2 = 'output2_3L_d100_C23.csv'; % Your CSV file
data2 = csvread(csvFile2);
subplot(1,2,2);

medianXIndices2 = zeros(nRows, 1); % Preallocate array for median x indices
for i = 1:nRows
    medianXIndices2(i) = find(cumsum(data2(i, :)) >= sum(data2(i, :)) / 2, 1);
end
xCoords2 = 2.15+(medianXIndices2 - 1) * xBinSize; % Convert indices to x-axis values
subplot(1,2,2);

imshow(data2, [], 'XData', x, 'YData', y);
colorbar; axis on; hold on;
plot(xCoords2, y, 'r', 'LineWidth', 2); % Red line with thickness 2
plot(xCoords, y, 'g', 'LineWidth', 1); % Red line with thickness 2
legend("with islands",'without islands',Location='southeast')

ylabel('time (ns)');
xlabel('energy (eV)');
daspect([4 100 1]);
axis tight;
xlim([2.18 2.33]);

title("2 Full Layers + 20% 3rd");
set(gcf, 'Position', [100, 100, 1000, 800]); % [x, y, width, height]