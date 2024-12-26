clear; close all;
% Parameters
xBinSize = 0.3/200; % Size of each bin along the x-axis
yBinSize = 5/200; % Size of each bin along the y-axis
csvFile = 'output.csv'; % Your CSV file

% Read the CSV file
data = csvread(csvFile);

% Compute the x and y axis ranges
[nRows, nCols] = size(data); % Get the dimensions of the image
x = 2.1+(0:nCols-1) * xBinSize; % x-axis values
y = (0:nRows-1) * yBinSize; % y-axis values

% Display the image with specified axes
imshow(data, [], 'XData', x, 'YData', y);
%colormap(gray); % Apply grayscale colormap
colorbar; % Optional: Display a color scale
axis on; % Turn on the axes
hold;
%
% MEDIANS

% Calculate the median x position for each row
[nRows, nCols] = size(data);
medianXIndices = zeros(nRows, 1); % Preallocate array for median x indices

for i = 1:nRows
    % Find the column index corresponding to the median of the row
    medianXIndices(i) = find(cumsum(data(i, :)) >= sum(data(i, :)) / 2, 1);
end

% Prepare the x and y coordinates for the line
xCoords = 2.1+(medianXIndices - 1) * xBinSize; % Convert indices to x-axis values
yCoords = y; % Y-axis values for each row

% Plot the line on top of the image
plot(xCoords, yCoords, 'r', 'LineWidth', 2); % Red line with thickness 2


ylabel('time (ns)');
xlabel('energy (eV)');
%title('Grayscale Image with Specified Axes');
daspect([4 80 1]);
axis tight;
set(gcf, 'Position', [100, 100, 1000, 800]); % [x, y, width, height]
