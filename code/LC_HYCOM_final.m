%%ReadMe

%This script was developed by Dan Holstein and Gaby Carpenter in the
%Seascape Ecology Lab @LSU for the manuscript "Mesoscale circulation and
%coral disease risk at isolated reefs"

%1. Download HYCOM data and extract a daily 0.17m SSH contour + clean
%results

    % HYCOM-TSIS 1/100º Gulf of Mexico Reanalysis (HYCOM-TSIS GOMb0.01)
    % https://www.hycom.org/data/gomb0pt01/gom-reanalysis

%2. Plot daily Loop Current (LC) contours - colored by LC length and
%max/min for each "disease period" of interest. Create histogram of loop
%current lengths


%% Step 1: Download HYCOM data and extract a daily 0.17M SSH contour

%Script to take the LC for the first hour that matches criteria

%There is a slight change in the url on day 152 in 2017. Need to download
%that year is 2 batches and adjust '/010_archv.' to '/023_archv.' on line
%44

options = weboptions('Timeout', 200);
LCresults = struct();
t = 0;
LatLim = 21; % Northern latitude limit for contour filtering
LonLim = -80.0; % Western longitude limit for contour filtering
for y = 2001:2023
    tic
    for d = 1:yeardays(y)  % Loop through each day of the year - can replace "yeardays" with number if needed
        t = t + 1;
        delete HYCOM.nc
        for h = 0:23  % Loop through each hour of the day
            
            tic
            Dt = datetime(y, 1, d, h, 0, 0); % Date with hour
            % Format the hour correctly (00 to 23)
            hourStr = sprintf('%02d', h);
            % Construct the URL with the correct hour
            url = strcat('http://ncss.hycom.org/thredds/ncss/datasets/GOMb0.01/reanalysis/data/', ...
                sprintf('%04d', y), '/010_archv.', sprintf('%04d', y), '_', ...
                sprintf('%03d', d), '_', hourStr, '_2d.nc?var=ssh&disableLLSubset=on&disableProjSubset=on&horizStride=1&time=', ...
                sprintf('%04d', y), '-', sprintf('%02d', month(Dt)), '-', ...
                sprintf('%02d', day(Dt)), 'T', hourStr, '%3A00%3A00Z&accept=netcdf4');
            % Download and read data
            temp = websave('HYCOM.nc', url, options);
            toc
            SSH = ncread('HYCOM.nc', 'ssh');
            Lon = ncread('HYCOM.nc', 'Longitude');
            Lat = ncread('HYCOM.nc', 'Latitude');
            % Compute contour lines
            C2 = contourc(double(Lon), double(Lat), SSH', [0.17 0.17]);
            myind = find(C2(1,:) == 0.17);
            myvects = struct();
            for i = 1:length(myind)
                myvects(i).X = C2(1, myind(i)+1 : myind(i)+C2(2, myind(i)));
                myvects(i).Y = C2(2, myind(i)+1 : myind(i)+C2(2, myind(i)));
                myvects(i).Length = linelength(geolineshape(myvects(i).Y, myvects(i).X));
                % plot(myvects(i).X, myvects(i).Y); hold on
            end
            [B, I] = sort([myvects.Length], 'descend');
            for l = 1:length(I)
                if any([myvects(I(l)).Y] < LatLim) && any([myvects(I(l)).X] > LonLim)
                    if any([myvects(I(l)).Y] > 23) % Ensure contour runs north of Cuba
                        LClongest = I(l);
                        break;
                    end
                end
            end
            % Store results in structure
            LCresults(t).Year = year(Dt);
            LCresults(t).Month = month(Dt);
            LCresults(t).Day = day(Dt);
            LCresults(t).Hour = hour(Dt);  % Store the hour
            LCresults(t).DayOfYear = d;
            
            if exist('LClongest', 'var')
                LCresults(t).X = myvects(LClongest).X;
                LCresults(t).Y = myvects(LClongest).Y;
                LCresults(t).Length = myvects(LClongest).Length;
                clear LClongest
                break % Dan thinks this will break the h loop and go back to the d loop
            end
            % Cleanup
            
            
            toc
        end
    end
    toc
end
save Results LCresults

% Loop to create a datetime field
for i = 1:numel(LCresults)
    LCresults(i).Date = datetime(LCresults(i).Year, LCresults(i).Month, LCresults(i).Day);
end

% Define the expected date range
startDate = datetime(2001,1,1);
endDate = datetime(2023,12,31);
expectedDates = startDate:endDate;

% Convert structure dates to datetime format
actualDates = datetime([LCresults.Year], [LCresults.Month], [LCresults.Day]);

% Find missing dates
missingDates = setdiff(expectedDates, actualDates);

% Find duplicate dates
[uniqueDates, ~, idx] = unique(actualDates);
duplicateCounts = accumarray(idx, 1);
duplicateDates = uniqueDates(duplicateCounts > 1);

% Display results
disp('Missing Dates:');
disp(missingDates);

disp('Duplicate Dates:');
disp(duplicateDates);

% Remove duplicate dates
[uniqueDates, uniqueIndices] = unique(actualDates, 'stable'); % Keep first occurrence
LCresults = LCresults(uniqueIndices); % Retain only unique entries

% Sort the structure based on dates
[sortedDates, sortIndices] = sort(actualDates);
LCresults = LCresults(sortIndices); % Sort the structure

% Count rows without X data
numMissingX = sum(arrayfun(@(s) isempty(s.X), LCresults));
fprintf('Number of rows without X data: %d\n', numMissingX);



%%Clean up data by cutting LC at YC and SF

% Cutting
for i = 1:length(LCresults)
    xtest = LCresults(i).X < -81.5; %Playing around with this because of things in N FL
    ytest = LCresults(i).Y > 21.7;
    LCresults(i).Xcut = LCresults(i).X(xtest & ytest);
    LCresults(i).Ycut = LCresults(i).Y(xtest & ytest);
    LCresults(i).CutLength = linelength(geolineshape(LCresults(i).Ycut, LCresults(i).Xcut));
end

%%Clean up data by removing outliers

% Initialize an empty array to store matching indices
matchingIndices = [];

% Loop through each row in the structure
for i = 1:numel(LCresults)
    % Extract the latitude and longitude matrix for this row
    latMatrix = LCresults(i).Ycut;
    lonMatrix = LCresults(i).Xcut;
    
    % Check if any value in the matrix satisfies the condition
    if any(lonMatrix(:) < -95 & latMatrix(:) < 25)  % Flatten to a vector using (:)
        matchingIndices = [matchingIndices, i];  % Store the row index
    end
end

% Extract the filtered rows from LCresults
filteredRows = LCresults(matchingIndices);

%Plot these to visualize

figure()
cmap = parula(length(filteredRows));
worldmap([18 32],[-100 -80]) 
colorbar
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
hold on
for i = 1:1   %length(filteredRows)
    if ~isempty(filteredRows(i).Xcut)
        plotm(filteredRows(i).Ycut,filteredRows(i).Xcut,'Color',cmap(i,:)); 
    end
end

hold off

%%REMOVE the rows that are out of the LC range
% Extract dates from filteredRows
filteredDates = [filteredRows.Date];  % Assuming there is a 'Date' field

% Extract all dates from LCresults
allDates = [LCresults.Date];

% Find indices of rows in LCresults that match any date in filteredRows
matchingIndices = ismember(allDates, filteredDates);

% Remove matching rows
LCresultsCleaned = LCresults(~matchingIndices);

LCresults = LCresultsCleaned

%ROUND 2 of data cleaning
%Find and remove rows with a length of zero
% Find the indices of rows where CutLength is not zero
zeros = [LCresults.CutLength] ~= 0;

% Keep only the valid rows
LCresults = LCresults(zeros);



%%FIND CONTOURS THAT DO NOT GO THROUGH SF

% Initialize an empty array to store matching indices
matchingIndices2 = [];

% Loop through each row in the structure
for i = 1:numel(LCresults)
    % Extract the latitude and longitude matrix for this row
    latMatrix = LCresults(i).Ycut;
    lonMatrix = LCresults(i).Xcut;
    
    % Check if any value in the matrix satisfies the condition
    if all(lonMatrix(:) <= -81.6)  % Flatten to a vector using (:)
        matchingIndices2 = [matchingIndices2, i];  % Store the row index
    end
end

% Extract the filtered rows from LCresults
filteredRows2 = LCresults(matchingIndices2);

%Plot these to visualize

figure()
cmap = parula(length(filteredRows2));
worldmap([18 32],[-100 -80]) 
colorbar
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
hold on
for i = 1:length(filteredRows2)
    if ~isempty(filteredRows2(i).Xcut)
        plotm(filteredRows2(i).Ycut,filteredRows2(i).Xcut,'Color',cmap(i,:)); 
    end
end

hold off

%%REMOVE the rows 
% Extract dates from filteredRows
filteredDates = [filteredRows2.Date];  % Assuming there is a 'Date' field

% Extract all dates from LCresults
allDates = [LCresults.Date];

% Find indices of rows in LCresults that match any date in filteredRows
matchingIndices2 = ismember(allDates, filteredDates);

% Remove matching rows
LCresultsCleaned = LCresults(~matchingIndices2);

LCresults = LCresultsCleaned

%%REMOVE the rows that are below CUBA 

% Initialize an empty array to store matching indices
matchingIndices3 = [];

% Loop through each row in the structure
for i = 1:numel(LCresults)
    % Extract the latitude and longitude matrix for this row
    latMatrix = LCresults(i).Ycut;
    lonMatrix = LCresults(i).Xcut;
    
    % Check if any value in the matrix satisfies the condition
    if any(lonMatrix(:) < -91 & latMatrix(:) < 22)  % Flatten to a vector using (:)
        matchingIndices3 = [matchingIndices3, i];  % Store the row index
    end
end

% Extract the filtered rows from LCresults
filteredRows3 = LCresults(matchingIndices3);

%%REMOVE the rows 
% Extract dates from filteredRows
filteredDates = [filteredRows3.Date];  % Assuming there is a 'Date' field

% Extract all dates from LCresults
allDates = [LCresults.Date];

% Find indices of rows in LCresults that match any date in filteredRows
matchingIndices3 = ismember(allDates, filteredDates);

% Remove matching rows
LCresultsCleaned = LCresults(~matchingIndices3);

LCresults = LCresultsCleaned

%%REMOVE the rows that start in GOM 

% Initialize an empty array to store matching indices
matchingIndices4 = [];

% Loop through each row in the structure
for i = 1:numel(LCresults)
    % Extract the latitude and longitude matrix for this row
    latMatrix = LCresults(i).Ycut;
    lonMatrix = LCresults(i).Xcut;
    
    % Check if any value in the matrix satisfies the condition
    if any(lonMatrix(:) > -84 & latMatrix(:) < 21.8)  % Flatten to a vector using (:)
        matchingIndices4 = [matchingIndices4, i];  % Store the row index
    end
end

% Extract the filtered rows from LCresults
filteredRows4 = LCresults(matchingIndices4);

%Plot these to visualize

figure()
cmap = parula(length(filteredRows4));
worldmap([18 32],[-100 -80]) 
colorbar
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
hold on
for i = 1:length(filteredRows4)
    if ~isempty(filteredRows4(i).Xcut)
        plotm(filteredRows4(i).Ycut,filteredRows4(i).Xcut,'Color',cmap(i,:)); 
    end
end

hold off

%%THESE SHOULD NOT BE REMOVED


%LOOK FOR OUTLIERS

% Extract CutLength values from the structure
cutLengthValues = [LCresults.CutLength];
% Create a box-and-whisker plot
figure;
boxplot(cutLengthValues);
histogram(cutLengthValues)

% Identify outliers in the CutLength field based on the mean
outlierIndices = isoutlier([LCresults.CutLength], 'mean');

% Create a subset of the outliers for review
outlierSubset = LCresults(outlierIndices);

%PLOT outliers
figure()

worldmap([18 32],[-100 -80]) 

geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
hold on
for i = 15:15   %length(outlierSubset)
    if ~isempty(outlierSubset(i).Xcut)
        plotm(outlierSubset(i).Ycut,outlierSubset(i).Xcut,'Color',cmap(i,:)); 
    end
end

hold off

%After examining outliers and some dates around, going to drop everything
%that is greater than 3,000,000 rather than all outliers

% Find indices where CutLength is greater than 3,000,000
highCutIndices = [LCresults.CutLength] > 3e6;

% Create a subset for review
highCutSubset = LCresults(highCutIndices);

% Keep only the rows that do not match the highCutIndices
LCresults = LCresults(~highCutIndices);




%% Step 2: Plot daily LC for study period (2001 - 2023)

% Figure 1: Map of LC colored by length - No smoothing - try colors by
% season, ramp by year

%Cyclic color scale assigning a color to each day of the year - maybe
%sample so we only have 1 a week

%% PLOT LOOP CURRENT CONTOURS WITH ADJUSTMENTS

% Extract CutLength and Date values from your results structure
lineLengths = [LCresults.CutLength];
contourDates = [LCresults.Date];

% --- Define special periods and colors ---
specialPeriods = [datetime(2005, 1, 2), datetime(2005, 4, 2);   % Q1 2005
                  datetime(2005, 7, 3), datetime(2005, 10, 1);  % Q3 2005
                  datetime(2016, 7, 3), datetime(2016, 10, 1);  % Q3 2016
                  datetime(2022, 7, 3), datetime(2022, 10, 1)]; % Q3 2022

% Define colors for special periods
color2005_Q1 = [247, 215, 255] / 255; % Light Lavender/Pink
color2005_Q3 = [251, 0, 116] / 255;   % Deep Pink
color2016_Q3 = [102, 0, 153] / 255;   % Deeper Purple
color2022_Q3 = [255, 0, 255] / 255;   % Fuchsia/Hot Pink

% Create the map figure
figure()
hold on

% Load and plot the coastline shapefile
if isfile('coastL1.shp')
    geoshow('coastL1.shp', 'FaceColor', '#D3D3D3', 'FaceAlpha', 1);
else
    warning('Coastline shapefile "coastL1.shp" not found.');
end

% Check for coral shapefile and plot
if isfile('coral_gom3.shp')
    % --- MODIFICATION: Set both FaceColor and EdgeColor to red ---
    geoshow('coral_gom3.shp', 'FaceColor', '#FF0000', 'EdgeColor', '#FF0000', 'FaceAlpha', 1);
    disp('Coral reef shapefile plotted successfully.');
else
    warning('Coral reef shapefile "coral_gom3.shp" not found. It will not be plotted.');
end

% Set up the colormap for the general contour plotting
numBins = 50;
cmap = parula(numBins);
minLength = 0;
maxLength = max(lineLengths);
colormap(cmap);
c = colorbar;
c.Label.String = 'Contour Length';
caxis([minLength maxLength]);

% Identify indices for each special period
is2005_Q1 = (contourDates >= specialPeriods(1,1) & contourDates <= specialPeriods(1,2));
is2005_Q3 = (contourDates >= specialPeriods(2,1) & contourDates <= specialPeriods(2,2));
is2016_Q3 = (contourDates >= specialPeriods(3,1) & contourDates <= specialPeriods(3,2));
is2022_Q3 = (contourDates >= specialPeriods(4,1) & contourDates <= specialPeriods(4,2));

% Get indices that sort the lines by length in descending order
[~, sortedIndices] = sort(lineLengths, 'descend');

% Plot all contours from longest to shortest
disp('Plotting all contours from longest to shortest...');
for i = sortedIndices
    colorIdx = round((lineLengths(i) - minLength) / (maxLength - minLength) * (size(cmap,1) - 1)) + 1;
    colorToUse = cmap(colorIdx, :);
    plot(LCresults(i).Xcut, LCresults(i).Ycut, 'Color', [colorToUse, 0.5]);
end

% Function to plot the longest and shortest lines for special periods
plot_extreme_lines = @(indices, color) ...
    plot(LCresults(indices).Xcut, LCresults(indices).Ycut, 'Color', color, 'LineWidth', 2);

% Process each special period to highlight the extreme contours
disp('Highlighting extreme contours for special periods...');
specialPeriodIndices = {is2005_Q1, is2005_Q3, is2016_Q3, is2022_Q3};
specialPeriodColors = {color2005_Q1, color2005_Q3, color2016_Q3, color2022_Q3};

for j = 1:length(specialPeriodIndices)
    idx = find(specialPeriodIndices{j});
    if ~isempty(idx)
        [~, maxIdx] = max(lineLengths(idx));
        [~, minIdx] = min(lineLengths(idx));
        
        plot_extreme_lines(idx(maxIdx), specialPeriodColors{j});
        plot_extreme_lines(idx(minIdx), specialPeriodColors{j});
    end
end

% Add an updated legend for the special periods
h1 = plot(nan, nan, 'Color', color2005_Q1, 'LineWidth', 2);
h2 = plot(nan, nan, 'Color', color2005_Q3, 'LineWidth', 2);
h3 = plot(nan, nan, 'Color', color2016_Q3, 'LineWidth', 2);
h4 = plot(nan, nan, 'Color', color2022_Q3, 'LineWidth', 2);
legend([h1, h2, h3, h4], {'Q1 2005', 'Q3 2005', 'Q3 2016', 'Q3 2022'}, 'Location', 'southoutside')

% Set map limits and finalize plot
xlim([-100 -80])
ylim([18 32])
title('Loop Current Contours with Special Period Highlights')
xlabel('Longitude')
ylabel('Latitude')
hold off
disp('Plotting complete.')

f=gcf
exportgraphics(f, 'Lcmap_nodiseaseyear_NEW.png', 'Resolution',500)

%OLD CODE WITHOUT HYPOXIA PERIOD
% %% PLOT DISEASE TIME PERIODS WITH COLOR-MAPPED LINE LENGTHS
% % Extract CutLength values
% lineLengths = [LCresults.CutLength];
% 
% % Extract Dates (Assuming LCresults.Date exists as datetime)
% contourDates = [LCresults.Date];
% 
% % Define colormap based on line lengths (50 colors to match histogram)
% numBins = 50;
% cmap = parula(numBins);
% 
% % Define time periods for special coloring
% specialPeriods = [datetime(2005, 1, 2), datetime(2005, 4, 2); % Jan-Feb 2005
%                   datetime(2005, 7, 3), datetime(2005, 10, 1); % Jul-Oct 2005
%                   datetime(2022, 7, 3), datetime(2022, 10, 1)]; % 2022 period
% 
% % Create world map
% figure()
% hold on
% 
% % Load and plot the coastline shapefile
% geoshow('coastL1.shp', 'FaceColor', '#D3D3D3', 'FaceAlpha', 1)
% 
% % Load and plot the additional coral shapefile with the requested color -
% % changed to red
% geoshow('coral_gom3.shp', 'FaceColor', '#FF0000', 'FaceAlpha', 1, 'EdgeColor', '#F0A9DD')
% 
% % Normalize CutLength values for colormap indexing
% minLength = 0;
% maxLength = max(lineLengths);
% colormap(cmap)
% c = colorbar;
% c.Label.String = 'Line Length';
% caxis([minLength maxLength])
% 
% % Identify special period indices
% is2005_1 = (contourDates >= specialPeriods(1,1) & contourDates <= specialPeriods(1,2));
% is2005_2 = (contourDates >= specialPeriods(2,1) & contourDates <= specialPeriods(2,2));
% is2022 = (contourDates >= specialPeriods(3,1) & contourDates <= specialPeriods(3,2));
% 
% % Define colors for longest and shortest lines in each period
% color2005_1 = [247, 215, 255] / 255; % Jan-Feb 2005
% color2005_2 = [251, 0, 116] / 255;   % Jul-Oct 2005
% color2022 = [214, 0, 255] / 255;     % 2022
% 
% alphaValue = 1; % Opacity
% 
% % **First, plot all lines (including disease periods)**
% for i = 1:length(LCresults)
%     colorIdx = round((lineLengths(i) - minLength) / (maxLength - minLength) * (size(cmap,1) - 1)) + 1;
%     colorToUse = cmap(colorIdx, :);
%     plot(LCresults(i).Xcut, LCresults(i).Ycut, 'Color', [colorToUse, 0.5]);
% end
% 
% % Function to plot longest and shortest lines
% plot_extreme_lines = @(indices, color) ...
%     plot(LCresults(indices).Xcut, LCresults(indices).Ycut, 'Color', color, 'LineWidth', 2);
% 
% % Process each disease period
% diseasePeriods = {is2005_1, is2005_2, is2022};
% diseaseColors = {color2005_1, color2005_2, color2022};
% meanLengths = zeros(1,3);
% 
% for j = 1:3
%     idx = find(diseasePeriods{j});
%     if ~isempty(idx)
%         [~, maxIdx] = max(lineLengths(idx));
%         [~, minIdx] = min(lineLengths(idx));
%         plot_extreme_lines(idx(maxIdx), diseaseColors{j}); % Longest
%         plot_extreme_lines(idx(minIdx), diseaseColors{j}); % Shortest
%         meanLengths(j) = mean(lineLengths(idx)); % Compute mean
%     end
% end
% 
% Add legend for disease periods
% h1 = plot(nan, nan, 'Color', color2005_1, 'LineWidth', 2);
% h2 = plot(nan, nan, 'Color', color2005_2, 'LineWidth', 2);
% h3 = plot(nan, nan, 'Color', color2022, 'LineWidth', 2);
% legend([h1, h2, h3], {'Jan-Feb 2005', 'Jul-Oct 2005', '2022'}, 'Location', 'southoutside')
% 
% 
% xlim([-100 -80])
% ylim([18 32])
% hold off
% 
% f=gcf
% exportgraphics(f, 'Lcmap_nodiseaseyear.png', 'Resolution',500)
% 
% %% Calculate Medians for Histogram
% 
% % Extract CutLength values and Dates
% lineLengths = [LCresults.CutLength];
% contourDates = [LCresults.Date];
% 
% % Define disease periods
% start2005_1 = datetime(2005, 1, 2); end2005_1 = datetime(2005, 4, 2);
% start2005_2 = datetime(2005, 7, 3); end2005_2 = datetime(2005, 10, 1);
% start2022 = datetime(2022, 7, 3); end2022 = datetime(2022, 10, 1);
% 
% % Initialize median values
% median2005_1 = NaN;
% median2005_2 = NaN;
% median2022 = NaN;
% medianAll = NaN;  % Median for all dates
% 
% % Calculate overall median
% if ~isempty(lineLengths)
%     medianAll = median(lineLengths);
% end
% 
% % Find and calculate medians using explicit if statements
% if any(contourDates >= start2005_1 & contourDates <= end2005_1)
%     median2005_1 = median(lineLengths(contourDates >= start2005_1 & contourDates <= end2005_1));
% end
% 
% if any(contourDates >= start2005_2 & contourDates <= end2005_2)
%     median2005_2 = median(lineLengths(contourDates >= start2005_2 & contourDates <= end2005_2));
% end
% 
% if any(contourDates >= start2022 & contourDates <= end2022)
%     median2022 = median(lineLengths(contourDates >= start2022 & contourDates <= end2022));
% end
% 
% % Display results
% disp('Median Line Lengths:')
% disp(['All Dates: ', num2str(medianAll)])
% disp(['Q1 2005: ', num2str(median2005_1)])
% disp(['Q3 2005: ', num2str(median2005_2)])
% disp(['Q3 2022: ', num2str(median2022)])

%% DATA PREPARATION AND MEDIAN CALCULATION
% This section calculates the medians needed for the plot.

% Extract data from your results structure
cutLengths = [LCresults.CutLength];
contourDates = [LCresults.Date];

% Define the special time periods to match the map
specialPeriods = [datetime(2005, 1, 2), datetime(2005, 4, 2);   % Q1 2005
                  datetime(2005, 7, 3), datetime(2005, 10, 1);  % Q3 2005
                  datetime(2016, 7, 3), datetime(2016, 10, 1);  % Q3 2016 (New)
                  datetime(2022, 7, 3), datetime(2022, 10, 1)]; % Q3 2022

% Create logical indices to filter data for each period
is2005_Q1 = (contourDates >= specialPeriods(1,1) & contourDates <= specialPeriods(1,2));
is2005_Q3 = (contourDates >= specialPeriods(2,1) & contourDates <= specialPeriods(2,2));
is2016_Q3 = (contourDates >= specialPeriods(3,1) & contourDates <= specialPeriods(3,2));
is2022_Q3 = (contourDates >= specialPeriods(4,1) & contourDates <= specialPeriods(4,2));

% Calculate the median length for each period and overall
median2005_Q1 = median(cutLengths(is2005_Q1));
median2005_Q3 = median(cutLengths(is2005_Q3));
median2016_Q3 = median(cutLengths(is2016_Q3)); % New median
median2022_Q3 = median(cutLengths(is2022_Q3));
medianAll = median(cutLengths); % Overall median

%% CUTLENGTH HISTOGRAM PLOTTING
% This section creates the histogram figure.

% Define bin edges and counts
numBins = 50;
edges = linspace(min(cutLengths), max(cutLengths), numBins + 1);
[N, ~] = histcounts(cutLengths, edges, 'Normalization', 'probability');

% Convert bin heights to percentages
N_percent = N * 100; 

% Create a colormap based on bin centers (length values)
binCenters = (edges(1:end-1) + edges(2:end)) / 2;
cmap = parula(numBins);

% Create the bar plot with colored bars
hist_fig = figure;
barHandles = bar(binCenters, N_percent, 'FaceColor', 'flat', 'EdgeColor', [0.5 0.5 0.5], 'BarWidth', 1);

% Assign colors to each bar based on its corresponding length value
for i = 1:numBins
    barHandles.CData(i, :) = cmap(i, :);
end

% Add colorbar and labels
c = colorbar;
c.Label.String = 'Contour Length (km)';
caxis([min(cutLengths) max(cutLengths)]);

% Adjust y-axis to show percentages
ylabel('Percentage of Total (%)');
xlabel('LC Length (km)');
title('Distribution of Loop Current Contour Lengths');

% --- MODIFICATION: Define colors to match the map plot ---
color2005_Q1 = [247, 215, 255] / 255; % Light Lavender/Pink
color2005_Q3 = [251, 0, 116] / 255;   % Deep Pink
color2016_Q3 = [102, 0, 153] / 255;   % Deeper Purple
color2022_Q3 = [255, 0, 255] / 255;   % Fuchsia/Hot Pink
colorOverall = [0, 0, 0];            % Black for overall median

% Plot vertical lines for medians
hold on;
p1 = xline(median2005_Q1, '--', 'Color', color2005_Q1, 'LineWidth', 2);
p2 = xline(median2005_Q3, '--', 'Color', color2005_Q3, 'LineWidth', 2);
p3 = xline(median2016_Q3, '--', 'Color', color2016_Q3, 'LineWidth', 2); % New period line
p4 = xline(median2022_Q3, '--', 'Color', color2022_Q3, 'LineWidth', 2);
p_all = xline(medianAll, '-', 'Color', colorOverall, 'LineWidth', 2);
hold off;

% --- MODIFICATION: Add a legend to identify the median lines ---
legend([p_all, p1, p2, p3, p4], ...
    {sprintf('Overall Median: %.1f km', medianAll), ...
     sprintf('Q1 2005 Median: %.1f km', median2005_Q1), ...
     sprintf('Q3 2005 Median: %.1f km', median2005_Q3), ...
     sprintf('Q3 2016 Median: %.1f km', median2016_Q3), ...
     sprintf('Q3 2022 Median: %.1f km', median2022_Q3)}, ...
    'Location', 'northeast')

f=gcf
exportgraphics(f, 'Lchist_NEW.pdf', 'Resolution',800)


% 
% OLD
% %% CUTLENGTH HISTOGRAM
% % Extract CutLength values
% cutLengths = [LCresults.CutLength];
% 
% % Define bin edges and counts
% numBins = 50; % Adjust as needed
% edges = linspace(min(cutLengths), max(cutLengths), numBins + 1); % Bin edges
% [N, edges] = histcounts(cutLengths, edges, 'Normalization', 'probability'); % Histogram counts
% 
% % Convert bin heights to percentages
% N = N * 100; 
% 
% % Create a colormap based on bin centers
% binCenters = (edges(1:end-1) + edges(2:end)) / 2; % Midpoints of bins
% colormap(parula); % Apply colormap
% cmap = parula(numBins); % Get colors from colormap
% 
% % Create a bar plot with colored bars and black borders
% hist = figure;
% barHandles = bar(binCenters, N, 'FaceColor', 'flat', 'EdgeColor', [0.5 0.5 0.5], 'BarWidth', 1); % Black borders
% 
% % Assign colors to each bar based on bin center value
% for i = 1:numBins
%     barHandles.CData(i, :) = cmap(i, :);
% end
% 
% % Add colorbar and label
% c = colorbar;
% c.Label.String = 'Cut Length';
% caxis([min(cutLengths) max(cutLengths)]); % Scale colorbar
% 
% % Adjust y-axis to show percentages
% yt = yticks;
% yticklabels(string(yt)); % Ensure y-axis shows percentages
% 
% % Add labels and title
% xlabel('LC Length (km)');
% ylabel('Percentage of Total (%)');
% title('Distribution of LC Lengths');
% 
% % Define disease period colors (matching the map)
% color2005_1 = [247, 215, 255] / 255; % Jan-Apr 2005
% color2005_2 = [251, 0, 116] / 255;   % Jul-Oct 2005
% color2022 = [214, 0, 255] / 255;     % 2022
% colorOverall = [0, 0, 0];            % Black for overall median
% 
% 
% % Plot vertical lines for medians
% hold on
% xline(median2005_1, '--', 'Color', color2005_1, 'LineWidth', 2);
% xline(median2005_2, '--', 'Color', color2005_2, 'LineWidth', 2);
% xline(median2022, '--', 'Color', color2022, 'LineWidth', 2);
% xline(medianAll, '-', 'Color', colorOverall, 'LineWidth', 2);
% 
% hold off



exportgraphics(hist, 'Lchist.pdf', 'Resolution',1200)


%% Calculation from histogram
% Define threshold
retracted = 1e6; % 1,000,000

% Find bins where LC Length < 1,000,000
binsBelowThreshold = binCenters < retracted;

% Calculate cumulative percentage
cumulativePercent = sum(N(binsBelowThreshold));

% Display result
disp(['Cumulative percentage of LC lengths below 1,000,000: ', num2str(cumulativePercent), '%']);

% Define threshold
extended = 2e6; % 1,000,000

% Find bins where LC Length > 2,000,000
binsabove = binCenters > extended;

% Calculate cumulative percentage
cumulativePercent2 = sum(N(binsabove));

% Display result
disp(['Cumulative percentage of LC lengths below 1,000,000: ', num2str(cumulativePercent2), '%']);



%% Figure 4: LC lengths over study period

% First, create a datetime field in LCresults
for i = 1:numel(LCresults)
    LCresults(i).Date = datetime(LCresults(i).Year, LCresults(i).Month, LCresults(i).Day);
end

% Extract data
LClengths = [LCresults.CutLength];
LClengths(LClengths == 0) = NaN; % Replace zeros with NaN to avoid plotting them
dates = [LCresults.Date];

% **Sort data by date**
[dates, sortIdx] = sort(dates); % Sort dates in ascending order
LClengths = LClengths(sortIdx); % Reorder lengths accordingly

% Define shaded time periods for even years from 2002 to 2022
evenYears = 2002:2:2022; % Generate even years
shadePeriods = [datetime(evenYears,1,1)', datetime(evenYears,12,31)']; % Create start and end dates

% **Create Figure**
LClengthFig = figure;
hold on;

% **Dynamically set y-axis limits**
yLimits = [min(LClengths, [], 'omitnan'), max(LClengths, [], 'omitnan')];
if yLimits(1) == yLimits(2) % In case all values are the same
    yLimits = yLimits + [-1, 1]; % Expand slightly to avoid a flat line
end

% **Shade the even-year periods**
for i = 1:size(shadePeriods,1)
    xPatch = [shadePeriods(i,1) shadePeriods(i,2) shadePeriods(i,2) shadePeriods(i,1)];
    yPatch = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];
    patch(xPatch, yPatch, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

% **Plot the LC lengths over time (New Color, No Markers)**
plot(dates, LClengths, 'Color', [53/255, 42/255, 135/255], 'LineWidth', 1.5); % "#352a87"

% **Format x-axis to show years centered at July 1st**
years = year(dates);  % Extract years
uniqueYears = unique(years);  % Find unique years
xticks(datetime(uniqueYears, 7, 1));  % Set ticks at July 1st each year
datetick('x', 'yyyy', 'keepticks'); % Format x-axis as years only

% **Labels, Title, and Grid Adjustments**
xlabel('Year');
ylabel('LC Lengths');
title('LC Lengths Over Time');

% **Remove Grid Crosshairs**
ax = gca;
ax.XGrid = 'off'; % No vertical grid lines
ax.YGrid = 'off'; % No horizontal grid lines

% **Final Adjustments**
xlim([min(dates) max(dates)]); % Ensure limits fit data
ylim(yLimits); % Keep y-limits consistent
hold off;

%%INCLUDES MISSING DATES WITH RED DOT

% First, create a datetime field in LCresults
for i = 1:numel(LCresults)
    LCresults(i).Date = datetime(LCresults(i).Year, LCresults(i).Month, LCresults(i).Day);
end

% Extract data
LClengths = [LCresults.CutLength];
LClengths(LClengths == 0) = NaN; % Replace zeros with NaN to avoid plotting them
dates = [LCresults.Date];

% **Sort data by date**
[dates, sortIdx] = sort(dates); % Sort dates in ascending order
LClengths = LClengths(sortIdx); % Reorder lengths accordingly

% **Find missing dates**
allDates = (min(dates):max(dates))'; % Generate a complete date range
missingDates = setdiff(allDates, dates); % Identify missing dates

% **Find LC lengths for missing dates (use previous day's value)**
missingLengths = nan(size(missingDates)); % Preallocate
for i = 1:numel(missingDates)
    prevIdx = find(dates < missingDates(i), 1, 'last'); % Find last available date
    if ~isempty(prevIdx)
        missingLengths(i) = LClengths(prevIdx); % Use previous day's LC length
    end
end

% Define shaded time periods for even years from 2002 to 2022
evenYears = 2002:2:2022; % Generate even years
shadePeriods = [datetime(evenYears,1,1)', datetime(evenYears,12,31)']; % Create start and end dates

% **Create Figure**
LClengthFig = figure;
hold on;

% **Dynamically set y-axis limits**
yLimits = [min(LClengths, [], 'omitnan'), max(LClengths, [], 'omitnan')];
if yLimits(1) == yLimits(2) % In case all values are the same
    yLimits = yLimits + [-1, 1]; % Expand slightly to avoid a flat line
end

% **Shade the even-year periods**
for i = 1:size(shadePeriods,1)
    xPatch = [shadePeriods(i,1) shadePeriods(i,2) shadePeriods(i,2) shadePeriods(i,1)];
    yPatch = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];
    patch(xPatch, yPatch, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

% **Plot the LC lengths over time (New Color, No Markers)**
plot(dates, LClengths, 'Color', [53/255, 42/255, 135/255], 'LineWidth', 1.5); % "#352a87"

% **Plot missing dates as smaller red dots at previous day's LC length**
plot(missingDates, missingLengths, 'wo', 'MarkerSize', 2, 'MarkerFaceColor', 'w');

% **Format x-axis to show years centered at July 1st**
years = year(dates);  % Extract years
uniqueYears = unique(years);  % Find unique years
xticks(datetime(uniqueYears, 7, 1));  % Set ticks at July 1st each year
datetick('x', 'yyyy', 'keepticks'); % Format x-axis as years only

% **Labels, Title, and Grid Adjustments**
xlabel('Year');
ylabel('LC Lengths');
title('LC Lengths Over Time (Missing Dates Indicated)');

% **Remove Grid Crosshairs**
ax = gca;
ax.XGrid = 'off'; % No vertical grid lines
ax.YGrid = 'off'; % No horizontal grid lines

% **Final Adjustments**
xlim([min(dates) max(dates)]); % Ensure limits fit data
ylim(yLimits); % Keep y-limits consistent
hold off;


% Set figure size (12 inches wide, 4 inches tall) before exporting
set(LClengthFig, 'Units', 'inches', 'Position', [1, 1, 15, 5]); % Position in inches
set(LClengthFig, 'PaperUnits', 'inches', 'PaperSize', [15, 5]); % Paper size for export

% Export as high-quality PDF (12" x 4" at 300 dpi)
print(LClengthFig, 'LC_Lengths_Over_Time.pdf', '-dpdf', '-r300', '-painters');


% **Export as High-Quality PDF (12 inches x 4 inches)**
set(LClengthFig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 24 4]); % Set figure size
print(LClengthFig, 'LC_Lengths_Over_Time.pdf', '-dpdf', '-r300'); % Export as 300 dpi PDF

%COLORS FROM DAN 

lineLengths = [LCresults.CutLength];
% Extract Dates (Assuming LCresults.Date exists as datetime)
contourDates = [LCresults.Date];
% Define colormap based on line lengths (50 colors to match histogram)
numBins = 50;
cmap = parula(numBins);
% Define time periods for special magenta coloring
specialPeriods = [datetime(2005, 1, 1), datetime(2005, 2, 28); % Jan-Feb 2005 period
                  datetime(2005, 7, 3), datetime(2005, 10, 1); % Jul-Oct 2005 period
                  datetime(2022, 7, 3), datetime(2022, 10, 1)]; % 2022 period
% Create world map
figure()
hold on

% Load and plot the coastline shapefile
geoshow('coastL1.shp', 'FaceColor', '#D3D3D3', 'FaceAlpha', 1)

% Load and plot the additional coral shapefile with the requested color
geoshow('coral_gom3.shp', 'FaceColor', '#F0A9DD', 'FaceAlpha', 1, 'EdgeColor', '#F0A9DD')


% Normalize CutLength values for colormap indexing
minLength = 0;
maxLength = max(lineLengths); % Get the longest line length
colormap(cmap) % Set colormap to match histogram
c = colorbar; % Add colorbar
c.Label.String = 'Line Length'; % Label colorbar
caxis([minLength maxLength]) % Scale colorbar from 0 to max length
% Identify special period indices
is2005_1 = (contourDates >= specialPeriods(1,1) & contourDates <= specialPeriods(1,2)); % Jan-Feb 2005 period
is2005_2 = (contourDates >= specialPeriods(2,1) & contourDates <= specialPeriods(2,2)); % Jul-Oct 2005 period
is2022 = (contourDates >= specialPeriods(3,1) & contourDates <= specialPeriods(3,2)); % 2022 period
isSpecial = is2005_1 | is2005_2 | is2022;
% Count special and non-special lines
num2005_1 = sum(is2005_1);
num2005_2 = sum(is2005_2);
num2022 = sum(is2022);
numRegular = sum(~isSpecial);
% Define new colors
color2005_1 = [247, 215, 255] / 255; % Jan-Feb 2005 (#d600ff)
color2005_2 = [251, 0, 116] / 255; % Jul-Oct 2005 (#fb00be)
color2022 = [214, 0, 255] / 255; % 2022 (#f7d7ff)
alphaValue = 1; % 50% opacity
% **First, plot non-special lines (regular colormap)**
for i = 1:length(LCresults)
    if ~isSpecial(i) % Plot only if NOT in special period
        % Normalize color index for colormap
        colorIdx = round((lineLengths(i) - minLength) / (maxLength - minLength) * (size(cmap,1) - 1)) + 1;
        colorToUse = cmap(colorIdx, :);
        plot(LCresults(i).Xcut, LCresults(i).Ycut, 'Color', [colorToUse, 0.5]);
    end
end
% **Then, plot special lines (semi-transparent magenta shades)**
for i = 1:length(LCresults)
    if isSpecial(i)
        if is2005_1(i)
            colorToUse = color2005_1; % Jan-Feb 2005 period
        elseif is2005_2(i)
            colorToUse = color2005_2; % Jul-Oct 2005 period
        else
            colorToUse = color2022; % 2022 period
        end
        
        % Using `patch` for transparency effect
        x = LCresults(i).Xcut;
        y = LCresults(i).Ycut;
        patch([x nan], [y nan], colorToUse, 'EdgeColor', colorToUse, 'FaceAlpha', alphaValue, 'EdgeAlpha', alphaValue);
    end
end
% **Add Legend for Special Periods**
h1 = plot(nan, nan, 'Color', color2005_1, 'LineWidth', 2); % Jan-Feb 2005 period
h2 = plot(nan, nan, 'Color', color2005_2, 'LineWidth', 2); % Jul-Oct 2005 period
h3 = plot(nan, nan, 'Color', color2022, 'LineWidth', 2); % 2022 period
legend([h1, h2, h3], {'Jan-Feb 2005 Disease Period', 'Jul-Oct 2005 Disease Period', '2022 Disease Period'}, 'Location', 'southoutside')
xlim([-100 -80])
ylim([18 32])
hold off

%OLD MAPS WITH ALL DISEASE PERIOD LINES

% Extract CutLength values
lineLengths = [LCresults.CutLength];

% Extract Dates (Assuming LCresults.Date exists as datetime)
contourDates = [LCresults.Date];

% Define colormap based on line lengths (50 colors to match histogram)
numBins = 50; 
cmap = parula(numBins); 

% Define time periods for special magenta coloring
specialPeriods = [datetime(2005, 1, 1), datetime(2005, 2, 28); % Jan-Feb 2005 period
                  datetime(2005, 7, 3), datetime(2005, 10, 1); % Jul-Oct 2005 period
                  datetime(2022, 7, 3), datetime(2022, 10, 1)]; % 2022 period

% Create world map
figure()
hold on

% Load and plot the coastline shapefile
geoshow('coastL1.shp', 'FaceColor', '#D3D3D3', 'FaceAlpha', 1)

% Load and plot the additional coral shapefile with the requested color
geoshow('coral_gom3.shp', 'FaceColor', '#F0A9DD', 'FaceAlpha', 1, 'EdgeColor', '#F0A9DD')

% Normalize CutLength values for colormap indexing
minLength = 0; 
maxLength = max(lineLengths); % Get the longest line length

colormap(cmap) % Set colormap to match histogram
c = colorbar; % Add colorbar
c.Label.String = 'Line Length'; % Label colorbar
caxis([minLength maxLength]) % Scale colorbar from 0 to max length

% Identify special period indices
is2005_1 = (contourDates >= specialPeriods(1,1) & contourDates <= specialPeriods(1,2)); % Jan-Feb 2005 period
is2005_2 = (contourDates >= specialPeriods(2,1) & contourDates <= specialPeriods(2,2)); % Jul-Oct 2005 period
is2022 = (contourDates >= specialPeriods(3,1) & contourDates <= specialPeriods(3,2)); % 2022 period
isSpecial = is2005_1 | is2005_2 | is2022;

% Count special and non-special lines
num2005_1 = sum(is2005_1);
num2005_2 = sum(is2005_2);
num2022 = sum(is2022);
numRegular = sum(~isSpecial);

% Define new colors
color2005_1 = [252, 137, 97] / 255; % Jan-Feb 2005 (#d600ff)
color2005_2 = [183, 55, 121] / 255; % Jul-Oct 2005 (#fb00be)
color2022 = [81, 18, 124] / 255; % 2022 (#f7d7ff)
alphaValue = 1; % 50% opacity

% **First, plot non-special lines (regular colormap)**
for i = 1:length(LCresults)
    if ~isSpecial(i) % Plot only if NOT in special period
        % Normalize color index for colormap
        colorIdx = round((lineLengths(i) - minLength) / (maxLength - minLength) * (size(cmap,1) - 1)) + 1;
        colorToUse = cmap(colorIdx, :);
        plot(LCresults(i).Xcut, LCresults(i).Ycut, 'Color', [colorToUse, 0.5]);
    end
end

% **Then, plot special lines (semi-transparent magenta shades)**
for i = 1:length(LCresults)
    if isSpecial(i)
        if is2005_1(i)
            colorToUse = color2005_1; % Jan-Feb 2005 period
        elseif is2005_2(i)
            colorToUse = color2005_2; % Jul-Oct 2005 period
        else
            colorToUse = color2022; % 2022 period
        end
        
        % Using `patch` for transparency effect
        x = LCresults(i).Xcut;
        y = LCresults(i).Ycut;
        patch([x nan], [y nan], colorToUse, 'EdgeColor', colorToUse, 'FaceAlpha', alphaValue, 'EdgeAlpha', alphaValue);
    end
end

% **Add Legend for Special Periods**
h1 = plot(nan, nan, 'Color', color2005_1, 'LineWidth', 2); % Jan-Feb 2005 period
h2 = plot(nan, nan, 'Color', color2005_2, 'LineWidth', 2); % Jul-Oct 2005 period
h3 = plot(nan, nan, 'Color', color2022, 'LineWidth', 2); % 2022 period
legend([h1, h2, h3], {'Jan-Feb 2005 Disease Period', 'Jul-Oct 2005 Disease Period', '2022 Disease Period'}, 'Location', 'southoutside')

xlim([-100 -80])
ylim([18 32])
hold off

%%TESTING WITH OPACITY
%% PLOT DISEASE TIME PERIODS WITH COLOR-MAPPED LINE LENGTHS

% Extract CutLength values and corresponding Dates
lineLengths = [LCresults.CutLength];
contourDates = [LCresults.Date];

% Define colormap based on line lengths (50 colors to match histogram)
numBins = 50; 
cmap = parula(numBins); 

% Define time periods for disease (fully opaque, thick lines)
diseasePeriods = [datetime(2005, 1, 1), datetime(2005, 2, 28); % Jan-Feb 2005
                  datetime(2005, 7, 3), datetime(2005, 10, 1); % Jul-Oct 2005
                  datetime(2022, 7, 3), datetime(2022, 10, 1)]; % 2022

% Create world map
figure()
hold on

% Load and plot the coastline shapefile
geoshow('coastL1.shp', 'FaceColor', '#D3D3D3', 'FaceAlpha', 1)

% Load and plot the additional coral shapefile with the requested color
geoshow('coral_gom3.shp', 'FaceColor', '#F0A9DD', 'FaceAlpha', 1, 'EdgeColor', '#F0A9DD')

% Normalize CutLength values for colormap indexing
minLength = 0; 
maxLength = max(lineLengths);

colormap(cmap) % Set colormap to match histogram
c = colorbar; % Add colorbar
c.Label.String = 'Line Length'; % Label colorbar
caxis([minLength maxLength]) % Scale colorbar from 0 to max length

% Identify disease period indices
isDisease = false(size(contourDates));
for i = 1:size(diseasePeriods,1)
    isDisease = isDisease | (contourDates >= diseasePeriods(i,1) & contourDates <= diseasePeriods(i,2));
end

% Set opacity and line width values
opacityRegular = 0.1;  % 30% opacity for non-disease times
opacityDisease = 1.0;  % 100% opacity for disease periods
lineWidthRegular = 0.5;  % Normal line width
lineWidthDisease = 0.5;  % Double thickness for disease periods

% **Sort indices by decreasing line length (longest first)**
[~, sortIdx] = sort(lineLengths, 'descend');

% **Plot all lines in order of longest to shortest**
for i = sortIdx
    % Normalize color index for colormap
    colorIdx = round((lineLengths(i) - minLength) / (maxLength - minLength) * (size(cmap,1) - 1)) + 1;
    colorToUse = cmap(colorIdx, :);
    
    % Set opacity and line width based on disease period
    if isDisease(i)
        alphaValue = opacityDisease; % Fully opaque for disease periods
        lineWidth = lineWidthDisease; % Thick lines
    else
        alphaValue = opacityRegular; % 30% opacity for non-disease periods
        lineWidth = lineWidthRegular; % Normal thickness
    end
    
    % Plot the line
    plot(LCresults(i).Xcut, LCresults(i).Ycut, 'Color', [colorToUse alphaValue], 'LineWidth', lineWidth);
end

% **Add Legend**
h1 = plot(nan, nan, 'k', 'LineWidth', lineWidthDisease); % Thick black line for disease periods
legend(h1, {'Disease Periods (Thicker & Fully Opaque)'}, 'Location', 'southoutside')

xlim([-100 -80])
ylim([18 32])
hold off





f=gcf
exportgraphics(f, 'Lcmap_300.png', 'Resolution',300)
print('DANOUTEXAMPLE','-dpdf','-vector')

print(LClengthFig, 'LC_Lengths_Over_Time.pdf', '-dpdf', '-r300', '-painters');


%%OTher
% Convert structure to table
LCresultsTable = struct2table(LCresults);

% Sort table based on CutLength field in descending order
LCresultsTable = sortrows(LCresultsTable, 'CutLength', 'descend');

% Convert back to structure
LCresults = table2struct(LCresultsTable);



% No smoothing
figure()
%cmap = parula(length(LCresults));
worldmap([18 32],[-100 -80]) 
colorbar
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
hold on
for i = 8270:8283
    if ~isempty(LCresults(i).Xcut)
        plotm(LCresults(i).Ycut,LCresults(i).Xcut,'Color',cmap(i,:)); 
    end
end

hold off

%Other stuff for normalization
figure;
histogram([LCresults.CutLength])

test = rescale([LCresults.CutLength],0,5);

for i = 1:length(LCresults)
    LCresults(i).norm = (LCresults(i).CutLength - min([LCresults.CutLength]))/(max([LCresults.CutLength])-min([LCresults.CutLength]))*5;
end

figure;
histogram(test,[0:.1:5],'Normalization','percentage')
