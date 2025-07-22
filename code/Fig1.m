
%Add bathy
lonb = ncread('Hycom_depth.nc', 'lon');
latb = ncread('Hycom_depth.nc','lat');
bathy = ncread('Hycom_depth.nc','bathy');

latbv = latb(1,:);
lonbv = lonb(:,1);


imagesc(lonbv, latbv, bathy');set(gca,'YDir','normal');
axis xy; 


%Make map

[LatGrid, LonGrid] = meshgrid(latbv, lonbv); % Flipped to match bathy orientation

f1 = figure();
hold on

worldmap([18 32], [-100 -80]);

% Correctly oriented bathymetry
geoshow(LatGrid, LonGrid, bathy, 'DisplayType', 'texturemap');

% Custom grayscale colormap
nColors = 256;
% minGray = 1;
% maxGray = 0.1;
% customGray = linspace(minGray, maxGray, nColors)';
% cmap = repmat(customGray, 1, 3);
%cmap = cmocean('deep')
colormap(sky);
colorbar;

% Optional: color limits
caxis([min(bathy(:)) max(bathy(:))]);

geoshow('coastL1.shp', 'FaceColor', '#D3D3D3', 'FaceAlpha', 1);

g(1) = geoshow('coral_gom3.shp', 'FaceColor', '#FF0000', ...
               'FaceAlpha', 1, 'EdgeColor', '#FF0000');



legend([g(1)], {'Reef areas'}, 'Location', 'eastoutside');


%scaleruler('MajorTick',0:100:200)

northarrow("Latitude",22,"Longitude",-99,"LineWidth",1.5);

scaleruler on
setm(handlem('scaleruler1'), 'XLoc',0,'YLoc',1,'MajorTick',0:200:400,'RulerStyle','patches')




exportgraphics(f1,'fig_1_base_red.png','Resolution',1200) 