clear all
addpath(genpath(pwd))

% Import Global Centroid Moment Tensor Catalogue ('gCMT'), one of the 
% most commonly used global earthquake catalogues
% The provided magnitudes are Moment Magnitudes 'Mw'
catname = 'dat/jan76_dec10_plusMonthlies_jan11_dec14_lat-90_90_lon-360_360.mat';
cat     = read_GMT_cat(catname);
cat

% Plot simple seismicity map
clf; hold on; grid on; box on;
plot(cat.lon,cat.lat,'.k')

% Plot seismicity map with colour proporational to hypocentral depth
clf; hold on; grid on; box on;
scatter(cat.lon,cat.lat,20,cat.depth,'filled')


% Plot seismicity map in 3D 
clf; hold on; grid on; box on;
scatter3(cat.lon,cat.lat,cat.depth,4*cat.magnitude,cat.depth,'filled')

set(gca,'zdir','reverse','zlim',[0 1000])
cb = colorbar;
cb.Label.String = 'Hypocentral Depth [km]';
set(gca,'view',[-20 45])

% Zoom into Japan region
set(gca,'xlim',[125 155],'ylim',[25 55])


% Plot time versus magnitude
m = cat.magnitude;
mkSize = 1+(m-min(m)).^3;
clf; hold on; grid on; box on;
scatter(cat.t0,cat.magnitude,mkSize,cat.magnitude,'filled','markerEdgeColor',[.4 .4 .4])
cb = colorbar;
cb.Label.String = 'Moment Magnitude';
ylabel('Moment Magnitude')

% Repeat, but with colour proportional to depth
clf; hold on; grid on; box on;
scatter(cat.t0,cat.magnitude,mkSize,cat.depth,'filled','markerEdgeColor',[.4 .4 .4])
cb = colorbar;
cb.Label.String = 'Hypocentral Depth [km]';
ylabel('Moment Magnitude')










