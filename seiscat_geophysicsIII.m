% Generic catalog overview; should work well with all catalogs imported 
% with read_eq_catalog.m
% Modified version of seismo/matlab/projects/seismicity/analyse_eq_catalog.m 
% for BULGG xyz coordinate system
%
% Need to change some paths to generalmodel/ dir
%
% Note: some clean-ups / updates in share/seismicity/analyse_eq_catalog_EH.m

clear all
addpath(genpath(pwd)); % adds all subfolders
addpath(genpath('../../../../geophysical-characterization/BULG-3D-model/'))
addpath(genpath('../../../../geophysical-characterization/BULG-3D-model/01-Data/'))
addpath(genpath('../../../tools/coordinateSystems/'))



fprintf(1,'check and document time of all catalogs.\n')
%also incorporate latest GES catalog, from ~210518
%pause



% Outstanding Tasks
% -----------------
% x need to convert template locations to xyz
% x TK's TM cats don't have locations; assume same location as template loc
% x BD's locations
% x add SED cat events
% x make flip plots to see overlap among cats
% x need consistent depth & altittude definitions: 
%     alt = altitude above sea level in m
%     dep = depth below sea level in m
%     z   = depth in local CS with z = CS origin depth at tunnel level (?) in [m]
% . need consistent CS definitions

% Questions 
% ---------
% . what's resolution of TK templates depth of 1.0km?
% .


% OVERVIEW ================================ 
% . Load catalog(s)
% . Select reference point
% . Select sub-catalog
% . Plot overview
% 
% Older stuff
% . Plot seismicity maps
% . Plot seismicity profiles
% . Plot seismicity rates
% . Stem plot of m vs t
% . Plot cluster evolution with increasing time windows
% . Pairplot of everything with everything ...
% . Plot FMDs
%   . in time windows
%   . in depth windows
% . Plot seismicity rates
% . Save video avi-file from VideoWriter
% =========================================





%% Load catalogs



% GES catalogs
% ------------
% Load and concatenate all GES catalogs
% Where is csv file for largest events described in pptx slides?
catGES1 = load_GES_catalogs;
catGES2 = read_eqcat_Bedretto('dat/ges/GES_largest_events_from_ppt_q.csv','tmpGES');
catGES  = merge_GES_cats(catGES1,catGES2);
%catGES = load_single_GES_catalog;



% TKTM catalogs
% -------------
% Load TKs TM detection catalog
catFullName = 'dat/tk/detection-catalog-combo_Mag.dat';
catTMdet    = read_eqcat_Bedretto(catFullName,'TKTMdet_2105');

% Load TKs TM detection catalog (only 5 templates)
catFullName = 'dat/tk/templates_mam.dat';
catTMtpl    = read_eqcat_Bedretto(catFullName,'TKTMtpl_2105');

% Write coordinates of TMtmpl into TMdet cat, and concatenate
catTM = combine_TM_cats(catTMdet,catTMtpl);


% SED catalogs
% ------------
catFullName = 'dat/tk/Bedretto_SEDevents.csv';
catSED    = read_eqcat_Bedretto(catFullName,'SED_2105');



% Transfer lat/lon/alt coordinates into local BULGG x/y/z CS
[catSED.y,  catSED.x,  catSED.z]   = WGS84_to_local_Bedretto_coords(catSED.lat,   catSED.lon,   catSED.alt);
[catTM.y,   catTM.x,   catTM.z]    = WGS84_to_local_Bedretto_coords(catTM.lat,    catTM.lon,    catTM.alt);
[catTMtpl.y,catTMtpl.x,catTMtpl.z] = WGS84_to_local_Bedretto_coords(catTMtpl.lat, catTMtpl.lon, catTMtpl.alt);






%% Plot 3d seismicity with tunnel & boreholes (copy from BULG-3D-model/)
% close all
hf = plot_BULG_3D_model_tunnel_and_boreholes_210515; hold on; box on;
set(gca,'xlim',[-500 800],'ylim',[-1600 200],'zlim',[180 1500])

h0 = scatter3(catSED.x,   catSED.y,   catSED.alt,   100, 'filled','markerFaceColor','k','markerEdgeColor',[.4 .4 .4]);
h1 = scatter3(catTMtpl.x, catTMtpl.y, catTMtpl.alt,  80, 'filled','markerFaceColor','r','markerEdgeColor',[.4 .4 .4]);
h2 = scatter3(catTM.x,    catTM.y,    catTM.alt,     20, 'filled','markerFaceColor','b','markerEdgeColor',[.4 .4 .4]);
h3 = scatter3(catGES.x,   catGES.y,   catGES.alt,    20, 'filled','markerFaceColor','y','markerEdgeColor',[.4 .4 .4]);

legend([h0,h1,h2,h3],'SED','TK''s TM templates','TK''s TM detections','GES') 

set(gca,'view',[0 90])
set(gca,'view',[0 0])
set(gca,'view',[-50 15])
set(gca,'view',[270 0])









%% Plot m-vs-t for catGES & catTM
[hf,figName]     = plot_2cats_stem_plot(catTM,catGES);
hf.PaperPosition = [0 0 15 4];
print('-depsc2',figName)
print('-dpng'  ,figName)

[hf,figName] = plot_2cats_stem_plot(catTM,catTMtpl);
hf.PaperPosition = [0 0 15 4];
print('-depsc2',figName)
print('-dpng'  ,figName)





%% Find overlap between catalogs
cat1 = catGES;
cat2 = catTM;

dmmax   = 3;
dtmax   = seconds(2);
drmax   = 1000;
matches = find_cat_overlap(cat1,cat2,dtmax,dmmax,drmax)

% Go back to 3d seismicity map
hold on;
ii = 1;
i1 = matches.i1{ii};
i2 = matches.i2{ii}(1);
p1 = plot3(cat1.x(i1),cat1.y(i1),cat1.alt(i1),'xr','markerSize',20);
p2 = plot3(cat2.x(i2),cat2.y(i2),cat2.alt(i2),'xr','markerSize',20);
delete([p1, p2])

% First match, agrees with TK
print_cat_entry(cat1,matches.i1{1})
print_cat_entry(cat2,matches.i2{1})

% Second match, agrees with TK
print_cat_entry(cat1,matches.i1{2})
print_cat_entry(cat2,matches.i2{2}(2))

% Third match, agrees with TK
print_cat_entry(cat1,matches.i1{3})
print_cat_entry(cat2,matches.i2{3})

% Fourth match, agrees with TK
print_cat_entry(cat1,matches.i1{4})
print_cat_entry(cat2,matches.i2{4})

% Overwrite double entry
matches.i2{2} = matches.i2{2}(2);

% How many of the larger events do not have a match?
i1 = cell2mat(matches.i1);
i2 = cell2mat(matches.i2);

sort(cat2.m)
sort(cat2.m(i2))
sum(cat2.m>-1.84)
% >> 45 events with m larger than the smallest match.





%% Select main catalog for further analysis
cat = catGES;
%cat = catTM;




%% Reference point
ref.x   = 0;
ref.y   = 0;
ref.z   = 0;
ref.alt = 1480;
ref.m   = 0;
%ref.t   = cat.t0(find(cat.m>-1.4));
ref.t   = cat.t0(find(cat.m==max(cat.m)));
ref.t   = dateshift(ref.t, 'start', 'day');
ref.i   = [];
cat.ref = ref;


% Distance wrt reference event
cat.rhyp_ref = sqrt( (cat.x-ref.x).^2 +(cat.y-ref.y).^2 +(cat.alt-ref.alt).^2 );
cat.repi_ref = sqrt( (cat.x-ref.x).^2 +(cat.y-ref.y).^2 );
cat.dm       = cat.m-min(cat.m(~isnan(cat.m)))+.1;
cat.dt_days  = days(cat.t0 - ref.t);
%cat.azi   = azimuth_xyz(cat.lat(iref),cat.lon(iref),cat.lat,cat.lon);
cat0      = cat;





%% Plot seismicity overview (so far without tunnel & BH)
cat.plt.x.lim = [-300 -50];
cat.plt.y.lim = [-300 -50];
cat.plt.z.lim = [100 300];
cat.plt.m.lim = [1.1*min(cat.m) 0.9*max(cat.m)];
cat.plt.s.val = 40*cat.dm; 

isOut = cat.x<cat.plt.x.lim(1) | cat.x>cat.plt.x.lim(2) | ...
        cat.y<cat.plt.y.lim(1) | cat.y>cat.plt.y.lim(2) | ...
        cat.z<cat.plt.z.lim(1) | cat.z>cat.plt.z.lim(2);
fprintf(1,sprintf('%i events outside plotting box. \n',sum(isOut)));

cat.plt.c.val = cat.z;
cat.plt.c.lim = cat.plt.z.lim;
cat.plt.c.lim = [150 250];
cat.plt.c.lab = 'Hypocentral depth [m]';
cat.plt.c.n   = 20;
%fg2 = plot_seiscat_overview(cat,0); %2d
fg3 = plot_seiscat_overview(cat,1); %3d

%print('-dpng'  ,[cat.str.figdir,'/seiscat_overview_zhyp'],'-painters')
%print('-depsc2',[cat.str.figdir,'/seiscat_overview_zhyp'],'-painters')
% fg3 = plot_seiscat_overview_other_params(cat,1); %3d










% Magnitudes



































%% OLDER STUFF
plot_seismicity_map_simple
% print('-dpng',sprintf('%s/stem',fbasename))
% print('-dpng',sprintf('%s/stem_woFS',fbasename))

%clf; plot(cat.t0,cat.magnitude,'.k'); title(cat.params.catName)




%% Select sub-catalogue
% cat = cat0;
if false
    isInBox = cat.lon>=-116.55 & cat.lon<-116.45 & cat.lat>=33.47 & cat.lat<33.57;  % Only around this sequence
    cat      = select_subcat(cat,find(isInBox));
    cat      = select_subcat(cat,find(isInBox & istimely));
    
    cat      = select_subcat(cat,find(cat.magnitude>4.7));
    cat      = select_subcat(cat,find(cat.magnitude<2 & isInBox));
    
    tstart    = cat.t0(cat.refeq.idx)-days(1);
    tstart    = cat.t0(cat.refeq.idx)-years(3);
    tend      = cat.t0(end); %tend      = datetime(2019,12,31,'InputFormat','yy,mm,dd');
    istimely  = datenum(cat.t0)>=datenum(tstart) & datenum(cat.t0)<datenum(tend);
    cat      = select_subcat(cat,find(istimely));
end



%% Preparations

% Plot parameters
% ===============
cgrey  = [.2 .2 .2];    % Noble grey
chist  = [.8 .8 .8];    % Shade of grey for historic seismicity
ftSize = 12;

ms.t0        = cat.t0(ims);
ms.magnitude = cat.magnitude(ims);
ms.lat       = cat.lat(ims);
ms.lon       = cat.lon(ims);
ms.depth     = cat.depth(ims);

fs.t0        = cat.t0(i6);
fs.magnitude = cat.magnitude(i6);
fs.lat       = cat.lat(i6);
fs.lon       = cat.lon(i6);
fs.depth     = cat.depth(i6);


% From 190711
ms.ehreloc.lat   = 35.7738;
ms.ehreloc.lon   = -117.5958;
ms.ehreloc.depth = 1.74;
fs.ehreloc.lat   = 35.7043;
fs.ehreloc.lon   = -117.4990;
fs.ehreloc.depth = 9.8900;




% Compute distance to n'th event
% ==============================
ieq  = 10;
ilat = cat.lat(ieq);
ilon = cat.lon(ieq);
iz   = cat.depth(ieq);
hypDist = hypoDistance(ilat,ilon,iz*1e3,cat.lat,cat.lon,cat.depth*1e3);   % Input in [m] & [deg], output in [m]
cat.rh = hypDist/1e3;                                                       % Save as km
epiDist = hypoDistance(ilat,ilon,0     ,cat.lat,cat.lon,0);              
cat.re = epiDist/1e3;                                                     


% Compute angle & distance wrt point
% ==================================
% cat.lon0 = -117.5; % Point where the two main faults meet
% cat.lat0 = 35.68;
cat.lon0 = fs.lon; % Point where the two main faults meet
cat.lat0 = fs.lat;
cat.z0   = 5; 

hypDist  = hypoDistance(cat.lat0,cat.lon0,cat.z0*1e3,cat.lat,cat.lon,cat.depth*1e3);   % Input in [m] & [deg], output in [m]
%epiDist  = hypoDistance(ilat,ilon,0     ,cat.lat,cat.lon,0);              
cat.rh0 = hypDist/1e3;                                                       % Save as km
%cat.re0 = epiDist/1e3;                                                     
cat.azi = azimuth(cat.lat0,cat.lon0,cat.lat,cat.lon);









%% SEISMICITY MAPS
% ================




% Plot seismicity map 3 (cluster evolution)
% =========================================
xlm = [-117.9 -117.3]; ylm = [35.5 36];                                   % Full area but smaller
make_time_vectors

% tstart    = datetime(2019,07,04,16,7,0,'InputFormat','yy,mm,dd,HH,mm,ss');
% tend      = datetime(2019,07,06,16,0,0,'InputFormat','yy,mm,dd,HH,mm,ss');
% dt        = minutes(10);
% %dt        = hours(1)/5;
% tt        = tstart:dt:tend;
% tt        = tstart + minutes([86:10:156, 186:60:1500]');
% tt        = tstart + minutes([86:2:98, 100:5:150, 180:60:1500]');
plot_seismicity_evo_map(cat,tt,xlm,ylm)




% Plot seismicity map 3 (in increments of depth, time, ...)
% =========================================================
zedges = 0:1:15; zedges(end)=100;
nbins  = numel(zedges)-1; 
for ibin = 1:nbins
    
    useme = cat.depth>=zedges(ibin) & cat.depth<zedges(ibin+1);
    catX = select_subcat(cat,useme);
    
    clr = days(catX.t0-ms.t0       ); clm = [-1 4];   clab = 'Days since MS';               nc = 10;
    %sz  = 5*catX.magnitude.^2.5;
    sz  = .3*(catX.dm).^3.5; 

    plot_seismicity_map_scatter
    plot_1907SearlesValley_surface_rupture;
    tstringBin = sprintf('z=%1.0f-%1.0fkm',zedges(ibin),zedges(ibin+1));
    title([tstringBin,': ',tstring])
    
    set(gcf,'PaperPositionMode','auto')
    figname = sprintf('%s/eq_map_z_%02d_%02dkm.png','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new',zedges(ibin),zedges(ibin+1));
    %print('-dpng',figname)
    pause
end


% Plot rotated seismicity
% =======================
catX = cat;

% Rotation clockwise around [x0,y0] by alpha
y0         = 35.87;
x0         = -117.7;
alpha      = -43*pi/180;
mminProf   = -5;        % Minimum magnitude for depth profile
stressdrop = 3e6;       % Stress drop in Pa, for plotting source radii   

% Rotate to see short leg
x0    = -117.64;
y0    = 35.56;
alpha = 48*pi/180;

% Area around main fault
xlmr = [-.2 .6];
ylmr = [-.3 .3];
plot_rotated_seismicity_profile

set(gcf,'PaperPositionMode','auto')
%print('-dpng'  ,'~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot.png')
%print('-depsc2','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot.eps')
%print('-depsc2','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot_ZRTM.png')
%print('-dpng'  ,'~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot_ZRTM.png')
%print('-dpng'  ,'~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot_ZRTM_shortLeg.png')




% Near Garlock Fault
% ==================
plot_rotated_historic_seismicity_Garlock

set(gcf,'PaperPositionMode','auto')
%print('-dpng'  ,'~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot_Garlock.png')
%print('-depsc2','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/maps/new/seismicity_rot_Garlock.eps')





% Plot profile sequence at different dy (fault-normal distances)
% ==============================================================
%clf; histogram(cat.yr)
dyedges = -.1:.02:.1;
ndy     = numel(dyedges)-1;

plot_rotated_seismicity_profile_slices



% Plot profile sequence as function of time
% =========================================
tedges = fs.t0-minutes(2):hours(1):ms.t0+minutes(60);                   % Around FS
%tedges = min(cat.t0):minutes(30):max(cat.t0);
%tedges = fs.t0-minutes(40):minutes(1):fs.t0+minutes(60);                   % Around FS
tedges = fs.t0-minutes(40):minutes(1)/3:fs.t0+minutes(60);                   % Around FS
%tedges = ms.t0-minutes(40):minutes(1):ms.t0+minutes(60);                   % Around MS
nt     = numel(tedges )-1;

% Max normal distance from fault for depth profile
dymax  = 0.03;

% Rotate to see short leg
x0    = -117.64;
y0    = 35.56;
alpha = 48*pi/180;

xall=[]; yall=[]; xnall=[]; ynall=[]; znall=[]; snall=[];

hf = figure(1017); clf; 

for it = 1:nt
    
    print_iteration_numbers(it,nt,'tens')
    tlo   = tedges(it);
    tup   = tedges(it+1);
    useme = cat.t0>=tlo & cat.t0<tup;
    catX = select_subcat(cat,useme);

    plot_rotated_seismicity_profile_anim
    
    %drawnow
    F(it) = getframe(hf);
end


vidName = 'SearlesValley_rot48_nt100';
video = VideoWriter(vidName);
video.FrameRate=1;
open(video)
writeVideo(video,F)
close(video)





%% TIME SEQUENCES 
% ===============

% Stem plot of <y> vs t
% =====================
mmin  = -5;
useme = cat.magnitude>mmin;
%catX = select_subcat(cat,useme);
catX = cat;

% Not accurate, needs correction
%catX.yrkm = catX.yr*deg2km(1);
%catX.dm   = catX.magnitude-min(catX.magnitude)+.1;

x   = catX.t0;        xlab = 'Date';
y   = catX.magnitude; ylab = 'Catalog Magnitude';
%y   = catX.depth;     ylab = 'Catalog Depth [km]';                          
%clr = catX.azi;       clab = 'Azimuth wrt/ Mw6.4 epicenter [°]';  clm = [0 360]; nc = 12;
%clr = catX.yrkm;      clab = 'approximate fault-normal distance [km]';  clm = [-5 5]; nc = 20;
%clr = catX.re0;       clab = 'Distace wrt/ hinge point [km]'; clm = [0 20]; nc = 10;
clr = catX.depth;     clab = 'Catalog Depth [km]';           clm = [5 15]; nc = diff(clm)
%clr = catX.depth;     clab = 'Catalog Depth [km]'; clm = [5 15];
%clr = datenum(catX.t0); clab = 'Date';               clm = [5 7];
%clr = days(cat.t0-cat.t0(1)); clm = [130 150]; clab = 'Days since 19/01/15';
%clr = catX.lat;     clab = 'Latitude [°]'; clm = [35.5 35.8];
%clr = catX.lon;     clab = 'Longitude [°]'; clm = [-117.6 -117.5];
%clr = catX.rh;      clab = 'Hypocentral distance from MS [km]'; clm = [0 15];
%clr = catX.re;      clab = 'Epicentral distance from MS [km]'; clm = [0 15];
s   = 5*(catX.magnitude+2).^2.5; 
s   = .04*catX.dm.^4.5; 
%s   = 50*catX.depth; 

ylm = [-1 30];
ylm = [-1 8];
plot_seismicity_stem_plot




%% SEISMICITY AND MOMENT RATES
% ============================

% Plot seismicity and cumulative moment rates
% ===========================================
tstart    = datetime(2019,07,02,'InputFormat','yy,mm,dd');
tend      = datetime(2019,09,04,'InputFormat','yy,mm,dd');
dt        = days(1)/4;
dt        = hours(1);
%dt        = calweeks(1);
%dt        = calmonths(1);
T         = tstart:dt:tend;

plot_seismicity_rates(cat,T)

set(gcf,'PaperPositionMode','auto')
%print('-dpng',sprintf('%s_srates',fbasename))



% Plot cumulative moment & event numbers
% ======================================
plot_cum_Mo_and_neq



% Plot cluster evolution with increasing time windows
% ===================================================
printme = 1;
% convert -delay 2 -loop 0 GR*.png SearlesValley_Mw6p4_clusterEvo.gif

tstart    = datetime(2019,07,04,16,7,0,'InputFormat','yy,mm,dd,HH,mm,ss');
tend      = datetime(2019,07,06,16,0,0,'InputFormat','yy,mm,dd,HH,mm,ss');
%dt        = days(1);
dt        = hours(1)/5;
%dt        = calweeks(1);
%dt        = calmonths(1);
T         = tstart:dt:tend;
% T         = tstart + minutes([86:5:151, 160:60:1500]');
% T         = tstart + minutes([86,87:94,95:5:115, 120:60:1500]');           % 1min increments after MS
T         = tstart + minutes([86:10:156, 186:60:1500]');
T         = tstart + minutes([86:2:98, 100:5:150, 180:60:1500]');

plot_cluster_evolution



%% Pairplot of everything with everything
% =======================================
%  Add fields: angle wrt/ angle point; distance to MS, ...
dat = [cat.magnitude, cat.depth, datenum(cat.t0), cat.re, cat.lat, cat.lon];
lab = {'M','z [km]', 't_0', 'r_e', 'lat', 'lon'}';
grp = repmat({''},n,1);
clr = [.3 .3 .3];
hp = figure(239); clf; hold on;
pairplot_mam(dat, lab, grp, clr, 'histogram')

annotation('textbox', [.15 .88 1 0.1], ... 
    'String', tstring, ... 
    'EdgeColor', 'none', ... 
    'backgroundColor','none', ... 
    'HorizontalAlignment', 'left', ... 
    'FitBoxToText','on', ... 
    'fontSize',12)

% clear all
% addpath(genpath('~/programs/mapSeis/'))

% Also see scripts study_juice_catalog.m and study_scedc_catalog.m in 
% projects/stf/deconvSTF






%% FMDs
% =====

fmd = calc_FMD_mam(cat.magnitude,0.1);

figure(1002); clf; hold on; box on; grid on;
li = plot(fmd.mCenters,fmd.inc,'ok');
lc = plot(fmd.mCenters,fmd.cum,'r');
set(gca,'yscale','log','ylim',[.9 1.2*max(fmd.cum)])
legend([lc,li],'cumulative FMD','incremental FMD')
ylabel('No. of cases')
xlabel('Magnitude')
title({tstring;'ANNS cat'})

%plot_FMD_in_time_windows
%plot_FMD_in_time_windows_variable % outdated
%plot_FMD_in_depth_windows         % outdated

cat.t0(cat.magnitude>6)

set(gcf,'PaperPositionMode','auto')
%print('-dpng',sprintf('%s/fmd',fbasename))




% Plot FMD in neq windows
% =======================
neq  = numel(cat.magnitude);
dneq = 10;
is  = 1;
ie  = 0;
tstring = sprintf('%s --- neq 1 : %i : %i',cat.params.catName,dneq,neq);

plot_fmd_in_neq_increments

set(gca,'xlim',[0 7])
set(gcf,'PaperPositionMode','auto')
figname = '~/programs/seismo/fig/i39/seismicity/fmd/new/FMD_neqwindows';
%figname = '~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/fmd/new/FMD_neqwindows';
%print('-dpng',figname); %print('-depsc2',figname)



% Plot FMD in time windows
% ========================
tstart = min(cat.t0);
tend   = max(cat.t0);
dt     = caldays(1);
dt     = hours(2);
%t = tstart:caldays(1):tend;
t  = tstart:dt:tend;
nt = numel(t);
tstring = sprintf('QTM --- %s : %s : %s',datestr(tstart),dt,datestr(tstart));

plot_fmd_in_time_increments

set(gca,'xlim',[-1 7])
set(gcf,'PaperPositionMode','auto')
%figname = '~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/fmd/new/FMD_timewindows';
figname = '~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/qtm/new/FMD_timewindows';
%print('-dpng',figname) %print('-depsc2',figname)



% Plot FMDs between three main Eqs
% ================================
% cat = select_subcat(cat,useme);
% cat.params.catName = [cat.params.catName,'; starting at M6.4'];

% FMDs up to main events
fmd1 = calc_FMD_mam(cat.magnitude(1:i6),0.1);
fmd2 = calc_FMD_mam(cat.magnitude(1:i5),0.1);
fmd3 = calc_FMD_mam(cat.magnitude(1:i7),0.1);
fmd4 = calc_FMD_mam(cat.magnitude(1:end),0.1);
lgdStrings = {'19/05/01 - M6.4','19/05/01 - M5.4','19/05/01 - M7.1','19/05/01 - 19/11/14'};

% FMDs in between main events
fmd1 = calc_FMD_mam(cat.magnitude( 1:i6),0.1);
fmd2 = calc_FMD_mam(cat.magnitude(i6:i5),0.1);
fmd3 = calc_FMD_mam(cat.magnitude(i5:i7),0.1);
fmd4 = calc_FMD_mam(cat.magnitude(i7:end ),0.1);
lgdStrings = {'19/05/01 - M6.4','M6.4 - M5.4','M5.4 - M7.1','19/05/01 - 19/11/14'};

figure(1002); clf; hold on; box on; grid on;
set(gca,'yscale','log','xlim',[.1 7.2],'ylim',[.9 1.2*max(fmd4.cum)])
ylabel('N( m\geqm'' )')
xlabel('Magnitude')
title(cat.params.catName)

lc1 = plot(fmd1.mCenters,fmd1.cum,'r','lineWidth',2); legend([lc1]            ,lgdStrings(1));   print('-dpng','FMD_ridgecrest1907_1')
lc2 = plot(fmd2.mCenters,fmd2.cum,'k','lineWidth',2); legend([lc1;lc2]        ,lgdStrings(1:2)); print('-dpng','FMD_ridgecrest1907_2')
lc3 = plot(fmd3.mCenters,fmd3.cum,'b','lineWidth',2); legend([lc1;lc2;lc3]    ,lgdStrings(1:3)); print('-dpng','FMD_ridgecrest1907_3')
lc4 = plot(fmd4.mCenters,fmd4.cum,'m','lineWidth',2); legend([lc1;lc2;lc3;lc4],lgdStrings(1:4)); print('-dpng','FMD_ridgecrest1907_4')

legend([lc1,lc2,lc3,lc4],lgdStrings)

%print('-dpng','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/fmd/scedc/new/FMD_between_mainEqs_ETAS_woIncompleteness')
%print('-dpng','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/fmd/scedc/new/FMD_upto_mainEqs')

cat.t0(i7) - cat.t0(i5)
cat.t0(i5) - cat.t0(i6)



% Plot FMD in depth windows
% =========================
% Sort catalogs by depth
[val,idx] = sort(cat.depth,'ascend');
fieldList = fieldnames(cat);
nfd       = numel(fieldList);
for ifd = 1:nfd
    cat.(fieldList{ifd}) = cat.(fieldList{ifd})(idx); 
end


% Plot FMDs in depth increments with fixed numbers of events
nev = 2500; % Number of earthquakes per depth increment

clf; hold on; box on; grid on; 
set(gca,'yscale','log')
[mFMDC, mFMD] = calc_FMD_mam(cat.magnitude);
plot(mFMDC(1,:),mFMDC(2,:),'lineWidth',2);

neq        = numel(cat.magnitude);
indexList  = (0:nev:neq)+1;
nbins      = numel(indexList)-1;
zmedian    = zeros(nbins,1);
lc         = zeros(nbins,1);
lgdStrings = cell (nbins,1);
for ibin = 1:nbins
    idx           = indexList(ibin):indexList(ibin+1);
    m             = cat.magnitude(idx);
    zmedian(ibin) = median(cat.depth(idx));
    [mFMDC, mFMD] = calc_FMD_mam(m);
    lc(ibin)         = plot(mFMDC(1,:),mFMDC(2,:));
    lgdStrings{ibin} = sprintf('z ~ %3.1fkm',zmedian(ibin));
end

hl = legend(lc,lgdStrings);
xlabel('Magnitude')
ylabel('N(m>=m'')')

set(gcf,'PaperPositionMode','auto')
%print('-dpng','FMD_depthwindows.png')




% FMD with shifted magnitudes at m<m'
% ===================================
figure(22); clf; hold on; box on; grid on; 
set(gca,'yscale','log')
xlabel('Magnitude')
ylabel('N(m>=m'')')
title(cat.params.catName)

fmd0 = calc_FMD_mam(cat.magnitude);
lc   = plot(fmd0.mCenters,fmd0.cum,'r','lineWidth',2);        

mPrime = 3.5;

m2            = cat.magnitude;
m2(m2<mPrime) = m2(m2<mPrime)-.3;
fmd           = calc_FMD_mam(m2);
lc            = plot(fmd.mCenters,fmd.cum,'k','lineWidth',1);        

m3            = cat.magnitude;
m3(m3<mPrime) = m3(m3<mPrime)-.6;
fmd           = calc_FMD_mam(m3);
lc            = plot(fmd.mCenters,fmd.cum,'k','lineWidth',1);        
%print('-dpng','~/programs/seismo/fig/i39/specialEQs/ridgecrest1907/seismicity/fmd/new/fmd_const_mShift')






%% APPENDIX
% =========
% % Plot seismicity rates
% ts  = tstart-days(1);
% neq = []; 
% T   = [];
% while ts<tend
%     ts  = ts + days(1);
%     te  = ts + days(1);
%     n   = sum(cat.t0>ts & cat.t0<=te);
%     neq = [neq; n]; 
%     T   = [T  ; ts];
% end
% % ceil(ts)
% % datetime('','format','yyyy')
% 
% clf; hold on; grid on; box on;
% plot(T,neq)
% ylabel('No. of events per time increment')
% set(gcf,'PaperPositionMode','auto')
% %print('-dpng','seismicity_rates.png')
