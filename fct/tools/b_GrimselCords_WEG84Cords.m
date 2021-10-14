clc
close all
clear all


%% General initialization
if isunix
    idcs   = strfind(pwd,'/');
    curdir = pwd;
    datadir = [curdir(1:idcs(end)) '01_Data/'];
    resultdir = [curdir(1:idcs(end)) '03_Result/'];   
    if exist(datadir, 'dir')
        addpath('../01_Data/');
    else
        mkdir(datadir);
    end
    if exist(resultdir, 'dir')
        addpath('../01_Data/');
    else
        mkdir(resultdir);
    end
elseif ispc
    idcs   = strfind(pwd,'\');
    curdir = pwd;
    datadir = [curdir(1:idcs(end)) '01_Data\'];
    resultdir = [curdir(1:idcs(end)) '03_Result\'];   
    if exist(datadir, 'dir')
        addpath('..\01_Data\');
    else
        mkdir(datadir);
    end
    if exist(resultdir, 'dir')
        addpath('..\01_Data\');
    else
        mkdir(resultdir);
    end
end

pic_width = 1600/1.5;
pic_hight = 900/2.5;



%% Code
%%Load Grimsel Cords
load('fct/tools/cords_hyd_stimualtion_v02.mat')

%%Set reference new, cause you wan't all depths to be positive

z_n = 50;
% x cord, y cord center in borehole array
xy_n = mean(HS_cords.GMuG_cords_norm(16:23, 2:3));

%%Reference Grimsel Origin in LV03 cords
ref_grimsel_LV03 = [667400, 158800, 1700]; % East, North, Height [m, m, m]
new_grimsel_ref  = ref_grimsel_LV03 + [xy_n, z_n];
new_grimsel_ref  = [667449.35949, 158904.02465, 1750.00000];

%%Translation in WEG 84 Coordinates (https://www.swisstopo.admin.ch/de/karten-daten-online/calculation-services/navref.html)
ref_grimsel_WEG84     = [8.317911007, 46.57707751,  1752.287]; % Long, Lat, Height [°, °, m]
ref_grimsel_WEG84_new = [8.318570126, 46.578008239, 1802.292];


spheroid = wgs84Ellipsoid; % Reference ellipsoid for World Geodetic System 1984 [m]
fid_out = fopen('Grimsel_receiver_locations.csv','w'); %open file and create fid
% fprintf(fid_out,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','Receiver', 'Longitude[deg]', 'Latitude[deg]', 'Elevation[m]', 'Azimuth[m]', 'Dip[deg]', 'OnTime[yyyy]', '[mm]', '[dd]', 'OffTime[yyyy]','[mm]', '[dd]');

% Cords of interest
cords_HS = HS_cords.GMuG_cords_norm(1:32,:);
cords_HS(:,2:4) = cords_HS(:,2:4) - [xy_n(1), xy_n(2), z_n];

for i = 1:length(cords_HS)
    [lat,lon,h] = ned2geodetic(cords_HS(i,3) * 1000,cords_HS(i,2) * 1000, -cords_HS(i,4) * 1000, ref_grimsel_WEG84_new(2), ref_grimsel_WEG84_new(1), ref_grimsel_WEG84_new(3), spheroid);
    fprintf(fid_out,'%s\t%0.4f\t%0.4f\t%0.3f\n', [ 'R' sprintf('%03d', cords_HS(i,1))], lon, lat, h/1000);     % I want depth as positive, therefore I need a "-" again
end

fclose(fid_out); %close file


