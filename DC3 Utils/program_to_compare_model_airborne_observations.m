%clc
%clear all
%close all

%add paths for folders that contain DC3 observations or GEOS-Chem Lat/Lon

addpath E:\DC3_SEAC4RS\GEOS_Chem_Utils%C:\Users\Benjamin\Desktop\DC3_SEAC4RS\MPN_prod_satellite\GEOS_Chem_Utils
%Path for GEOS Lat/Lon
load('geos_lat.mat')
load('geos_lon.mat')
%Definition of NA region
tlon=find(lon >= -120 & lon <= -65);tlat = find(lat >= 20 & lat <= 60);
NA_lon = lon(tlon); NA_lat = lat(tlat);

%addpath C:/Users/Benjamin/Desktop/DC3_SEAC4RS/DC3_SEAC4RS_INTEX_Paper
%Path for DC3 Obs

%load_dc3 %Loads all DC3 observations into 1 matrix

%clear dc3_0518 dc3_0519 dc3_0521 dc3_0525 dc3_0526 dc3_0529 dc3_0530 dc3_0601 ...
 %   dc3_0602 dc3_0605 dc3_0606 dc3_0607 dc3_0611 dc3_0615 dc3_0616 dc3_0616 ...
  %  dc3_0617 dc3_0621 dc3_0622

%Defines observations to compare against for NOx, O3, and HNO3 during DC3

%hno3_dc3 = nanmean([dc3_merges(:,77) dc3_merges(:,261)]')';

%tut_dc3 = find(dc3_merges(:,9) < 350 & dc3_merges(:,9) > 200 & ...(dc3_merges(:,71).*1e3+dc3_merges(:,73))./hno3_dc3 < 5 & ...
 %   ((dc3_merges(:,7)-360 > -90 & dc3_merges(:,1)./3600 - 4 >= 16 & dc3_merges(:,1)./3600 - 4 <= 20) | ...
  %  (dc3_merges(:,7)-360 <= -90 & dc3_merges(:,7)-360 > -105 & dc3_merges(:,1)./3600 - 5 >= 16 & dc3_merges(:,1)./3600 - 5 <= 20) | ...
   % (dc3_merges(:,7)-360 <= -105 & dc3_merges(:,1)./3600 - 6 >= 16 & dc3_merges(:,1)./3600 - 6 <= 20)));... & ...
    %(dc3_merges(:,71).*1e3+dc3_merges(:,73))./(hno3_dc3+dc3_merges(:,71).*1e3+dc3_merges(:,73)+dc3_merges(:,75)) <= .6 & ...
    %(dc3_merges(:,71).*1e3+dc3_merges(:,73)) <= 1000 dc3_merges(:,72)./dc3_merges(:,63) < 1.5 & ...
    %dc3_merges(:,72) < 150 &);
    %(dc3_merges(:,71).*1e3+dc3_merges(:,73))./((hno3_dc3)) > 0 );%& ...
    %(dc3_merges(:,71).*1e3 + dc3_merges(:,73))./(dc3_merges(:,71).*1e3 + dc3_merges(:,73) +...
    %dc3_merges(:,75) + dc3_merges(:,76) + .95.*(dc3_merges(:,77))) < .6);% & ...
    %dc3_merges(:,7)-360 > -105 );
%pause
%dc3_merges_ut = dc3_merges(tut_dc3,:);


what_to_compare = input('Which model do you want to compare against? [Base Base_33 MPN N2O5 PNA HNO3 Dissert Dissert_33 fk_ans] ','s');

load_model_results

remove_strat_bl_results

no2_col_dens_calc

%match_flight_model

%match_results_from_flight

%find_model_results_at_same_press