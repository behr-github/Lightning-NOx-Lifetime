function [  ] = save_mpn_gc_netcdfs( )
%Save the GEOS-Chem .mat files as netCDF files

save_dir = '/Volumes/share2/USERS/LaughnerJ/MPN_Project/Workspaces/netCDFs';
cases_mat_names = {'JPLnoMPN-Pickp0', 'JPLnoMPN_N2O5-Pickp0','JPLnoMPN_PNA-Pickp0','HendnoMPN-Pickp0','HendwMPN-PNA-N2O5-Pickp0','HendwMPN-PNA-N2O5-Pickp33'};
cases_save_names = {'BaseCase', 'N2O5Case','PNACase','HNO3Case','UpdatedCase','UpdatedPlus33Case'};
dc3_variables = {'HNO3','NO2_IJAVG','NO_IJAVG','O3','PSURF','TP'};
dc3_descriptions = {'HNO3 averaged between 1600 and 2000 local time';
    'NO2 averaged between 1600 and 2000 local time';
    'NO averaged between 1600 and 2000 local time';
    'O3 averaged between 1600 and 2000 local time';
    'Bottom grid box pressure averaged between 1600 and 2000 local time';
    'Model level at which the tropopause resides (may be fractional), averaged between 1600 and 2000 local time'};

% Ancillary parameters
base_dir = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization';
server_dir = '/Volumes/share2/USERS/LaughnerJ/MPN_Project/Workspaces/Pickering_parameterization/OMI_obs';

% Production data
fprintf('Loading Area data\n');
Area = load_and_extract(fullfile(base_dir,'Daily24hrs','All2012Daily','Area.mat'));
Area.description = 'GEOS-Chem grid cell areas';
fprintf('Saving Area data\n');
gcstruct2ncdf(fullfile(save_dir, 'GEOS-Chem-Areas.nc'), 'do_tedge', Area);

fprintf('Loading Base case production data\n');
BaseCaseProd = load_and_extract(fullfile(base_dir,'Daily24hrs','All2012Daily','JPLnoMPN-Pickp0-Prod.mat'));
BaseCaseProd = prod_descriptions(BaseCaseProd, 'Base Case');
fprintf('Saving Base case production data\n');
gcstruct2ncdf(fullfile(save_dir, 'BaseCaseNOProduction.nc'), 'do_tedge', BaseCaseProd);

UpdatedCaseProf = load_and_extract(fullfile(base_dir,'Daily24hrs','All2012Daily','HendwMPN-Pickp0-Prod.mat'));
UpdatedCaseProf = prod_descriptions(UpdatedCaseProf, 'Updated Case');
fprintf('Saving Updated case production data\n');
gcstruct2ncdf(fullfile(save_dir, 'UpdatedCaseNOProduction.nc'), 'do_tedge', UpdatedCaseProf);

UpdatedP33CaseProd = load_and_extract(fullfile(base_dir,'Daily24hrs','All2012Daily','HendwMPN-Pickp33-Prod.mat'));
UpdatedP33CaseProd = prod_descriptions(UpdatedP33CaseProd, 'Updated +33% Case');
gcstruct2ncdf(fullfile(save_dir, 'UpdatedPlus33CaseNOProduction.nc'), 'do_tedge', UpdatedP33CaseProd);

fprintf('Loading ancillary data\n')
Data = load_and_extract(fullfile(base_dir, 'DailyOMI', 'All2012Daily', 'JPLnoMPN-Pickp0-OMI_aux.mat'));
Data(1).description = 'Pressure at the bottom of each box';
Data(2).description = 'Height of model grid box in meters';
Data(3).description = 'Temperature of the model grid box';
Data(4).description = 'Number density of air in the model grid box';
Data(5).description = 'Model level at which the tropopause resides (may be fractional)';
Data(5).tVec = floor(Data(5).tVec);
fprintf('Saving ancillary data\n');
gcstruct2ncdf(fullfile(save_dir, 'AncillaryModelParameters.nc'), Data);

%OMI NO2 for all cases
for a=1:numel(cases_mat_names);
   fprintf('Loading %s\n', cases_mat_names{a});
   Data = load_and_extract(fullfile(server_dir,sprintf('%s-OMI-TIME-SER_NO2.mat', cases_mat_names{a})));
   Data.description = 'NO2 concentrations averaged between 1200 and 1400 local time';
   save_name = sprintf('%s-OMI_NO2_1200-1400.nc', cases_save_names{a});
   fprintf('Writing %s\n', save_name);
   gcstruct2ncdf(fullfile(save_dir, save_name), Data);
end

% DC3 data
[xx,yy] = misc_gc_plotting_fxns('region_xxyy','na');
[glon,glat] = geos_chem_centers('2x25');
glon = glon(xx);
glat = glat(yy);
for a=1:numel(cases_mat_names)
    fprintf('Working on %s\n', cases_mat_names{a});
    this_name = strrep(cases_mat_names{a}, '_','-'); % subfolders followed different convention
    this_name = strrep(this_name, 'Pick','');
    data_cell = cell(1,numel(dc3_variables));
    for b=1:numel(dc3_variables)
        file_name = sprintf('%s-DC3_%s.mat', cases_mat_names{a}, dc3_variables{b});
        fprintf('  Loading %s\n', file_name);
        data_cell{b} = load_and_extract(fullfile(base_dir, 'DailyDC3', this_name, file_name));
        data_cell{b}.description = dc3_descriptions{b};
        data_cell{b}.dataBlock = data_cell{b}.dataBlock(xx,yy,:,:); % should be fine whether it's 3d or 4d
    end
    
    Data = cat(1,data_cell{:});
    save_name = sprintf('%s-DC3_1600-2000.nc', cases_save_names{a});
    fprintf('Saving %s\n', save_name);
    gcstruct2ncdf(fullfile(save_dir, save_name), 'lon', glon, 'lat', glat, Data);
end



end

function Data = load_and_extract(fname)
D = load(fname);
fns = fieldnames(D);
if numel(fns) ~= 1
    error('load:multiple_vars','Not exactly one variable')
end
Data = D.(fns{1});
for i=1:numel(Data)
    if iscolumn(Data(i).modelName)
        Data(i).modelName = strtrim(Data(i).modelName');
    end
    Data(i).tVec = floor(Data(i).tVec);
end
end

function S = prod_descriptions(S, case_name)
for a=1:numel(S)
    S(a).description = sprintf('%s %s production', case_name, S(a).fullCat);
end
end
