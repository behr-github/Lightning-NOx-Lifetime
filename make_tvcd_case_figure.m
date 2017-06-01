function [ Base, Cases ] = make_tvcd_case_figure(  )
%MAKE_TVCD_CASE_FIGURE - Produce the supplemental figure looking at how each kinetics change altered the model tVCDs

mydir = fileparts(mfilename('fullpath'));
aux_file = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization/DailyOMI/All2012Daily/JPLnoMPN-Pickp0-OMI_aux.mat';

server_dir = '/Volumes/share2/USERS/LaughnerJ/MPN_Project/Workspaces/Pickering_parameterization/OMI_obs';
Base.file = 'JPLnoMPN-Pickp0-OMI_NO2.mat';

Cases(1).file = 'JPLwMPN-Pickp0-OMI-TIME-SER_NO2.mat';
Cases(1).title = '(a) % Change in tVCD_{NO2}, MPN Case';
Cases(2).file = 'JPLnoMPN_N2O5-Pickp0-OMI-TIME-SER_NO2.mat';
Cases(2).title = '(b) % Change in tVCD_{NO2}, N_2O_5 Case';
Cases(3).file = 'JPLnoMPN_PNA-Pickp0-OMI-TIME-SER_NO2.mat';
Cases(3).title = '(c) % Change in tVCD_{NO2}, PNA Case';
Cases(4).file = 'HendnoMPN-Pickp0-OMI-TIME-SER_NO2.mat';
Cases(4).title = '(d) % Change in tVCD_{NO2}, HNO_3 Case';

cb_limits = [-15 15];

% Load the auxiliary file and cut down the times contained in it
aux_struct = load_and_extract_var(aux_file); % loads JPLnoMPN_Pickp0_OMI_aux
aux_struct = cut_times(aux_struct);

% Do the same for the base file, also calculate the columns
Base.gc = load_and_extract_var(fullfile(server_dir, Base.file));
Base.gc = cut_times(Base.gc);
Base.gc = calc_vcd(Base.gc, aux_struct);

% and for each case
[gloncorn, glatcorn]=geos_chem_corners;
[xx,yy] = misc_gc_plotting_fxns('region_xxyy','na');
figure;
for a=1:numel(Cases)
    Cases(a).gc = load_and_extract_var(fullfile(server_dir, Cases(a).file));
    Cases(a).gc = cut_times(Cases(a).gc);
    Cases(a).gc = calc_vcd(Cases(a).gc, aux_struct);
    
    % Also go ahead and compute the relative difference to Base
    Cases(a).rdel = reldiff(Cases(a).gc.Columns, Base.gc.Columns)*100;
    
    % And make the plots - no need for a separate for loop
    subplot(2,2,a);
    pcolor(gloncorn(yy,xx), glatcorn(yy,xx), nanmean(Cases(a).rdel(xx,yy,:),3)');
    colorbar;
    colormap jet;
    caxis(cb_limits);
    state_outlines('k');
    title(Cases(a).title);
end

savefig(fullfile(mydir, 'GEOS-Chem-Cases-Figure.fig'));
saveas(gcf,fullfile(mydir,'GEOS-Chem-Cases-Figure.png'));

end

function S = calc_vcd(S, aux_struct)
try
    S = rmfield(S,'Columns');
catch err
    if strcmp(err.identifier,'MATLAB:rmfield:InvalidFieldname')
        fprintf('No "Columns" field to remove');
    else
        rethrow(err);
    end
end

S_tmp = cat(2, S, aux_struct);
S_tmp = integrate_geoschem_profile(S_tmp, 1e-9);
S = S_tmp(1);

end

function var_out = load_and_extract_var(filename)
fprintf('Loading %s\n', filename);
D = load(filename);
fns = fieldnames(D);
if numel(fns) > 1
    error('loading_var:too_may_vars','More than one variable in file %s', filename);
end
var_out = D.(fns{1});
end

function S = cut_times(S)
E = JLLErrors;
dd = S(1).tVec >= datenum('2012-05-01') & S(1).tVec < datenum('2012-06-30'); % match the DC3 times, which end on June 29
dd2 = dd;
dd2(find(dd2,1,'last')+1) = true;
for a=1:numel(S)
    if ndims(S(a).dataBlock) == 4
        S(a).dataBlock = S(a).dataBlock(:,:,:,dd);
    elseif ndims(S(a).dataBlock) == 3
        % tropopause level is only 3D, not 4D
        S(a).dataBlock = S(a).dataBlock(:,:,dd);
    else
        E.notimplemented('Do not know how to cut down dataBlock with ndims == %d', ndims(S(a).dataBlock));
    end
    S(a).tVec = S(a).tVec(dd);
    S(a).tEdge = S(a).tEdge(dd2);
end
end
