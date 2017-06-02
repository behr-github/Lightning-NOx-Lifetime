%Loads results from model runs
ans_dissert(1:23,1:21,1:47,1:60)=NaN;
switch what_to_compare
    case 'fk_ans'
        load('fkans_o3.mat')
        load('fkans_no.mat')
        load('fkans_no2_sat.mat')
        load('fkans_no2.mat')
        load('fkans_fkans.mat')
        load('fkans_hno3.mat')
        
        hno3 = fkans_hno3;clear fkans_hno3
        nox = fkans_no + fkans_no2;clear fkans_no
        no2 = fkans_no2; clear fkans_no2
        no2_sat = fkans_no2_sat; clear fkans_no2
        o3 = fkans_o3; clear fkans_o3
        ans_dissert = fkans_fkans; clear fkans_fkans
        
    
    case 'Base'
        load('hno3_wompn_base.mat');
        load('no2_wompn_base.mat');
        load('no_wompn_base.mat');
        load('o3_wompn_base.mat');
        load('h2o_dc3.mat');
        load('n2o5_night.mat');
        %load('no2_dc3.mat');
        load('oh_dc3.mat');
        load('no2_wompn_base_sat.mat');
        
        hno3 = hno3_wompn_base;clear hno3_wompn_base
        nox = no2_wompn_base+no_wompn_base;
            no2 = no2_wompn_base;clear no2_wompn_base no_wompn_base
        o3 = o3_wompn_base;clear o3_wompn_base
        h2o = h2o_dc3;clear h2o_dc3
        n2o5 = n2o5_night;clear n2o5_night
        no2_sat = no2_wompn_base_sat;clear no2_dc3;
        oh = oh_dc3;clear oh_dc3;
        
        
    case 'Base_33'
        load('hno3_wompn_sens.mat');
        load('no2_wompn_sens.mat');
        load('no_wompn_sens.mat');
        load('o3_wompn_sens.mat');
        
        hno3 = hno3_wompn_sens;clear hno3_wompn_sens
        nox = no2_wompn_sens + no_wompn_sens;
            no2 = no2_wompn_sens;clear no2_wompn_sens no_wompn_sens
        o3 = o3_wompn_sens;clear o3_wompn_sens
        
    case 'MPN'
        load('hno3_wmpn_jpl.mat');
        load('no2_wmpn_jpl.mat');
        load('no_wmpn_jpl.mat');
        load('o3_wmpn_jpl.mat');
        load('no2_wmpn_jpl_sat.mat');
        
        hno3 = hno3_wmpn_jpl;clear hno3_wmpn_jpl
        nox = no2_wmpn_jpl+no_wmpn_jpl;
            no2 = no2_wmpn_jpl;clear no2_wmpn_jpl no_wmpn_jpl
        o3 = o3_wmpn_jpl;clear o3_wmpn_jpl
        no2_sat = no2_wmpn_jpl_sat;clear no2_wmpn_jpl_sat;

        
    case 'PNA'
        load('hno3_wompn_slow_pna.mat');
        load('no2_wompn_slow_pna.mat');
        load('no_wompn_slow_pna.mat');
        load('o3_wompn_slow_pna.mat');
        load('no2_wompn_slow_pna_sat.mat');
        
        hno3 = hno3_wompn_slow_pna; clear hno3_wompn_slow_pna
        nox = no2_wompn_slow_pna + no_wompn_slow_pna;
            clear no2_wompn_slow_pna no_wompn_slow_pna
        o3 = o3_wompn_slow_pna; clear o3_wompn_slow_pna
        no2_sat = no2_wompn_slow_pna_sat; clear no2_wompn_slow_pna_sat
    
    case 'N2O5'
        load('hno3_wompn_slow_n2o5.mat');
        load('no2_wompn_slow_n2o5.mat');
        load('no_wompn_slow_n2o5.mat');
        load('o3_wompn_slow_n2o5.mat');
        load('no2_wompn_slow_n2o5_sat.mat');
        
        hno3 = hno3_wompn_slow_n2o5; clear hno3_wompn_slow_n2o5
        nox = no2_wompn_slow_n2o5 + no_wompn_slow_n2o5;
            clear no2_wompn_slow_n2o5 no_wompn_slow_n2o5
        o3 = o3_wompn_slow_n2o5; clear o3_wompn_slow_n2o5
        no2_sat = no2_wompn_slow_n2o5_sat;clear no2_wompn_slow_n2o5_sat;
        
    case 'HNO3'
        load('hno3_wompn_hend.mat');
        load('no2_wompn_hend.mat');
        load('no_wompn_hend.mat');
        load('o3_wompn_hend.mat');
        load('no2_wompn_hend_sat.mat');
        
        hno3 = hno3_wompn_hend; clear hno3_wompn_hend
        nox = no2_wompn_hend + no_wompn_hend;
            clear no2_wompn_hend no_wompn_hend
        o3 = o3_wompn_hend; clear o3_wompn_hend
        no2_sat = no2_wompn_hend_sat;clear no2_wompn_hend_sat;
    
    case 'Dissert'
        load('hno3_wmpn_base.mat');
        load('no2_wmpn_base.mat');
        load('no_wmpn_base.mat');
        load('o3_wmpn_base.mat');
        load('no2_wmpn_base_sat.mat');
        
        hno3 = hno3_wmpn_base;clear hno3_wmpn_base
        nox = no2_wmpn_base+no_wmpn_base;
            no2 = no2_wmpn_base;clear no2_wmpn_base no_wmpn_base
        o3 = o3_wmpn_base;clear o3_wmpn_base
        no2_sat = no2_wmpn_base_sat;clear no2_wmpn_base_sat;
        
    case 'Dissert_33'
        load('hno3_wmpn_sens.mat');
        load('no2_wmpn_sens.mat');
        load('no_wmpn_sens.mat');
        load('o3_wmpn_sens.mat');
        load('h2o_dc3.mat');
        load('n2o5_night.mat');
        load('no2_dc3.mat');
        load('oh_dc3.mat');
        load('ans_dissert.mat');
        
        hno3 = hno3_wmpn_sens;clear hno3_wmpn_sens
        nox = no2_wmpn_sens+no_wmpn_sens;
            no2 = no2_wmpn_sens;clear no2_wmpn_sens no_wmpn_sens
        o3 = o3_wmpn_sens;clear o3_wmpn_sens
        h2o = h2o_dc3;clear h2o_dc3;
        n2o5 = n2o5_night;clear n2o5_night
        no2_sat = no2_dc3;clear no2_dc3;
        oh = oh_dc3;clear oh_dc3;
        
end

load('numden.mat');
load('bxheight.mat');
load('press.mat');
load('trop_hght.mat');