%Matches results from model results with observation

pns_obs(1:23,1:21,1:90)=NaN;
pns_ct(1:23,1:21,1:90)=NaN;
co_obs(1:23,1:21,1:90)=NaN;
co_ct(1:23,1:21,1:90)=NaN;
hno3_obs(1:23,1:21,1:90)=NaN;
nox_obs(1:23,1:21,1:90)=NaN;
o3_obs(1:23,1:21,1:90)=NaN;
press_obs(1:23,1:21,1:90)=NaN;
nox_ct(1:23,1:21,1:90)=0;
hno3_ct(1:23,1:21,1:90)=0;
o3_ct(1:23,1:21,1:90)=0;
press_ct(1:23,1:21,1:90)=0;
no2_obs(1:23,1:21,1:90)=NaN;
no2_ct(1:23,1:21,1:90)=0;
jday_obs(1:23,1:21,1:90)=NaN;
jday_ct(1:23,1:21,1:90)=0;
h2o_obs(1:23,1:21,1:90)=NaN;
h2o_ct(1:23,1:21,1:90)=0;
ans_obs(1:23,1:21,1:90)=NaN;
ans_ct(1:23,1:21,1:90)=0;

tday = dc3_merges_ut(1,2);
press_avg = 0;hno3_avg = 0;nox_avg = 0;o3_avg = 0;
for d = 1:90
    
    while d+121 == tday
        tday = (dc3_merges_ut(ct,2));
        tlondc3 = (dc3_merges_ut(ct,7)-360);
        tlatdc3 = (dc3_merges_ut(ct,6));
                
        for cnt = 2:23;
            if tlondc3 <= NA_lon(cnt) && tlondc3 >= NA_lon(cnt-1);
                obs_lon_num = cnt;obs_lon_num;
            end
        end
        
        for c = 2:21;
            if tlatdc3 <= NA_lat(c) && tlatdc3 >= NA_lat(c-1);
                obs_lat_num = c;obs_lat_num;
            end
        end
        
        %press_obs(obs_lon_num-1,obs_lat_num-1,d)=dc3_merges_ut
        
        if isnan(h2o_obs(obs_lon_num-1,obs_lat_num-1,d))==1
            h2o_obs(obs_lon_num-1,obs_lat_num-1,d)=dc3_merges_ut(ct,66);
            h2o_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            h2o_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(h2o_obs(obs_lon_num-1,obs_lat_num-1,d)) + ...
                dc3_merges_ut(ct,66));
            h2o_ct(obs_lon_num-1,obs_lat_num-1,d)=h2o_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        if isnan(ans_obs(obs_lon_num-1,obs_lat_num-1,d))==1
            ans_obs(obs_lon_num-1,obs_lat_num-1,d)=dc3_merges_ut(ct,76);
            ans_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            ans_obs(obs_lon_num-1,obs_lat_num-1,d) = squeeze(ans_obs(obs_lon_num-1,obs_lat_num-1,d)) + ...
                dc3_merges_ut(ct,76);
            ans_ct(obs_lon_num-1,obs_lat_num-1,d)=ans_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        if isnan(press_obs(obs_lon_num-1,obs_lat_num-1,d))==1
            press_obs(obs_lon_num-1,obs_lat_num-1,d) = dc3_merges_ut(ct,9);
            press_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            press_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(press_obs(obs_lon_num-1,obs_lat_num-1,d)) + ...
                dc3_merges_ut(ct,9));
            press_ct(obs_lon_num-1,obs_lat_num-1,d)=press_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        if isnan(jday_obs(obs_lon_num-1,obs_lat_num-1,d))==1
            jday_obs(obs_lon_num-1,obs_lat_num-1,d) = dc3_merges_ut(ct,2);
            jday_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            jday_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(jday_obs(obs_lon_num-1,obs_lat_num-1,d)) + ...
                dc3_merges_ut(ct,2));
            jday_ct(obs_lon_num-1,obs_lat_num-1,d)=jday_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        
        if isnan(hno3_obs(obs_lon_num-1,obs_lat_num-1,d)) == 1
            hno3_obs(obs_lon_num-1,obs_lat_num-1,d) = hno3_dc3(ct);
            if isnan(hno3_dc3(ct))==0
                hno3_ct(obs_lon_num-1,obs_lat_num-1,d)=hno3_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
            end
        else
            hno3_obs(obs_lon_num-1,obs_lat_num-1,d) = (nansum([squeeze(hno3_obs(obs_lon_num-1,obs_lat_num-1,d))+...
                hno3_dc3(ct)]));
            if isnan(hno3_dc3(ct))==0
                hno3_ct(obs_lon_num-1,obs_lat_num-1,d)=hno3_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
            end
        end
        
        if isnan(no2_obs(obs_lon_num-1,obs_lat_num-1,d)) == 1
            no2_obs(obs_lon_num-1,obs_lat_num-1,d)=dc3_merges_ut(ct,73);
            no2_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            no2_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(no2_obs(obs_lon_num-1,obs_lat_num-1,d))+...
                dc3_merges_ut(ct,73));
            no2_ct(obs_lon_num-1,obs_lat_num-1,d)=no2_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        if isnan(co_obs(obs_lon_num-1,obs_lat_num-1,d)) == 1
            co_obs(obs_lon_num-1,obs_lat_num-1,d)=dc3_merges_ut(ct,63);
            co_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            co_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(co_obs(obs_lon_num-1,obs_lat_num-1,d))+...
                dc3_merges_ut(ct,63));
            co_ct(obs_lon_num-1,obs_lat_num-1,d)=co_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        if isnan(pns_obs(obs_lon_num-1,obs_lat_num-1,d)) == 1
            pns_obs(obs_lon_num-1,obs_lat_num-1,d)=dc3_merges_ut(ct,75);
            pns_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            pns_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(pns_obs(obs_lon_num-1,obs_lat_num-1,d))+...
                dc3_merges_ut(ct,75));
            pns_ct(obs_lon_num-1,obs_lat_num-1,d)=pns_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        %hno3_obs(cnt,c,d)
        %ct
        %d
        %pause
        if isnan(o3_obs(obs_lon_num-1,obs_lat_num-1,d)) == 1
            o3_obs(obs_lon_num-1,obs_lat_num-1,d) = dc3_merges_ut(ct,72);
            o3_ct(obs_lon_num-1,obs_lat_num-1,d)=1;
        else
            o3_obs(obs_lon_num-1,obs_lat_num-1,d) = (squeeze(o3_obs(obs_lon_num-1,obs_lat_num-1,d))+...
                dc3_merges_ut(ct,72));
            o3_ct(obs_lon_num-1,obs_lat_num-1,d)=o3_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
        end
        
        if isnan(nox_obs(obs_lon_num-1,obs_lat_num-1,d)) == 1
            nox_obs(obs_lon_num-1,obs_lat_num-1,d) = dc3_merges_ut(ct,71).*1000 + dc3_merges_ut(ct,73);
            if isnan((dc3_merges_ut(ct,71).*1000+dc3_merges(ct,73))) == 0
                nox_ct(obs_lon_num-1,obs_lat_num-1,d)=nox_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
            end
        else
            nox_obs(obs_lon_num-1,obs_lat_num-1,d) = (nansum([squeeze(nox_obs(obs_lon_num-1,obs_lat_num-1,d)),...
                (dc3_merges_ut(ct,71).*1000+dc3_merges(ct,73))]));
            if isnan((dc3_merges_ut(ct,71).*1000+dc3_merges(ct,73))) == 0
                nox_ct(obs_lon_num-1,obs_lat_num-1,d)=nox_ct(obs_lon_num-1,obs_lat_num-1,d)+1;
            end
        end
        
        ct = ct + 1;
        
        %d + 121
        %tday
        %pause
        if ct == length(dc3_merges_ut(:,1))+1
            %ct
            break
        end
    end
            if ct == length(dc3_merges_ut(:,1))+1
            %ct
            break
        end
end


%thno3cts=find(hno3_ct < 3);hno3_ct(thno3cts)=NaN;
tnoxcts = find(nox_ct < 10);nox_ct(tnoxcts)=NaN;
size(tnoxcts)
thno3cts = find(hno3_ct < 10);hno3_ct(thno3cts)=NaN;
to3cts = find(o3_ct < 5);o3_ct(to3cts)=NaN;
th2octs = find(h2o_ct < 10);h2o_ct(th2octs)=NaN;
tno2cts = find(no2_ct < 10);no2_ct(tno2cts)=NaN;
tanscts = find(ans_ct < 10);ans_ct(tanscts)=NaN;
ans_obs = ans_obs./ans_ct;
nox_obs = nox_obs./nox_ct;
hno3_obs = hno3_obs./hno3_ct;
o3_obs = o3_obs./o3_ct;
press_obs = press_obs./press_ct;
no2_obs = no2_obs./no2_ct;
jday_obs = jday_obs./jday_ct;
pns_obs = pns_obs./pns_ct;
h2o_obs = h2o_obs./h2o_ct;