%Gets model results at the same pressure levels as flight

model_hno3(1:23,1:21,1:90)=NaN;
model_nox(1:23,1:21,1:90)=NaN;
model_o3(1:23,1:21,1:90)=NaN;
model_press(1:23,1:21,1:90)=NaN;
model_no2(1:23,1:21,1:90)=NaN;
model_ans(1:23,1:21,1:90)=NaN;
%model_h2o(1:23,1:21,1:90)=NaN;
%model_n2o5(1:23,1:21,1:90)=NaN;
%model_no2_sat(1:23,1:21,1:90)=NaN;
%model_oh(1:23,1:21,1:90)=NaN;

for cnt = 1:23
    for ct = 1:21
        %for d = 1:46
            for c = 2:90

                    
                    if isnan(hno3_obs(cnt,ct,c)) == 0
                        model_hno3(cnt,ct,c) = squeeze(mean_hno3(cnt,ct,c));
                    end
                    
                    if isnan(ans_obs(cnt,ct,c)) == 0
                        model_ans(cnt,ct,c) = squeeze(mean_ans_dissert(cnt,ct,c));
                    end
%                    if isnan(hno3_obs(cnt,ct,c)) == 0
%                        model_no2_sat(cnt,ct,c) = squeeze(mean_no2_sat(cnt,ct,c));
%                    end
%                    if isnan(hno3_obs(cnt,ct,c)) == 0
%                        model_oh(cnt,ct,c) = squeeze(mean_oh(cnt,ct,c));
%                    end
                    if isnan(o3_obs(cnt,ct,c)) == 0
                        model_o3(cnt,ct,c) = squeeze(mean_o3(cnt,ct,c));
                    end
                    if isnan(nox_obs(cnt,ct,c)) == 0
                        model_nox(cnt,ct,c) = squeeze(mean_nox(cnt,ct,c));
                    end
                    if isnan(no2_obs(cnt,ct,c)) == 0
                        model_no2(cnt,ct,c) = squeeze(mean_no2(cnt,ct,c));
                    end
%                    if isnan(h2o_obs(cnt,ct,c)) == 0
%                        model_h2o(cnt,ct,c) = squeeze(mean_h2o(cnt,ct,c));
%                    end
                    
%                    if isnan(hno3_obs(cnt,ct,c)) == 0
%                        model_n2o5(cnt,ct,c) = squeeze(mean_n2o5(cnt,ct,c));
%                    end
                    %model_press(cnt,ct,c) = squeeze(press(cnt,ct,d,c));
                    %c
                    %pause
                    
                %end
            %end
            end
    end
end


