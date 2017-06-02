%removes results in the stratosphere & below prescribed pressure levels


%removes results above tropopause
for cnt = 1:23
    for ct = 1:21
        for d = 1:90
            hno3(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
            nox(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
            o3(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
            no2(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
            no2_sat(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
            switch what_to_compare
                case 'Dissert_33'
                    ans_dissert(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
                case 'fk_ans'
                    ans_dissert(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
            end
%            h2o(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
%            n2o5(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
%            no2_sat(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
%            oh(cnt,ct,floor(trop_hght(cnt,ct,d)):end,d)=NaN;
                
        end
    end
end



to3=find(o3 > 150);nox(to3)=NaN;hno3(to3)=NaN;o3(to3)=NaN;
%tnox = find(nox./hno3 > 5);nox(tnox)=NaN;hno3(tnox)=NaN;o3(tnox)=NaN;    
%removes results not in ut
for cnt = 1:23
    for ct = 1:21
        for c = 1:47
            for d = 1:90
                
                if press(cnt,ct,c,d) > 350
                    hno3(cnt,ct,c,d) = NaN;
                    nox(cnt,ct,c,d) = NaN;
                    o3(cnt,ct,c,d) = NaN;
                    no2(cnt,ct,c,d) = NaN;
                    switch what_to_compare 
                        case 'Dissert_33'
                            ans_dissert(cnt,ct,c,d)=NaN;
                        case 'fk_ans'
                            ans_dissert(cnt,ct,c,d)=NaN;
                    end
%                    h2o(cnt,ct,c,d) = NaN;
%                    n2o5(cnt,ct,c,d) = NaN;
%                    no2_sat(cnt,ct,c,d) = NaN;
%                    oh(cnt,ct,c,d) = NaN;
                end
                if press(cnt,ct,c,d) < 200
                    hno3(cnt,ct,c,d) = NaN;
                    nox(cnt,ct,c,d) = NaN;
                    o3(cnt,ct,c,d) = NaN;
                    no2(cnt,ct,c,d) = NaN;
                    switch what_to_compare 
                        case 'Dissert_33'
                            ans_dissert(cnt,ct,c,d)=NaN;
                        case 'fk_ans'
                            ans_dissert(cnt,ct,c,d)=NaN;
                    end
 %                   h2o(cnt,ct,c,d) = NaN;
 %                   n2o5(cnt,ct,c,d) = NaN;
 %                   no2_sat(cnt,ct,c,d) = NaN;
 %                   oh(cnt,ct,c,d) = NaN;
                end
                

                
            end
        end
    end
end
%hi = 1
mean_nox = squeeze(nanmean(nox,3));
mean_hno3 = squeeze(nanmean(hno3,3));
mean_o3 = squeeze(nanmean(o3,3));
mean_no2 = squeeze(nanmean(no2,3));
switch what_to_compare
    case 'Dissert_33'
        mean_ans_dissert = squeeze(nanmean(ans_dissert,3));
    case 'Base'
        mean_ans_dissert(1:23,1:21,1:90)=NaN;
    case 'MPN'
        mean_ans_dissert(1:23,1:21,1:90)=NaN;
    case 'PNA'
        mean_ans_dissert(1:23,1:21,1:90)=NaN;
    case 'HNO3'
        mean_ans_dissert(1:23,1:21,1:90)=NaN;
    case 'Dissert'
        mean_ans_dissert(1:23,1:21,1:90)=NaN;
    case 'N2O5'
        mean_ans_dissert(1:23,1:21,1:90)=NaN;
    case 'fk_ans'
        mean_ans_dissert = squeeze(nanmean(ans_dissert,3));
end
%mean_h2o = squeeze(nanmean(h2o,3));
%mean_n2o5 = squeeze(nanmean(n2o5,3));
%mean_no2_sat = squeeze(nanmean(no2_sat,3));
%mean_oh = squeeze(nanmean(oh,3));