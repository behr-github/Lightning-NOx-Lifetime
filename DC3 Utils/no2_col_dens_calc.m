%Calculates NO2 column density

calc_no2_sat(1:23,1:21,1:90)=0;

for cnt = 1:23;
    for ct = 1:21;
        for c = 1:90;
            for d = 1:47
                
                if isnan(no2_sat(cnt,ct,d,c)) == 1
                    no2_sat(cnt,ct,d,c) = 0;
                end
                
                calc_no2_sat(cnt,ct,c) = calc_no2_sat(cnt,ct,c) + (no2_sat(cnt,ct,d,c).*...
                    1e-9.*numden(cnt,ct,d,c).*bxheight(cnt,ct,d,c).*(100));
            end
        end
    end
end