%Takes a different, probably too long, approach to match observations to
%model results

q = 1;a = 1; b = 1;

for cnt = 1:22
    for ct = 1:20
        for c = 1:46
            for d = 1:90

                for cntr = 1:length(dc3_merges_ut(:,1));

                if dc3_merges_ut(cntr,2)-91 == d
                    if dc3_merges(cntr,9) < press(cnt,ct,c,d) && ...
                            dc3_merges(cntr,9) > press(cnt,ct,c+1,d)
                    if dc3_merges_ut(cntr,7)-360 > NA_lon(cnt) && ...
                            dc3_merges_ut(cntr,7)-360 < NA_lon(cnt+1)
                        if dc3_merges_ut(cntr,6) > NA_lat(ct) && ...
                                dc3_merges_ut(cntr,6) < NA_lat(ct+1)

                            if isnan(dc3_merges_ut(cntr,77)) == 0

                                model_results_hno3(a,1) = hno3(cnt,ct,c,d);
                                observ_results_hno3(a,1) = dc3_merges_ut(cntr,77);

                                a = a + 1;

                            end

                            if isnan(dc3_merges_ut(cntr,71).*1e3+dc3_merges_ut(cntr,73))==0

                                model_results_nox(b,1) = nox(cnt,ct,c,d);
                                observ_results_nox(b,1) = dc3_merges_ut(cntr,71).*1e3+dc3_merges(cntr,73);

                                b = b + 1;

                            end

                            if isnan(dc3_merges_ut(cntr,72))==0

                                model_results_o3(q,1) = o3(cnt,ct,c,d);
                                observ_results_o3(q,1) = dc3_merges_ut(cntr,72);

                                q = q + 1;
                                %d
                            end
                        end
                    end
                    end
                end
                end

            end
        end
    end
end