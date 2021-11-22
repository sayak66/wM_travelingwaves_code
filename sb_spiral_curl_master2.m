%% Master code that inputs data files and outputs coefficient data

clearvars -except tin
clc
tic

for master = 1:2 %sessions
    clearvars -except tin master
    clc
    close all
    

    switch master
         case 1
            load('tin_Drake_0611_all.mat');
            fnm = 'D0611_justCurl_ALL_reg4_PN2_win20_corr_';
            reg = 4;
         case 2
            load('tin_Drake_0612_all.mat');
            fnm = 'D0612_justCurl_ALL_reg4_PN2_win20_corr_';
            reg = 4;
    end

for lc = [4 8 12]
    if(lc == 4)
        hc = 8;
    elseif(lc == 8)
        hc = 12;
    elseif(lc == 12)
        hc = 30;
    else
        hc = 120;
    end

clearvars -except tin lc hc master fnm fnmS reg asf
close all

mi_gap = 99;
mi_range = 1;

t_rng = 1:7000;


mi_counter = 1;
for mi = mi_range
    
    clearvars -except tin mi mi_gap mi_range rn0_gap rn_a t_rng mi_counter non_sp_store sp_P_store sp_PN_store wv_mc_sp_PN wv_mc_sp_P wv_mc_nonsp ...
                    cc_spT lc hc all_sp_store fnm fnmS reg all_waves_loc cp_T asf cc_pvT
    
    cc=1;
    tout = tin;
    for i = 1:length(tout)
        ids = tout(i,1).outcome{1};
            if(strcmp(ids,'correct')) % correct trials only
                tout_d(cc,:) = tout(i,reg);
                cc=cc+1;
            end
    end



    tim = t_rng;


        c_rng = 1;
        for rng1 = 1:length(tout_d)

            mi
            rng1
            
            % code to calculate circular-circular correlations
            [cc_spA,pv] = spiral_nonspiral_mark_each_use_cpALL(tout_d,t_rng,rng1,reg,lc,hc,4);
            cc_spT(mi_counter,c_rng,:,:) = cc_spA;
            cc_pvT(mi_counter,c_rng,:,:) = pv;

            c_rng = c_rng + 1;
        end
    mi_counter = mi_counter + 1;
end
%%
fnme = strcat(fnm,'_',num2str(lc),'-',num2str(hc),'Hz');
save(fnme,'mi_range','mi_gap','cc_spT','cc_pvT','-v7.3');


end

end

toc
