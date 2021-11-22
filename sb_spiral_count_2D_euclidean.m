load('T0810_T0815_T0820_T0823_T0828_3freq_5curls_dPFC_all_left.mat');  %load coefficients
fnm = 'euclidean_wavetypecounts';

clear cc_sp_frT
cc_sp_frT(:,:,:,1) = cc_sp_fr44(:,:,:);
cc_sp_frT(:,:,:,2) = cc_sp_fr14(:,:,:);
cc_sp_frT(:,:,:,3) = cc_sp_fr41(:,:,:);
%
tic
clearvars -except cc_sp_frT a_44_14 cc_sp_fr44 cc_sp_fr14 cc_sp_fr41 cc_sp_fr48 fnm
clc
rn_frng = [1 501 1001 1501 2001 2501 3001 3501 4001 4501 5001 5501 6001 6501];

load('curl_planar01_sp24r_sp42r_sp24_sp42_noiseavg_fix_5curls'); %simulated data

a_cp3(:,1) = a44t;
a_cp3(:,2) = a84t;
a_cp3(:,3) = a48t;
a_44_14 = a_cp3;

ar = 0.1;

aout_b = zeros(3,size(cc_sp_frT,2),size(a_44_14,1));
aout_u = zeros(3,length(rn_frng),size(cc_sp_frT,2),size(a_44_14,1));


for fr = 1:3

    for tr = 1:size(cc_sp_frT,2)
        [fr tr]
    
        for i = 1:length(rn_frng)
            
            cc_use1 = squeeze(cc_sp_frT(fr,tr,rn_frng(i):rn_frng(i)+499,:));
            
            for j = 1:length(cc_use1)
                px = cc_use1(j,1);
                py = cc_use1(j,2);
                pz = cc_use1(j,3);
                if(px>0.3 || px<-0.3 || py>0.3 || py<-0.3 || pz>0.3 || pz<-0.3)
                    x1 = a_44_14(:,1);
                    y1 = a_44_14(:,2);
                    z1 = a_44_14(:,3);
                    ec = sqrt((x1-px).^2 + (y1-py).^2 + (z1-pz).^2);
                    idx = find(ec==min(ec));
                    aout_u(fr,i,tr,idx) = aout_u(fr,i,tr,idx) + 1;
                end
            end
        end
    end
end
    

toc

save(fnm,'aout_u','rn_frng');