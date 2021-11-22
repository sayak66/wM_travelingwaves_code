%% Find spike rates
clearvars -except cc_sp_fr14 cc_sp_fr44 cc_sp_fr41 phT
clc

tic


for i = 1

filen = strcat('Tiergan-DMTS-20180810'); %load session

mf = matfile(strcat(fnm,filen));


ui(i) = mf.unitInfo;
eeD(i) = mf.electrodeInfo;
ar(i).dat = ui(i).area;

aloc1 = ui(i).gridLoc;
spT = mf.spikeTimes; % spike times

    ttc = tt(i).corr;
    ci = 1;
    sp_out = zeros(length(ttc),size(spT,2),7001);
    for k = ttc
        cj = 1;
        for j = 1:size(spT,2)
            t1 = spT{k,j};
            t2 = round(t1*1000) + 2001;
            sp_out(ci,cj,t2) = 1; % 1 if spike at a particular instant
            cj = cj + 1;
        end
        ci = ci + 1;
    end
    
    sp_total = sp_out;


%electrode numbers to analyze
st = 157;
en = 228;


u1 = ui(i).unitID;
cut = 1;
clear use
for ut = st(i):en(i)
    u2 = u1{ut};
    u3 = u2(end-1:end);
    if(~strcmp(u3,'00'))
        use(cut) = ut;
        cut = cut + 1;
    end
end


    a1 = sp_total;
    a2 = a1(:,use,:);
    a22 = squeeze(mean(a2,2));%/1e-3; %spike rate
    wn = 100;
    a3 = 0*a22;
    for ai = 1:(size(a22,2)-wn)
        ai1 = a22(:,ai:(ai+wn));
        ai2 = sum(ai1,2)/(wn*1e-3);
        a3(:,ai) = ai2;
    end
    if(i==1)
        sp1 = a3;
    else
        sp1 = cat(1,sp1,a3);
    end
    
end

toc
%% Spikes vs wave type
tic
clearvars -except sp1 phT cc_sp_fr44 cc_sp_fr14 cc_sp_fr41 ...
    sp_check_D2_beta sp_check_D2_theta sp_check_D4_beta sp_check_D4_theta ...
    sp_check_T1_beta sp_check_T1_theta sp_check_T3_beta sp_check_T3_theta
clc

load('curl_planar01_sp24r_sp42r_sp24_sp42_noiseavg_fix_5curls');

a_cp3(:,1) = a44t;
a_cp3(:,2) = a14t;
a_cp3(:,3) = a41t;
a_44_14 = a_cp3;

fr = 3;
c14 = squeeze(cc_sp_fr14(fr,:,:));
c41 = squeeze(cc_sp_fr41(fr,:,:));
c44 = squeeze(cc_sp_fr44(fr,:,:));
sp_planar = nan*sp1;
sp_rotate = nan*sp1;

sp_check = nan*zeros(7764693,72);
ct = ones(1,72);
% max-ct T3 fr3 - 6064693
% max-ct T1 fr3 - 7005370
% sp_check = zeros(1,72);
% ct = zeros(1,72);

for tr = 1:size(c14,1)
    tr
    for i = 1:size(c14,2)
        px = c44(tr,i);
        py = c14(tr,i);
        pz = c41(tr,i);
        if(px>0.3 || px<-0.3 || py>0.3 || py<-0.3 || pz>0.3 || pz<-0.3)
            x1 = a_44_14(:,1);
            y1 = a_44_14(:,2);
            z1 = a_44_14(:,3);
            ec = sqrt((x1-px).^2 + (y1-py).^2 + (z1-pz).^2);
            idx = find(ec==min(ec));
%             sp_check(idx) = sp_check(idx) + sp1(tr,i);
            sp_check(ct(idx),idx) = sp1(tr,i);
            ct(idx) = ct(idx) + 1;
            if(idx<=32)
                sp_planar(tr,i) = sp1(tr,i);
            else
                sp_rotate(tr,i) = sp1(tr,i);
            end
        end
    end
end
figure;
% plot(sp_check./ct);
hold on;xline(33);hold on;xline(38);
hold on;xline(43);hold on;xline(48);
hold on;xline(53);hold on;xline(58);
hold on;xline(63);hold on;xline(68);hold on;plot(100*ct/sum(ct),'r');
x = 1:72;
y = nanmean(sp_check,1);
err = nanstd(sp_check,0,1)./(ct.^0.5);
hold on;
errorbar(x,y,err);
figure;
dwn = 1;
a1 = squeeze(sp_rotate(:,1:dwn:5500));
c1 = squeeze(sp_planar(:,1:dwn:5500));
subplot(1,2,1);
stdshade(sp1,size(sp1,1),0.4,'r');
hold on;yline(1.3);
ylim([0 8]);
subplot(1,2,2);
stdshade(a1,size(a1,1),0.4,'r');
hold on;
stdshade(c1,size(c1,1),0.4,'g');
hold on;yline(1.3);
%% Plot spikes vs high/low rotation wavelengths
clc

aa = sp_check;

kn = [33 38 43 48 53 58 63 68];
lwv = sort([kn kn+1]); % low wavelength waves
hwv = sort([kn+3 kn+4]); % high wavelength waves

aLWV = aa(:,lwv);
aLWV(find(isnan(aLWV))) = [];
aHWV = aa(:,hwv);
aHWV(find(isnan(aHWV))) = [];

figure;
x = 1:2;
y = [mean(aLWV) mean(aHWV)];
ys = [std(aLWV)/(length(aLWV)^0.5) std(aHWV)/(length(aHWV)^0.5)];
bar(x,y);
hold on;
er = errorbar(x,y,ys);
er.Color = [0 0 0];
er.LineStyle = 'none';

%% Smooth and plot (norm or otherwise)
clear sp
figure;
sp(1).dat = sp1;

norm = 1;
for i = 1%:4
    sp1 = sp(i).dat;
    sp2 = 0*sp1;
    for j = 1:size(sp1,1)
        [i j]
        sp2(j,:) = smooth(sp1(j,:),0.01);
    end
    rng = 1500:4550;
    bas = 1750:1950;
    if(norm)
        bb = mean(sp2(:,bas),2);
        bs = std(sp2(:,bas),0,2);
%         sp3 = (sp2-bb)./bb;
        sp3 = (sp2-bb);%./bs;
    else
        sp3 = sp2;
    end
%     aa = ind2sub([size(sp3,1) size(sp3,2)],find(sp3==Inf));
%     sp3(aa) = nan;
    
    subplot(2,2,i)
    stdshade(sp3,size(sp3,1),0.4,'r');
    hold on;
    xline(2000);hold on;xline(4500);
    hold on;yline(0);
    xlim([rng(1) rng(end)]);
%     ylim([1 4]);
end
