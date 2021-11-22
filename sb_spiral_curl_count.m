%% Analyze the coefficient data obtained from master code

tic

clearvars -except tin
clc

ci = 1;
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


for rec = 1:3

clearvars -except tin cc_sp_fr44 cc_sp_fr14 cc_sp_fr41 cc_sp_fr48 rec cc_sp_fr ...
    ccsp1 ccsp2 ccsp3 ccsp4 ccsp5 lc hc ci rec cc_sp_fr_t3 cc_sp_fr_t1 cc_sp_fr_d2
% clc 
        
    [rec ci]
    
% load all coefficient data from master code - all sessions
load(strcat('D0611_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cc_spT(:,:,rec,:);
load(strcat('D0612_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0731_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0605_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0517_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0515_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0807_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0814_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));
load(strcat('D0816_justCurl_ALL_reg4_PN2_win20_corr_',num2str(lc),'-',num2str(hc),'Hz.mat'));
cc_spt1 = cat(2,cc_spt1,cc_spT(:,:,rec,:));


cc_spT = cc_spt1;

cc_sp_store = cc_spT;


%%
cc_sp_fr(rec,ci,:,:) = cc_sp_store;

end

ci = ci + 1;

end

toc
%% The three coefficients

cc_sp_fr44 = squeeze(cc_sp_fr(1,:,:,:));
cc_sp_fr14 = squeeze(cc_sp_fr(2,:,:,:));
cc_sp_fr41 = squeeze(cc_sp_fr(3,:,:,:));

%% Plot trends
figure;
fr = 3;
nbe = -1:0.01:1;
idx = randi([1 5000]);
cuse = squeeze(cc_sp_fr44(fr,:,:));
rng = [500 1500 2000 2500 3000 3500 4000 4500 5000];
for i = 1:length(rng)
    bas = cuse(:,1:401);
    use = cuse(:,rng(i):(rng(i)+400));
    subplot(3,3,i)
    histogram(bas,nbe,'Normalization','probability','FaceColor',[1 0 0]);
    hold on;
    histogram(use,nbe,'Normalization','probability','FaceColor',[0 0 1]);
    hold on;
    xlim([-1 1]);ylim([0 0.015]);
    title(num2str(rng(i)));
end
%% Histcounts and signficance
cc_sp_fr = cc_sp_fr44; %coefficient to check

thr = 0.3;
nb1 = 50;
nbe = -1:0.04:1;
ymax = 16e4;
rn_frng = [1 1001 1501 2001 2501 3001 3501 4001 4501 5001 5501 6001 6500];
ct1 = 3;
ct2 = length(rn_frng);
rni = 1;
nfr = 3;

tr_rng = 1:size(cc_sp_fr,2);

n1 = zeros(nfr,length(tr_rng),length(rn_frng),length(nbe)-1);
b1 = zeros(nfr,length(tr_rng),length(rn_frng),length(nbe));
n2 = zeros(nfr,length(tr_rng),length(rn_frng),length(nbe)-1);
b2 = zeros(nfr,length(tr_rng),length(rn_frng),length(nbe));
n3 = zeros(nfr,length(tr_rng),length(rn_frng),length(nbe)-1);
b3 = zeros(nfr,length(tr_rng),length(rn_frng),length(nbe));
chsq_n = zeros(nfr,length(tr_rng),length(rn_frng));


for fr = 1:3
ctr = 1;
for tr = tr_rng

cc_bas = cc_sp_fr(fr,tr,501:1001);
cc_bas1 = cc_bas(:);

cc_sam = cc_sp_fr(fr,tr,2001:2501);
cc_sam1 = cc_sam(:);
ci = 1;

for rn_fr = rn_frng
    cc_sc = cc_sp_fr(fr,tr,rn_fr:(rn_fr+500));
    cc_sc1 = cc_sc(:);

    [n1(fr,ctr,ci,:),b1(fr,ctr,ci,:)] = histcounts(cc_bas1,nbe);
    [n2(fr,ctr,ci,:),b2(fr,ctr,ci,:)] = histcounts(cc_sc1,nbe);
    [n3(fr,ctr,ci,:),b3(fr,ctr,ci,:)] = histcounts(cc_sam1,nbe);

    bb = squeeze(b1(1,1,1,1:end-1));
    rng_bb = [1:find(bb<-thr,1,'last') find(bb>thr,1,'first'):length(bb)];
    rng_bbP = find(bb>thr,1,'first'):length(bb);
    rng_bbN = 1:find(bb<-thr,1,'last');

    chsq_n0 = ((n2(fr,ctr,ci,rng_bb) - n1(fr,ctr,ci,rng_bb)).^2)./(n1(fr,ctr,ci,rng_bb)+n2(fr,ctr,ci,rng_bb));
    chsq_n0(find(isnan(chsq_n0))) = 0;
    chsq_n(fr,ctr,ci) = sum(chsq_n0,4);

    chsq_n0S = ((n2(fr,ctr,ci,rng_bb) - n3(fr,ctr,ci,rng_bb)).^2)./(n3(fr,ctr,ci,rng_bb)+n2(fr,ctr,ci,rng_bb));
    chsq_n0S(find(isnan(chsq_n0S))) = 0;
    chsq_nS(fr,ctr,ci) = sum(chsq_n0S,4);

    chsq_n0P = ((n2(fr,ctr,ci,rng_bbP) - n1(fr,ctr,ci,rng_bbP)).^1)./(n1(fr,ctr,ci,rng_bbP)+n2(fr,ctr,ci,rng_bbP));
    chsq_n0P(find(isnan(chsq_n0P))) = 0;
    chsq_nP(fr,ctr,ci) = sum(chsq_n0P,4);

    chsq_n0N = ((n2(fr,ctr,ci,rng_bbN) - n1(fr,ctr,ci,rng_bbN)).^1)./(n1(fr,ctr,ci,rng_bbN)+n2(fr,ctr,ci,rng_bbN));
    chsq_n0N(find(isnan(chsq_n0N))) = 0;
    chsq_nN(fr,ctr,ci) = sum(chsq_n0N,4);

    ci = ci + 1;
end
ctr = ctr + 1;
end
end
%
bb = squeeze(b1(1,1,1,1:end-1));
figure;
ci = 1;
for fr = 1:3
for i = 1:size(n1,3)
subplot(nfr,size(n1,3),ci);
np1 = squeeze(n1(fr,:,i,:));
np2 = squeeze(n2(fr,:,i,:));
%     plot(bb,mean(np1),'r');
stdshade(np1,size(np1,1),0.4,'r',bb);
hold on;
stdshade(np2,size(np1,1),0.4,'g',bb);
%     plot(bb,mean(np2),'g');
hold on;
rng_bb = [1:find(bb<-thr,1,'last') find(bb>thr,1,'first'):length(bb)];
for ij = rng_bb
d1 = np1(:,ij);
d2 = np2(:,ij);
if(~isnan(ttest(d1,d2)) && ttest(d1,d2))
    hold on
    if(fr == 1)
        plot(bb(ij),30,'b*');
    elseif(fr == 2)
        plot(bb(ij),30,'b*');
    else
        plot(bb(ij),30,'b*');
    end

end
if(fr == 1)
    ylim([0 30]);
elseif(fr == 2)
    ylim([0 30]);
else
    ylim([0 30]);
end
end
hold on;
xline(thr);hold on;xline(-thr);
title(strcat('rng-',num2str(rn_frng(i))));
%     set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

ci = ci + 1;
end

end
%
figure;
for i = 1:3
subplot(3,1,i)
ci = squeeze(chsq_nP(i,:,:));
m1 = mean(ci,1);
s1 = std(ci,0,1)/(size(ci,1).^0.5)';
x = 1:size(ci,2);
y = m1;
err = s1;
errorbar(x,y,err,'-ro');
hold on
ci = squeeze(chsq_nN(i,:,:));
m1 = mean(ci,1);
s1 = std(ci,0,1)/(size(ci,1).^0.5)';
x = 1:size(ci,2);
y = m1;
err = s1;
errorbar(x,y,err,'-bo');
hold on; yline(0);
ylim([-2 2]);
for ti = 1:size(chsq_nP,3)
d1 = squeeze(chsq_nP(i,:,ti));
d2 = squeeze(chsq_nN(i,:,ti));
if(~isnan(ttest(d1)) && ttest(d1))
    plot(ti,-1.9,'r*');
end
if(~isnan(ttest(d2)) && ttest(d2))
    plot(ti,-1.6,'b*');
end
end

xticks(x);
xticklabels({'b1','b2','f','s','d1','d2','d3','d4','t','t1','t2','t3','t4'});
end
