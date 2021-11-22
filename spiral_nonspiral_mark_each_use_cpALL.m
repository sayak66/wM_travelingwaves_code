function [cc_sp,cc_pv] = spiral_nonspiral_mark_each_use_cpALL(tin,rng,tin_i,reg,lc,hc,filo)

tout = tin;
outw = tout(tin_i).outw;


outw = outw(rng,:,:);
outwN = outw(rng,:,:);
outwN(find(isnan(outwN))) = 0;
% Bandpass filter
% lc = 4;
% hc = 8;
Fs = 1000;
% filo = 4;
outw1 = permute(outwN,[2 3 1]);
out1 = bandpass_filter(outw1,lc,hc,filo,Fs);

% Hilbert transform
outa = analytic_signal(out1);  % X now contains the "analytic signal"
amp = abs( outa );  			  % a contains the "amplitude envelope" at each point in time
ph = angle( outa );  		  % p now contains the "phase maps"

ph1 = ph;
amp1 = amp;
% Filter
for i = 1:size(ph1,3)
    ph1s(:,:,i) = imgaussfilt(ph1(:,:,i));
end
% Interpolate missing electrodes
aid = find(isnan(permute(outw,[2 3 1])));
ph1s(aid) = NaN;
ph1(aid) = NaN;
for i = 1:size(ph1,3)
    ph1sa = squeeze(ph1s(:,:,i));
    ph1ma = fillmissing(ph1sa,'linear',2,'EndValues','nearest');
    ph1m(:,:,i) = ph1ma;
    
    ph1sa_unf = squeeze(ph1(:,:,i));
    ph1ma_unf = fillmissing(ph1sa_unf,'linear',2,'EndValues','nearest');
    ph1m_unf(:,:,i) = ph1ma_unf;
end

[pm,pd,dx,dy] = phase_gradient_complex_multiplication(outa,1,1);
mm = meshgrid(1:8,1:8);

cpM = [4 4; 1 4; 4 1; 4 8; 8 4]; % all choice points
for curli = 1:size(cpM,1)
    
cp = cpM(curli,:);
for tidx = rng
    ss = ph1m_unf(:,:,tidx);    
    aa = pd(:,:,tidx);
%     crl = curl(aa,mm);
aax = dx(:,:,tidx);
aay = dy(:,:,tidx);
aap = pm(:,:,tidx);

crl = curl(aap,mm);
[a,b] = max(crl(:));
[r,c] = ind2sub(size(crl),b);

    [cc,pv,cp1] = phase_correlation_rotation(ss,crl,cp,1);
    cc_sp(curli,tidx) = cc;
    cc_pv(curli,tidx) = pv;

end

end
