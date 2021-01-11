
% filter LSPS data
function fdata = LSPSfilter_LFP(header,data)

% apply low-pass filter to data
[b a] = butter(6,300/header.sampleRate,'low'); % frequency band is 250HZ?
fdata = filtfilt(b,a,data);
%h=fdesign.lowpass('N,F3dB',13,1000/header.sampleRate*2);
%d1 = design(h,'butter');
%fdata = filtfilt(d1.sosMatrix,d1.ScaleValues,data);

%f = [0 100/header.sampleRate*2 100/header.sampleRate*2 1];
%m = [1 1 0 0];
%[b,a] = yulewalk(8,f,m);
%fdata=filtfilt(b,a,data); %zero-phase filtering

% % move to mean of baseline
% fdata = bsxfun(@minus, fdata, mean(fdata));
% fdata = fdata+baseline;
 
% remove end points after filtering
% begs = fix(0.002*header.sampleRate);
% ends = fix(0.998*header.sampleRate);
% for i=1:header.nPts
%   fdata(1:begs,i) = baseline(1:begs,i);
%   fdata(ends:end,i) = baseline(ends:end,i);
% end
for i=1:1:3
    [b a]=butter(3, [i*60-1 i*60+1]*2/header.sampleRate,'stop');
    fdata = filtfilt(b,a,fdata);
end
return
