function [ output_args ] = AverageSpikeCount(folderpath1,folderpath2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
eventwindow=50;
SN1 = SpikeCount(folderpath1,eventwindow);
SN2 = SpikeCount(folderpath2,eventwindow);
SN1=deleteoutliers(SN1);
SN2=deleteoutliers(SN2);
[h,p]=StAnalysis(SN1,SN2,2);
figure()
cdfplot(SN1)
hold on
cdfplot(SN2)

title(sprintf('p=%i',p))
hold off
end

