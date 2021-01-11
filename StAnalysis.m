function [ h,p] =StAnalysis( x1,x2,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%type=1, ranksum; type=0, kstest
if lillietest(x1)+lillietest(x2)>=1
    if type==1
        [p,h]=ranksum(x1,x2);
    elseif type==2
        [h,p]=kstest2(x1,x2,'Tail','smaller');
    else[h,p]=kstest2(x1,x2);
    end
else [h,p]=ttest2(x1,x2);
end
end

