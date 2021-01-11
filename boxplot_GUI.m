function [ output_args ] =boxplot_GUI(Z,foname,zname,handles,axestouse)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L=length(Z(1,:));
if sum(~isnan(Z(:,end)))<5
    L=L-1;
end
Maxappeartime=handles.minmaps;
density_th=handles.cutoffdensity;
eventwindow=handles.eventwindow;
for i=1:L-1
    for j=i+1:L
        
       if lillietest(Z(:,j))+lillietest(Z(:,i))>1
           [p,h]=ranksum(Z(:,i),Z(:,j));
       else [h,p]=ttest(Z(:,i),Z(:,j));
       end
       C(i,j)=h;
       P(i,j)=p;
       
    end
end
set(gcf,'CurrentAxes',axestouse);
cla;
boxplot(Z)

title(zname)

foldername=fullfile(foname);
Tp=handles.avgTp;
        if (exist(foldername) == 0)
            mkdir (foldername);
        end

        if Tp==2
            typefolder='PTX';
        elseif Tp==3
            typefolder='APV_PTX';
        elseif Tp==1
            typefolder='highMg';
        elseif Tp==4
            typefolder='TTX';
        else
            typefolder=handles.additiondrugs;
        end
        f2=fullfile(foldername,typefolder);
        if (exist(f2) == 0)
            mkdir (f2);
        end
        f4=fullfile(f2,sprintf('AvgMap_disttopiaMin%iMax%iDens%i+%i',handles.mindist,handles.maxdist,Maxappeartime,density_th*100))
        if (exist(f4) == 0)
            mkdir (f4);
        end
   filesname1=fullfile(f4,sprintf('boxplot%seventwindow%iTp%iholdingPotential%i.eps',zname,eventwindow,handles.avgTp,handles.avgHoldingPotential));
    filesname2=fullfile(f4,sprintf('boxplot%seventwindow%iTp%iholdingPotential%i.fig',zname,eventwindow,handles.avgTp,handles.avgHoldingPotential));
    filesname3=fullfile(f4,sprintf('boxplot%seventwindow%iTp%iholdingPotential%i.png',zname,eventwindow,handles.avgTp,handles.avgHoldingPotential));
    saveImages=1;

     currfig= find_figure(zname);
     
        
       cla;
boxplot(Z)

title(zname)
 
print(currfig,'-depsc2',filesname1);
        saveas(currfig,filesname2,'fig');
        print(currfig,'-dpng',filesname3);
        close(currfig);


end

