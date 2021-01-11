function [ output_args ] =PercentDist(folderMeanpath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a=dir(fullfile(folderMeanpath,'*.mat'));
dx=40;
xdist=0:dx:1000;


for ii=1:size(a)
    filename=fullfile(folderMeanpath,a(ii).name);
    Cg= load(filename);
    
    if isfield(Cg,'Cellrecord')
        Cellrecord=Cg.Cellrecord;
        Num_cell=Cg.Num_cell;
        for fn=1:Cg.Num_cell+1
            
            flagD=Cellrecord{fn}.meanflagDttxInd60;
            
            ind=find(flagD);
            xaxis=Cellrecord{fn}.StimCoordinates(1,ind);
            yaxis=Cellrecord{fn}.StimCoordinates(2,ind);
            xsoma=Cellrecord{fn}.SomaCoordinates(1,1);
            ysoma=Cellrecord{fn}.SomaCoordinates(1,2);
            
            dist=sqrt((xaxis-xsoma).^2+(yaxis-ysoma).^2);
            
            Totalevent=length(ind);
            
            percD=zeros(length(xdist),1);
            for i=1:length(xdist)
                percD(i)=sum(dist<=xdist(i))/Totalevent*100;
            end
            Cellrecord{fn}.diameter=xdist;
            Cellrecord{fn}.percD=percD;
            
            
            flag=Cellrecord{fn}.meanflagttxInd60;
            
            ind=find(flag);
            xaxis=Cellrecord{fn}.StimCoordinates(1,ind);
            yaxis=Cellrecord{fn}.StimCoordinates(2,ind);
            xsoma=Cellrecord{fn}.SomaCoordinates(1,1);
            ysoma=Cellrecord{fn}.SomaCoordinates(1,2);
            
            dist=sqrt((xaxis-xsoma).^2+(yaxis-ysoma).^2);
            
            Totalevent=length(ind);
            
            perc=zeros(length(xdist),1);
            for i=1:length(xdist)
                perc(i)=sum(dist<=xdist(i))/Totalevent*100;
            end
            %       Cellrecord{fn}.diameter=xdist;
            Cellrecord{fn}.perc=perc;
            
            h=find_figure(sprintf('cdf_EventPerDistance%i',ii))
            subplot(5,1,fn)
            plot(xdist,perc,'r')
            hold on
            plot(xdist,percD,'k--')
            title(Cellrecord{fn}.MapnameflagttxInd60)
             
             
 
        end
         xlabel('distance')
    ylabel('%')
                    saveas(h,fullfile(folderMeanpath,sprintf('CdfPercDist%i_%i.fig',60,ii)),'fig');
print(h,'-depsc2',fullfile(folderMeanpath,sprintf('CdfPercDist%i_%i.eps',60,ii)));
    end
    Cg.Cellrecord=Cellrecord;
    save(filename,'Cellrecord','Num_cell')
    
    
end
end

