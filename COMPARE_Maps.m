function [ output_args ] = COMPARE_Maps(Control_path, Object_path,savepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%type=1; %use ranksum
folderMeanpath{1}=Control_path;
folderMeanpath{2}=Object_path;






for cmpind=1:2
    a=dir(fullfile(folderMeanpath{cmpind},'*.mat'));
    Cellfile=size(a);
    
    DirArea=[];
    Ag=[];
    Bnd=[];
    n=0;
    N=[];
    
    Xaxis=cell(Cellfile(1,1).*5,1);
    Yaxis=cell(Cellfile(1,1).*5,1);
    Zpk60avg=cell(Cellfile(1,1).*5,1);
    Zchg60avg=cell(Cellfile(1,1).*5,1);
    Zden60avg=cell(Cellfile(1,1).*5,1);
    Zpk10avg=cell(Cellfile(1,1).*5,1);
    Zchg10avg=cell(Cellfile(1,1).*5,1);
    Zden10avg=cell(Cellfile(1,1).*5,1);
    Zlat10avg=cell(Cellfile(1,1).*5,1);
    Zlat60avg=cell(Cellfile(1,1).*5,1);
    
    ZchgD60avg=cell(Cellfile(1,1).*5,1);
    ZdenD60avg=cell(Cellfile(1,1).*5,1);
    
    ZlatD60avg=cell(Cellfile(1,1).*5,1);
    for i= 1:1:Cellfile(1,1)
        Cg= load(fullfile(folderMeanpath{cmpind},a(i).name));
        
        
        if isfield(Cg,'Cellrecord')
            Cellrecord=Cg.Cellrecord;
            for fn=1:Cg.Num_cell+1
                
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')||isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
                        NL=length(Cellrecord{fn}.DistWidthPeakTp1Hp60);
                    else NL=length(Cellrecord{fn}.DistWidthPeakTp1Hp10);
                    end
                    Ag=[Ag;Cellrecord{fn}.age];
                    Bnd=[Bnd;Cellrecord{fn}.Boundry];
                    n=n+1;
                    
                    Soma= Cellrecord{fn}.SomaCoordinates;
                    flipimg=Cellrecord{fn}.flipimg;
                    pth=char(Cellrecord{fn}.Pth);
                    stimcoordinates=Cellrecord{fn}.StimCoordinates;
                    rotate_angle=Cellrecord{fn}.SpatialRotation;
                    [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                    X=stimCoordinates(1,:);
                    Y=stimCoordinates(2,:);
                    Xaxis{n}=X;
                    Yaxis{n}=Y;
                    Zpk60=zeros(size(X)).*NaN;
                    Zchg60=zeros(size(X)).*NaN;
                    Zden60=zeros(size(X));
                    Zpk10=zeros(size(X)).*NaN;
                    Zden10=zeros(size(X));
                    Zchg10=zeros(size(X)).*NaN;
                    Zlat60=zeros(size(X)).*NaN;
                    Zlat10=zeros(size(X)).*NaN;
                    ZlatD60=zeros(size(X)).*NaN;
                    ZchgD60=zeros(size(X)).*NaN;
                    ZdenD60=zeros(size(X));
                    
                    ZlatD10=zeros(size(X)).*NaN;
                    
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
                        
                        LdistPk60{n}=Cellrecord{fn}.LaydistPkTp1Hp60;
                        LdistChg60{n}=Cellrecord{fn}.LaydistChgTp1Hp60;
                        
                        Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
                        DirArea=[DirArea,length(Indden).*Cellrecord{fn}.PatternSpacing(1).^2];
                        Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg60));
                        Zpk60(Indpk)=Cellrecord{fn}.meanpeakHighMg60(Indpk);
                        Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg60));
                        Zchg60(Indchg)=Cellrecord{fn}.meanareaHighMg60(Indchg);
                        Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg60));
                        Zlat60(Indlat)=Cellrecord{fn}.meanlatencyHighMg60(Indlat);
                        Indden=find((Cellrecord{fn}.meanflagHighMg60>0));
                        Zden60(Indden)=Cellrecord{fn}.meanflagHighMg60(Indden);
                        
                        
                        Indchg=find(~isnan(Cellrecord{fn}.meanDirareaHighMg60));
                        ZchgD60(Indchg)=Cellrecord{fn}.meanDirareaHighMg60(Indchg);
                        Indlat=find(~isnan(Cellrecord{fn}.meanDirlatencyHighMg60));
                        ZlatD60(Indlat)=Cellrecord{fn}.meanDirlatencyHighMg60(Indlat);
                        Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
                        ZdenD60(Indden)=Cellrecord{fn}.meanflagDHighMg60(Indden);
                        
                        
                        
                    else Nanmaxtric= ones(1,NL).*NaN;
                        
                        LdistPk60{n}=NaN;
                        LdistChg60{n}=NaN;
                    end
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
                        
                        LdistPk10{n}=Cellrecord{fn}.LaydistPkTp1Hp10;
                        LdistChg10{n}=Cellrecord{fn}.LaydistChgTp1Hp10;
                        Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg10));
                        Zpk10(Indpk)=Cellrecord{fn}.meanpeakHighMg10(Indpk);
                        Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg10));
                        Zchg10(Indchg)=Cellrecord{fn}.meanareaHighMg10(Indchg);
                        Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg10));
                        Zlat10(Indlat)=Cellrecord{fn}.meanlatencyHighMg10(Indlat);
                        Indden=find((Cellrecord{fn}.meanflagHighMg10>0));
                        Zden10(Indden)=Cellrecord{fn}.meanflagHighMg10(Indden);
                        
                        
                    else
                        
                        LdistPk10{n}=NaN;
                        LdistChg10{n}=NaN;
                    end
                    
                    Zpk60avg{n}=Zpk60;
                    Zchg60avg{n}=Zchg60;
                    Zden60avg{n}=Zden60;
                    Zpk10avg{n}=Zpk10;
                    Zchg10avg{n}=Zchg10;
                    Zden10avg{n}=Zden10;
                    Zlat10avg{n}=Zlat10;
                    Zlat60avg{n}=Zlat60;
                    
                    ZchgD60avg{n}=ZchgD60;
                    ZdenD60avg{n}=ZdenD60;
                    
                    ZlatD60avg{n}=ZlatD60;
                    
                end
                
                
                
            end
        end
    end
    Cs{cmpind}.Zpk60avg=Zpk60avg;
    Cs{cmpind}.Zchg60avg=Zchg60avg;
    Cs{cmpind}.Zden60avg=Zden60avg;
    Cs{cmpind}.Zpk10avg=Zpk10avg;
    
    Cs{cmpind}.Zchg10avg=Zchg10avg;
    Cs{cmpind}.Zden10avg=Zden10avg;
    Cs{cmpind}.Zlat10avg=Zlat10avg;
    Cs{cmpind}. Zlat60avg= Zlat60avg;
    
    
    Cs{cmpind}.ZchgD60avg=ZchgD60avg;
    Cs{cmpind}.ZdenD60avg=ZdenD60avg;
    
    Cs{cmpind}.ZlatD60avg=ZlatD60avg;
    Cs{cmpind}.Yaxis=Yaxis;
    Cs{cmpind}.DirArea=DirArea;
    
    Cs{cmpind}.Xaxis=Xaxis;
    Cs{cmpind}.Bnd=Bnd;
    
    Cs{cmpind}.LdistPk60=LdistPk60;
    Cs{cmpind}.LdistPk10=LdistPk10;
    Cs{cmpind}.LdistChg60=LdistChg60;
    Cs{cmpind}.LdistChg10=LdistChg10;
    N=[N;n]
    
    
end

for objn=1:2
    Ncell=N(objn);
    
    
    dens=min(10,0.8*Ncell);
    Xavg=[];
    Yavg=[];
    dx=30;
    dy=30;
    Zpk60=[];
    Zchg60=[];
    Zden60=[];
    Zpk10=[];
    Zchg10=[];
    Zden10=[];
    Zlat60=[];
    Zlat10=[];
    
    ZpkD60=[];
    ZchgD60=[];
    ZdenD60=[];
    
    ZlatD60=[];
    
    Xaxis=Cs{objn}.Xaxis;
    Yaxis=Cs{objn}.Yaxis;
    Zpk60avg=Cs{objn}.Zpk60avg;
    Zchg60avg=Cs{objn}.Zchg60avg;
    Zden60avg=Cs{objn}.Zden60avg;
    Zlat60avg=Cs{objn}.Zlat60avg;
    Zpk10avg=Cs{objn}.Zpk10avg;
    Zchg10avg=Cs{objn}.Zchg10avg;
    Zden10avg=Cs{objn}.Zden10avg;
    Zlat10avg=Cs{objn}.Zlat10avg;
    ZchgD60avg=Cs{objn}.ZchgD60avg;
    ZdenD60avg=Cs{objn}.ZdenD60avg;
    ZlatD60avg=Cs{objn}.ZlatD60avg;
    Bnd=Cs{objn}.Bnd;
    
    
    xlim = [-330-0.1 960];
    ylim = [-480-0.1 480+0.1];
    
    xtick = xlim(1):dx:xlim(2);
    ytick = ylim(1):dy:ylim(2);
    xi = xtick(1:end-1)+dx/2;
    yi = ytick(1:end-1)+dy/2;
    xii=xi(1):1:xi(end);
    yii=yi(1):1:yi(end);
    for i=1:Ncell

        
        
        [yi1,xi1]=meshgrid(yi,xi);
        [yii,xii]=meshgrid(yii,xii);
        
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zden60avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('Singlecelldensity60_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        hold on
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        
        for ib=1:length(Bnd(1,:))-1
            ind{ib}=find(xii(:,1)>=Bnd(i,ib)&xii(:,1)<Bnd(i,ib+1));
            densM60{ib}=avg2(ind{ib},:);
            Xmax=max(xii(ind{ib},1));
            Xmin=min(xii(ind{ib},1));
            xcentre=(Xmax+Xmin)/2;
            Xax{ib}=(xii(ind{ib},:)-xcentre)./abs(Bnd(i,ib)-Bnd(i,ib+1)).*mean(abs(Bnd(:,ib)-Bnd(:,ib+1)));
            Yax{ib}=yii(ind{ib},:);
            
        end
        
        
        
        
        %     colormap(winter)
        colorbar
        %------------------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zpk60avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('Singlecellpeak60_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
        colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            PkM60{ib}=avg2(ind{ib},:);
            
        end
        
        %---------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zchg60avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('SinglecellChg60_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
        colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            ChgM60{ib}=avg2(ind{ib},:);
            
        end
        %---------------------------------------------------------------
        avg = gridavg2(Xaxis{i},Yaxis{i},Zlat60avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('SinglecellLatency60_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
        colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            LatM60{ib}=avg2(ind{ib},:);
            
        end
        %__________________________________________________________________
        
        
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zden10avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('Singlecelldensity10_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        hold on
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        
        for ib=1:length(Bnd(1,:))-1
            
            densM10{ib}=avg2(ind{ib},:);
            
            
        end
        
        
        
        
        %     colormap(winter)
        colorbar
        %------------------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zpk10avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('Singlecellpeak60_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
        colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            PkM10{ib}=avg2(ind{ib},:);
            
        end
        
        %---------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zchg10avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('SinglecellChg60_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
        colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            ChgM10{ib}=avg2(ind{ib},:);
            
        end
        %---------------------------------------------------------------
        avg = gridavg2(Xaxis{i},Yaxis{i},Zlat60avg{i},xtick,ytick,false,Ncell,0);
        
        currfig= find_figure(sprintf('SinglecellLatency10_%d',objn));
        subplot(10,6,i)
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
        axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
            'DataAspectRatio',[1 1 1]);
        
        imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
        colorbar
        
        for ib=1:length(Bnd(1,:))-1
            
            LatM10{ib}=avg2(ind{ib},:);
            
        end
        
        %__________________________________________________________________
        
        
        DensMLayer60{i}=densM60;
        PeakMLayer60{i}=PkM60;
        ChgMLayer60{i}=ChgM60;
        LatMlayer60{i}=LatM60;
         DensMLayer10{i}=densM10;
        PeakMLayer10{i}=PkM10;
        ChgMLayer10{i}=ChgM10;
        LatMlayer10{i}=LatM10;
        XaxisLayer{i}=Xax;
        YaxisLayer{i}=Yax;

    end
    LayerAvg{objn}.DensMLayer60=DensMLayer60;
    LayerAvg{objn}.PeakMLayer60=PeakMLayer60; 
     LayerAvg{objn}.ChgMLayer60=ChgMLayer60;
    LayerAvg{objn}.LatMlayer60=LatMlayer60; 
     LayerAvg{objn}.DensMLayer10=DensMLayer10;
    LayerAvg{objn}.PeakMLayer10= PeakMLayer10; 
     LayerAvg{objn}.ChgMLayer10=ChgMLayer10;
    LayerAvg{objn}.LatMlayer10=LatMlayer10;
      LayerAvg{objn}.XaxisLayer=XaxisLayer;
    LayerAvg{objn}.YaxisLayer=YaxisLayer;
    
end






    
    
    
