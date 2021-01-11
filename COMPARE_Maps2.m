function [ output_args ] = COMPARE_Maps(Control_path, Object_path,savepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%type=1; %use ranksum
folderMeanpath{1}=Control_path;
folderMeanpath{2}=Object_path;

close all;



 N=[];
for cmpind=1:2
    a=dir(fullfile(folderMeanpath{cmpind},'*.mat'));
    Cellfile=size(a);
    
    DirArea=[];
    Ag=[];
    Bnd=[];
    n=0;
    Nc=(Cellfile(1,1))*5;
    D60chg=[];
D10Pk=[];
D10chg=[];
D60Pk=[];
    Xaxis=cell(Nc,1);
    Yaxis=cell(Nc,1);
    Zpk60avg=cell(Nc,1);
    Zchg60avg=cell(Nc,1);
    Zden60avg=cell(Nc,1);
    Zpk10avg=cell(Nc,1);
    Zchg10avg=cell(Nc,1);
    Zden10avg=cell(Nc,1);
    Zlat10avg=cell(Nc,1);
    Zlat60avg=cell(Nc,1);
    
    ZchgD60avg=cell(Nc,1);
    ZdenD60avg=cell(Nc,1);
    
    ZlatD60avg=cell(Nc,1);
    for i= 1:1:Cellfile
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
                        
                        
                         D60chg=[D60chg;Cellrecord{fn}.DistWidthChargeTp1Hp60];
                         D60Pk=[D60Pk;Cellrecord{fn}.DistWidthPeakTp1Hp60];
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
                         
                    D60chg=[D60chg; Nanmaxtric];
                    D60Pk=[D60Pk; Nanmaxtric];
                        LdistPk60{n}=NaN;
                        LdistChg60{n}=NaN;
                    end
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
                        D10chg=[D10chg;Cellrecord{fn}.DistWidthChargeTp1Hp10];
                    D10Pk=[D10Pk;Cellrecord{fn}.DistWidthPeakTp1Hp10];
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
                        
                        
                    else Nanmaxtric= ones(1,NL).*NaN;
                        D10chg=[D10chg; Nanmaxtric];
                       D10Pk=[D10Pk; Nanmaxtric];
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
    Cs{cmpind}.Ag=Ag;
    Cs{cmpind}.LdistPk60=LdistPk60;
    Cs{cmpind}.LdistPk10=LdistPk10;
    Cs{cmpind}.LdistChg60=LdistChg60;
    Cs{cmpind}.LdistChg10=LdistChg10;
    
     Cs{cmpind}.D60chgavg=D60chg;
                   Cs{cmpind}.D10chgavg=D10chg;
                    Cs{cmpind}.D60pkavg=D60Pk;
                    Cs{cmpind}.D10pkavg=D10Pk;
    N=[N;n]
    
    
end

Dchg60=[];
Dchg10=[];
for cmpind=1:2
    Dchg60=[Dchg60;Cs{cmpind}.D60chgavg];
    Dchg10=[Dchg10;Cs{cmpind}.D10chgavg];
end
clusterN=2;
X=[Dchg60(:,1) Dchg60(:,3)]; CL=['k','b','r','m']
err=[];
Idcell=cell(1,100);
for ii=1:50
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};


find_figure('Scatter_group')
clf;
scatterhist(X(:,1),X(:,2),'Group',Id,'Location','SouthEast',...
    'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
    'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);



for cmpind=1:2
    D10y=Cs{cmpind}.D10chgavg(:,3);
     D10x=Cs{cmpind}.D10chgavg(:,1);
      D60y=Cs{cmpind}.D60chgavg(:,3);
       D60x=Cs{cmpind}.D60chgavg(:,1);
    sample=[Cs{cmpind}.D60chgavg(:,1) Cs{cmpind}.D60chgavg(:,3)];
    training=X;
    
    ids=classify(sample,training,Id);
    
    
    
    mC=[];vC=[];mBnd=[];
for ic=1:clusterN
    Id1=find(ids==ic);
mC=[mC;nanmean(D60y(Id1)) nanmean(D60y(Id1)) ];
vC=[vC;nanstd(D60y(Id1)) nanstd(D60y(Id1))];
mBnd=[mBnd;mean(Bnd(Id1,:),1)];
end

hold on
xa=[0.78 1.22];
for ic=2:clusterN
xa=[xa; xa(ic-1,:)+1];
end
find_figure(sprintf('Bargragh_ColumnarWidth%i',cmpind))
clf
errorbar(xa,mC,vC,'k*')
hold on
bar(mC)
hold off

find_figure(sprintf('Scatter_group%i',cmpind))
clf;
scatterhist(D60y,D60x,'Group',ids,'Location','SouthEast',...
    'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
    'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);

minage=0;
maxage=1000;
Ntotal=length(find(Ag>=minage&Ag<=maxage));
density=0.1;
Ncellpergroup=[];
savepath=sprintf('%s/group%i',folderMeanpath{cmpind},clusterN);
if ~isdir(savepath)
 mkdir(savepath);   
end

Ag=Cs{cmpind}.Ag;

Xaxis=Cs{cmpind}.Xaxis;
Yaxis=Cs{cmpind}.Yaxis;
Zpk60avg=Cs{cmpind}.Zpk60avg; 
Zchg60avg=Cs{cmpind}.Zchg60avg;
Zden60avg=Cs{cmpind}.Zden60avg;
Zlat60avg=Cs{cmpind}.Zlat60avg;
Zpk10avg=Cs{cmpind}.Zpk10avg; 
Zchg10avg=Cs{cmpind}.Zchg10avg;
Zden10avg=Cs{cmpind}.Zden10avg;
Zlat10avg=Cs{cmpind}.Zlat10avg;

ZchgD60avg=Cs{cmpind}.ZchgD60avg;
ZdenD60avg=Cs{cmpind}.ZdenD60avg;
ZlatD60avg=Cs{cmpind}.ZlatD60avg;

for c=1:clusterN
    
    Idc=find(Ag>=minage&Ag<=maxage&ids==c);
    N=length(find(Ag>=minage&Ag<=maxage&ids==c));
    Ncell=length(Idc);
    Ncellpergroup=[Ncellpergroup N./Ntotal]
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
    

    for i=1:length(Idc)
        Xavg=[Xavg Xaxis{Idc(i)}];
        Yavg=[Yavg Yaxis{Idc(i)}];
        Zpk60=[Zpk60 Zpk60avg{Idc(i)}];
        Zchg60=[Zchg60 Zchg60avg{Idc(i)}];
        Zden60=[Zden60 Zden60avg{Idc(i)}];
        Zpk10=[Zpk10 Zpk10avg{Idc(i)}];
        Zchg10=[Zchg10 Zchg10avg{Idc(i)}];
        Zden10=[Zden10 Zden10avg{Idc(i)}];
        Zlat10=[Zlat10 Zlat10avg{Idc(i)}];
        Zlat60=[Zlat60 Zlat60avg{Idc(i)}];
        
        ZchgD60=[ZchgD60 ZchgD60avg{Idc(i)}];
        ZdenD60=[ZdenD60 ZdenD60avg{Idc(i)}];
        
        ZlatD60=[ZlatD60 ZlatD60avg{Idc(i)}];

    end
    
    
    
   dx=30;
    dy=30;
xlim = [-420-0.1 960];
    ylim = [-480-0.1 480+0.1];
    
    xtick = xlim(1):dx:xlim(2);
    ytick = ylim(1):dy:ylim(2);
    xi = xtick(1:end-1)+dx/2;
    yi = ytick(1:end-1)+dy/2;
    xii=xi(1):1:xi(end);
    yii=yi(1):1:yi(end);
    
    
    
    % get rid of the places with few events
    [cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,Ncell,dens);
    
    
    sel = (cnt_exc_dir<=dens);
    isel = find(sel);
    selrm = [];
    for i=1:length(isel)
        sel2 = find(ind==isel(i)); selrm = [selrm sel2];
    end
    Xavg(selrm) = [];
    Zpk60(selrm) = [];
%     Zpksi(selrm) = [];
    Zchg10(selrm) =[];
    Zden60(selrm) = [];
%     Zdensi(selrm) = [];
    Yavg(selrm) = [];
    Zpk10(selrm)= [];
    Zchg60(selrm) = [];
%     Zchgsi(selrm) = [];
    Zden10(selrm) = [];
     Zlat10(selrm)= [];
    Zlat60(selrm) = [];
    ZchgD60(selrm) = [];
    ZdenD60(selrm) = [];
    ZlatD60(selrm) = [];
    [yi1,xi1]=meshgrid(yi,xi);
    [yii,xii]=meshgrid(yii,xii);
    avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
    Ind2=find(avg<=density); %avg(Ind2)=0;
    
    currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
    clf;
    avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avgdensity60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    
    
    %-------------------------------------------------direct density60
    
    avg = gridavg2(Xavg,Yavg,ZdenD60,xtick,ytick,false,Ncell,0);
    Ind2D=find(avg<=density); %avg(Ind2)=0;
    
    currfig= find_figure(sprintf('AvgDirectdensity60,Id==%i',c));
    clf;
    avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
    
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('AvgDirectdensity60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    %------------------------------------------------------------
    
    avg = gridavg2(Xavg,Yavg,Zpk60,xtick,ytick,false,Ncell,0);
    
    avg(Ind2)=0;
    currfig= find_figure(sprintf('AvgPk60,Id==%i',c));
    clf; avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    colorbar
    
    
    title(sprintf('#cell %i',Ncell))
    hold on
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avgpeak60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');



%-----------------------------------------------------------
    avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
    
    avg(Ind2)=0;
    currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
    clf; avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avgchg60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
      %__________________________________________________________________________________
     avg = gridavg2(Xavg,Yavg,ZchgD60,xtick,ytick,false,Ncell,0);
    
    avg(Ind2D)=0;
    currfig= find_figure(sprintf('AvgDirectchg60,Id==%i',c));
    clf; avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
    
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('AvgDirectchg60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
      %__________________________________________________________________________________ 
      
      
     avg =gridavg2(Xavg,Yavg,Zlat60,xtick,ytick,false,Ncell,0);
   avg(Ind2)=0;
    currfig= find_figure(sprintf('Avglatency60_2,Id==%i',c));
    clf;
    avg=avg';
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avglatency60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    
    %__________________________________________________________________________________
      
     avg =gridavg2(Xavg,Yavg,ZlatD60,xtick,ytick,false,Ncell,0);
   avg(Ind2D)=0;
    currfig= find_figure(sprintf('AvgDirectlatency60_2,Id==%i',c));
    clf;
    avg=avg';
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
    
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('AvgDirectlatency60Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    
    %__________________________________________________________________________________
    
    avg =gridavg2(Xavg,Yavg,Zden10,xtick,ytick,false,Ncell,0);
     Ind10=find(avg<=density);
%     avg(Ind10)=0;
    currfig= find_figure(sprintf('Avgdensity10_2,Id==%i',c));
    clf;
    avg=avg';
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
    
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avgdensity10Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
   %-----------------------------------------------------------------
   
    avg = gridavg2(Xavg,Yavg,Zpk10,xtick,ytick,false,Ncell,0); avg(Ind10)=0;
    currfig= find_figure(sprintf('AvgPk10,Id==%i',c));
    clf;
    avg=avg';
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.51560693641618 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avgpeak10Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    avg =gridavg2(Xavg,Yavg,Zchg10,xtick,ytick,false,Ncell,0);
    avg(Ind10)=0;
    currfig= find_figure(sprintf('Avgchg10,Id==%i',c));
    clf;
    avg=avg';
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
    
     for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avgchg10Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    
    
    %__________________________________________________________________________________
     avg =gridavg2(Xavg,Yavg,Zlat10,xtick,ytick,false,Ncell,0);
     avg(Ind10)=0;
    currfig= find_figure(sprintf('Avglatency10_2,Id==%i',c));
    clf;
    avg=avg';
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    hold on
    
    for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)],'LineWidth',2,'Color',[1 1 1])
    end
    
    plot(0,0,'wo','MarkerSize',15)
    
    hold(axes1,'all');
    xmat=sum(avg);
    ymat=sum(avg,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yi,xmat,'k-','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
    plot(xi,ymat,'k-','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('Avglatency10Id%iAge%ito%i.fig',c,minage,maxage)),'fig');
    

    
     avg60 = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
    Ind2=find(avg60<=density); avg60(Ind2)=0;
    
    avg10 = gridavg2(Xavg,Yavg,Zden10,xtick,ytick,false,Ncell,0);
    Ind2=find(avg10<=density); avg10(Ind2)=0;
    
    currfig= find_figure(sprintf('CompareDensity,Id==%i',c));
    clf;
    avg60=avg60';
    avg10=avg10';
   % C2=1-(~(avg60.*avg10));
   
    avg2_60=interp2(yi1,xi1,avg60,yii,xii,'cubic');
      avg2_10=interp2(yi1,xi1,avg10,yii,xii,'cubic');
%       C22=interp2(yi1,xi1,C2,yii,xii,'cubic');
      C=repmat(avg2_60,[1,1,3]);
      C1=repmat(avg2_60,[1,1,3]);
      C2=repmat(avg2_60,[1,1,3]);
      C2(:,:,1)=avg2_60;
      %[cid1,cid2]=find(avg2_60+avg2_10);
      C2(:,:,2)=avg2_60;
      %C(cid1,cid2,2)=0;
      C2(:,:,3)=0;
      C2(find(C2<0|isnan(C2)|C2>1))=0;
      
      C1(:,:,1)=0;
      %[cid1,cid2]=find(avg2_60+avg2_10);
      C1(:,:,2)=avg2_10;
      %C(cid1,cid2,2)=0;
      C1(:,:,3)=avg2_10;
      C1(find(C1<0|isnan(C1)|C1>1))=0;
      C=(1-C1).*(1-C2);
      bd=0.15;
%       for i=1:length(C(:,1,1))
%           for j=1:length(C(1,:,1))
%               
%                   if C(i,j,1)<bd&&C(i,j,2)<bd&&C(i,j,3)<bd
%                       C(i,j,1)=1;
%                       C(i,j,2)=1;
%                       C(i,j,3)=1;
%                   end
%           end
%       end
      
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),C); axis image;axis ij;
    
    
    title(sprintf('#cell %i',Ncell))
    hold on
    for ib=1:length(mBnd(1,:))
    line([-100 100],[mBnd(c,ib) mBnd(c,ib)])
    end
    plot(0,0,'wo','MarkerSize',15)
    
  % hold(axes1,'all');
    xmat60=sum(avg2_60);
    ymat60=sum(avg2_60,2);
    xmat10=sum(avg2_10);
    ymat10=sum(avg2_10,2);
    axes2 = axes('Parent',currfig,'XAxisLocation','top',...
        'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
    plot(yii(1,:),xmat60./max(xmat60),'r-','LineWidth',2)
    hold on
    plot(yii(1,:),xmat10./max(xmat10),'b','LineWidth',2)
    set(gca,'XLim',[min(ytick) max(ytick)])
    hold off
    box off
    axes3 = axes('Parent',currfig,'YDir','reverse',...
        'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
     plot(xii(:,1),ymat60./max(ymat60),'r-','LineWidth',2)
    hold on
    plot(xii(:,1),ymat10./max(ymat10),'b','LineWidth',2)
    set(gca,'XLim',[min(xtick) max(xtick)])
    
    view(90,90)
    saveas(currfig,fullfile(savepath,sprintf('CompareDensityId%iAge%ito%i.fig',c,minage,maxage)),'fig');
end








end




















 dx=30;
    dy=30;
xlim = [-420-0.1 960];
    ylim = [-480-0.1 480+0.1];
    
    xtick = xlim(1):dx:xlim(2);
    ytick = ylim(1):dy:ylim(2);
    xi = xtick(1:end-1)+dx/2;
    yi = ytick(1:end-1)+dy/2;
    dxint=2;
    xii=xi(1):dxint:xi(end);
    yii=yi(1):dxint:yi(end);
  [yi1,xi1]=meshgrid(yi,xi);
        [yii,xii]=meshgrid(yii,xii);
for objn=1:2
    Ncell=N(objn);
    
    
    dens=min(10,0.8*Ncell);
   
   
   
    
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
    
    
    
    
     DensMLayer60=cell(Ncell,1);
     ChgMLayer60=cell(Ncell,1);
     PeakMLayer60=cell(Ncell,1);
    LatMlayer60=cell(Ncell,1); 
    DensMLayer10=cell(Ncell,1);
     ChgMLayer10=cell(Ncell,1);
     PeakMLayer10=cell(Ncell,1);
    LatMlayer10=cell(Ncell,1); 
    XaxisLayer=cell(Ncell,1);
       YaxisLayer=cell(Ncell,1);
  k=1;
    for i=1:Ncell

        
        
        
        
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zden60avg{i},xtick,ytick,false,Ncell,0);
        
        
        

        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
       
%                 currfig= find_figure(sprintf('Singlecelldensity60_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
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
%         colorbar
        %------------------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zpk60avg{i},xtick,ytick,false,Ncell,0);
        
      
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
       
%           currfig= find_figure(sprintf('Singlecellpeak60_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij; axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            PkM60{ib}=avg2(ind{ib},:);
            
        end
        
        %---------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zchg60avg{i},xtick,ytick,false,Ncell,0);
        

        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
       
%                currfig= find_figure(sprintf('SinglecellChg60_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1) 
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            ChgM60{ib}=avg2(ind{ib},:);
            
        end
        %---------------------------------------------------------------
        avg = gridavg2(Xaxis{i},Yaxis{i},Zlat60avg{i},xtick,ytick,false,Ncell,0);
        

        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
       
%                 currfig= find_figure(sprintf('SinglecellLatency60_%d_%d',objn,k));
%        subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            LatM60{ib}=avg2(ind{ib},:);
            
        end
        %__________________________________________________________________
        
        
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zden10avg{i},xtick,ytick,false,Ncell,0);
        
      
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
        
%           currfig= find_figure(sprintf('Singlecelldensity10_%d_%d',objn,k));
%       subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        hold on
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        
        for ib=1:length(Bnd(1,:))-1
            
            densM10{ib}=avg2(ind{ib},:);
            
            
        end
        
        
        
        
        %     colormap(winter)
%         colorbar
        %------------------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zpk10avg{i},xtick,ytick,false,Ncell,0);
        
       
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
        
%          currfig= find_figure(sprintf('Singlecellpeak10_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            PkM10{ib}=avg2(ind{ib},:);
            
        end
        
        %---------------------------------------------------------------
        
        avg = gridavg2(Xaxis{i},Yaxis{i},Zchg10avg{i},xtick,ytick,false,Ncell,0);
        
       
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
       
%          currfig= find_figure(sprintf('SinglecellChg10_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            ChgM10{ib}=avg2(ind{ib},:);
            
        end
        %---------------------------------------------------------------
        avg = gridavg2(Xaxis{i},Yaxis{i},Zlat10avg{i},xtick,ytick,false,Ncell,0);
        
        
        avg=avg';
        
        avg2=interp2(yi1,xi1,avg,yii,xii,'linear');
     
%         currfig= find_figure(sprintf('SinglecellLatency10_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
        end
        %     colormap(winter)
%         colorbar
        
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
        if mod(i,9)==0
            k=k+1;
        end

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

objn=2;
Bnd=Cs{objn}.Bnd;
Ncell=N(objn);
density=0.1; 
%avg Layer1
for j=1:length(Bnd(1,:))-1
   Zdens60=[];
Zpk60=[];
Zchg60=[];
Zdens10=[];
Zpk10=[];
Zchg10=[];
X=[];
Y=[]; 
xtick=LayerAvg{objn}.YaxisLayer{1}{1}(1,:)-dxint/2;
ytick=-1*mean(abs(Bnd(:,j)-Bnd(:,j+1)))/2-dxint/2:dxint:mean(abs(Bnd(:,j)-Bnd(:,j+1)))/2-dxint/2;
xi=xtick+0.5;
yi=ytick+0.5;
for i=1:Ncell
  X=[X;LayerAvg{objn}.XaxisLayer{i}{j}(:)];
  Y=[Y;LayerAvg{objn}.YaxisLayer{i}{j}(:)];
  Zdens60=[Zdens60;LayerAvg{objn}.DensMLayer60{i}{j}(:)];
  Zchg60=[Zchg60;LayerAvg{objn}.ChgMLayer60{i}{j}(:)];
  Zpk60=[Zpk60;LayerAvg{objn}.PeakMLayer60{i}{j}(:)];
  Zdens10=[Zdens10;LayerAvg{objn}.DensMLayer10{i}{j}(:)];
  Zchg10=[Zchg10;LayerAvg{objn}.ChgMLayer10{i}{j}(:)];
  Zpk10=[Zpk10;LayerAvg{objn}.PeakMLayer10{i}{j}(:)];
end

avg = gridavg2(Y,X,Zdens60,xtick,ytick,false,Ncell,0);
Ind60=find(avg<density);

currfig= find_figure(sprintf('AvgDensity60_object%b',objn));
subplot(length(Bnd(1,:))-1,1,j)
imagesc(xi,yi,avg); axis image;axis ij;
currfig= find_figure(sprintf('Avgchg60_object%b',objn));
avg = gridavg2(Y,X,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind60)=0;
subplot(length(Bnd(1,:))-1,1,j)
imagesc(xi,yi,avg); axis image;axis ij;

currfig= find_figure(sprintf('Avgpk60_object%b',objn));
avg = gridavg2(Y,X,Zpk60,xtick,ytick,false,Ncell,0);
avg(Ind60)=0;
subplot(length(Bnd(1,:))-1,1,j)
imagesc(xi,yi,avg); axis image;axis ij;


avg = gridavg2(Y,X,Zdens10,xtick,ytick,false,Ncell,0);
Ind10=find(avg<density);
currfig= find_figure(sprintf('AvgDensity10_object%b',objn));
subplot(length(Bnd(1,:))-1,1,j)
imagesc(xi,yi,avg); axis image;axis ij;
currfig= find_figure(sprintf('Avgchg10_object%b',objn));
avg = gridavg2(Y,X,Zchg10,xtick,ytick,false,Ncell,0);
avg(Ind10)=0;
subplot(length(Bnd(1,:))-1,1,j)
imagesc(xi,yi,avg); axis image;axis ij;

currfig= find_figure(sprintf('Avgpk10_object%b',objn));
avg = gridavg2(Y,X,Zpk10,xtick,ytick,false,Ncell,0);
avg(Ind10)=0;
subplot(length(Bnd(1,:))-1,1,j)
imagesc(xi,yi,avg); axis image;axis ij;

end




    
    
    
