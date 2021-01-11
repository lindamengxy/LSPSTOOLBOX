folderMeanpath='/Users/lindameng/Documents/study/Myelin/wholecell/FVBNJ/meanValueeachcell_eventwindow50_direct8';
a=dir(fullfile(folderMeanpath,'*.mat'));
Cellfile=size(a);
Pchglayer60=[];
Ppeaklayer60=[];
Pchglayer10=[];
Ppeaklayer10=[]; 
Pchglayersi=[];
Ppeaklayersi=[];
chglayer60=[];
peaklayer60=[];
chglayer10=[];
peaklayer10=[];
chglayersi=[];
peaklayersi=[];
DistWidthPeak=[];
DistWidthCharge=[];
SkPeaklayer=[];
Skchglayer=[];
Skdenslayer=[];
DistWidthdens=[];
CorrHighMg60=[];
CorrHighMg10=[];
CorrPtx60=[];
AreaSi=[];
CorrPtx50=[];
arealayer10=[];
arealayersi=[];
arealayer60=[];
Nsi=0;
Pchgsilayer=[];
Ppeaksilayer=[];
DistWidthPeaksi=[];
DistWidthChargesi=[];
DirArea=[];
AR=[];
D60Pk=[];
Bnd=[];
D60chg=[];
D10Pk=[];
D10chg=[];
DsiPk=[];
Dsichg=[];
DistWidthdenssi=[];

Ag=[];
k60=0;k10=0;
for i= 1:1:Cellfile
    Cg= load(fullfile(folderMeanpath,a(i).name));
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
                Nsi=Nsi+1;
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Xaxis{Nsi}=X;
                Yaxis{Nsi}=Y;
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
%                 Zpksi=zeros(size(X)).*NaN;
%                 Zchgsi=zeros(size(X)).*NaN;
%                 Zdensi=zeros(size(X));
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
                    k60=k60+1;
                    D60chg=[D60chg;Cellrecord{fn}.DistWidthChargeTp1Hp60];
                    D60Pk=[D60Pk;Cellrecord{fn}.DistWidthPeakTp1Hp60];
                    Pchglayer60=[Pchglayer60;Cellrecord{fn}.PchgTp1Hp60];
                    Ppeaklayer60=[Ppeaklayer60;Cellrecord{fn}.PpeakTp1Hp60];
                    chglayer60=[chglayer60;Cellrecord{fn}.chgTp1Hp60];
                    peaklayer60=[peaklayer60;Cellrecord{fn}.peakTp1Hp60];
                    arealayer60=[arealayer60;Cellrecord{fn}.AreaTp1Hp60];
                    
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
                    
                    
                    StimcorrX60{k60}=X;
                    StimcorrY60{k60}=Y;
                    ChgCell60{k60}=Zchg60;
                   CorrHighMg60=[CorrHighMg60,Cellrecord{fn}.CorrflagHighMg60];
                    
                else Nanmaxtric= ones(1,NL).*NaN;
                    D60chg=[D60chg; Nanmaxtric];
                    D60Pk=[D60Pk; Nanmaxtric];
                    Nanmaxtric2=ones(1,5).*NaN;
                    Pchglayer60=[Pchglayer60;Nanmaxtric2];
                    Ppeaklayer60=[Ppeaklayer60;Nanmaxtric2];
                    chglayer60=[chglayer60;Nanmaxtric2];
                    peaklayer60=[peaklayer60;Nanmaxtric2];
                    arealayer60=[arealayer60;Nanmaxtric2];
                     CorrHighMg60=[CorrHighMg60,NaN];
                end
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
                    k10=k10+1;
                    D10chg=[D10chg;Cellrecord{fn}.DistWidthChargeTp1Hp10];
                    D10Pk=[D10Pk;Cellrecord{fn}.DistWidthPeakTp1Hp10];
                    Pchglayer10=[Pchglayer10;Cellrecord{fn}.PchgTp1Hp10];
                    Ppeaklayer10=[Ppeaklayer10;Cellrecord{fn}.PpeakTp1Hp10];
                     chglayer10=[chglayer10;Cellrecord{fn}.chgTp1Hp10];
                    peaklayer10=[peaklayer10;Cellrecord{fn}.peakTp1Hp10];
                    arealayer10=[arealayer10;Cellrecord{fn}.AreaTp1Hp10];
                    Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg10));
                    Zpk10(Indpk)=Cellrecord{fn}.meanpeakHighMg10(Indpk);
                    Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg10));
                    Zchg10(Indchg)=Cellrecord{fn}.meanareaHighMg10(Indchg);
                     Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg10));
                    Zlat10(Indlat)=Cellrecord{fn}.meanlatencyHighMg10(Indlat);
                    Indden=find((Cellrecord{fn}.meanflagHighMg10>0));
                    Zden10(Indden)=Cellrecord{fn}.meanflagHighMg10(Indden);
                     CorrHighMg10=[CorrHighMg10,Cellrecord{fn}.CorrflagHighMg10];
                
                    StimcorrX10{k10}=X;
                    StimcorrY10{k10}=Y;
                    ChgCell10{k10}=Zchg10;
                    
                else Nanmaxtric= ones(1,NL).*NaN;
                    D10chg=[D10chg; Nanmaxtric];
                    D10Pk=[D10Pk; Nanmaxtric];
                    
                    Nanmaxtric2=ones(1,6).*NaN;
                    Pchglayer10=[Pchglayer10;Nanmaxtric2];
                    Ppeaklayer10=[Ppeaklayer10;Nanmaxtric2];
                     arealayer10=[arealayer10;Nanmaxtric2];
                        chglayer10=[chglayer10;Nanmaxtric2];
                    peaklayer10=[peaklayer10;Nanmaxtric2];
                    CorrHighMg10=[CorrHighMg10,NaN];
                end
%                 if isfield(Cellrecord{fn},'DistWidthPeakTp2Hp50si')
%                     Dsichg=[Dsichg;Cellrecord{fn}.DistWidthChargeTp2Hp50si];
%                     DsiPk=[DsiPk;Cellrecord{fn}.DistWidthPeakTp2Hp50si];
%                     Pchglayersi=[Pchglayersi;Cellrecord{fn}.PchgTp2si];
%                     Ppeaklayersi=[Ppeaklayersi;Cellrecord{fn}.PpeakTp2si];
%                     chglayersi=[chglayersi;Cellrecord{fn}.chgTp2si];
%                     peaklayersi=[peaklayersi;Cellrecord{fn}.peakTp2si];
%                     arealayersi=[arealayersi;Cellrecord{fn}.AreaTp2si];
%                     Indsi=find(Cellrecord{fn}.Matrix_NMDAonly>0);
%                       AreaSi=[AreaSi;length(Indsi)];
%                     Zpksi(Indsi)=Cellrecord{fn}.meanpeakPtxInd50(Indsi);
%                     
%                     Zchgsi(Indsi)=Cellrecord{fn}.meanareaPtxInd50(Indsi);
%                     
%                     
%                     Zdensi(Indsi)=Cellrecord{fn}.Matrix_NMDAonly(Indsi);
%                 else Nanmaxtric= ones(1,NL).*NaN;
%                     Dsichg=[Dsichg; Nanmaxtric];
%                     DsiPk=[DsiPk; Nanmaxtric];
%                     AreaSi=[AreaSi;length(Indsi)];
%                     Nanmaxtric2=ones(1,4).*NaN;
%                     Pchglayersi=[Pchglayersi;Nanmaxtric2];
%                     Ppeaklayersi=[Ppeaklayersi;Nanmaxtric2];
%                     AreaSi=[AreaSi;NaN];
%                     chglayersi=[chglayersi;Nanmaxtric2];
%                     peaklayersi=[peaklayersi;Nanmaxtric2];
%                     arealayersi=[arealayersi;Nanmaxtric2];
%                 end
                Zpk60avg{Nsi}=Zpk60;
                Zchg60avg{Nsi}=Zchg60;
                Zden60avg{Nsi}=Zden60;
                Zpk10avg{Nsi}=Zpk10;
                Zchg10avg{Nsi}=Zchg10;
                Zden10avg{Nsi}=Zden10;
                Zlat10avg{Nsi}=Zlat10;
                Zlat60avg{Nsi}=Zlat60;
                  
                ZchgD60avg{Nsi}=ZchgD60;
                ZdenD60avg{Nsi}=ZdenD60;
              
                ZlatD60avg{Nsi}=ZlatD60;
%                 Zpksiavg{Nsi}=Zpksi;
%                 Zchgsiavg{Nsi}=Zchgsi;
%                 Zdensiavg{Nsi}=Zdensi;
            end
            
            
            
        end
    end
end



% PCA

D60x=D60Pk(:,1);
D60xext=D60Pk(:,2);
D60y=D60Pk(:,3);
D60yratio=D60Pk(:,4);
D10x=D10Pk(:,1);
D10xext=D10Pk(:,2);
D10yratio=D10Pk(:,4);
D10y=D10Pk(:,3);
% Dsix=DsiPk(:,1);
% Dsixext=DsiPk(:,4);
% Dsiy=DsiPk(:,5);
% Dsiyratio=DsiPk(:,6);
Ppk60=Ppeaklayer60(:,4);
Pchg60=Pchglayer60(:,4);
Ppk10=Ppeaklayer10(:,4);
Pchg10=Pchglayer10(:,4);
% Ppksi=Ppeaklayersi(:,4);
% Pchgsi=Pchglayersi(:,4);
pk60=peaklayer60(:,4);
chg60=chglayer60(:,4);
pk10=peaklayer10(:,4);
chg10=chglayer10(:,4);
% pksi=peaklayersi(:,4);
% chgsi=chglayersi(:,4);
% areasi=arealayersi(:,4);
area60=arealayer60(:,4);
area10=arealayer10(:,4);
clusterN=2;




%     iii=find(~isnan(Pchg60)&~isnan(Pchg10)&~isnan(Pchgsi));
% Pchg60=Pchg60(iii);
% Pchg10=Pchg10(iii);
% Pchgsi=Pchgsi(iii);
% xext=xext(iii);
% yy=yy(iii);
X=[D60xext D60yratio]; CL=['k','b','r','m']
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};


mC=[];vC=[];mBnd=[];
for ic=1:clusterN
    Id1=find(Id==ic)
mC=[mC;nanmean(D60y(Id1)) nanmean(D10y(Id1)) ];
vC=[vC;nanstd(D60y(Id1)) nanstd(D60y(Id1))];
mBnd=[mBnd;mean(Bnd(Id1,:),1)];
end

hold on
xa=[0.78 1.22];
for ic=2:clusterN
xa=[xa; xa(ic-1,:)+1];
end
find_figure('Bargragh_ColumnarWidth')
clf
errorbar(xa,mC,vC,'k*')
hold on
bar(mC)
hold off
xhist=-2500:1000:60000;

% Id1=find(Id==1);
% Id2=find(Id==2);
% Id3=find(Id==3);
% Id4=find(Id==4);
hist60A=cell(clusterN);
hist10A=cell(clusterN);
histsiA=cell(clusterN);
find_figure('hist_l56inputarea')
clf
for ic=1:clusterN
    Id1=find(Id==ic)
subplot(clusterN+1,2,1+(ic-1)*2)
[y1,x1]=hist(area60(Id1),xhist);
hist60A{ic}=area60(Id1);
y1=y1./sum(y1);
bar(x1,y1,CL(ic),'EdgeColor','none')
box off
subplot(clusterN+1,2,1+(clusterN)*2)
hold on
h=cdfplot(area60(Id1))

set(h,'Color',CL(ic))
 box off
 
subplot(clusterN+1,2,2+(ic-1)*2)
[y1,x1]=hist(area10(Id1),xhist)
hist10A{ic}=area10(Id1);
y1=y1./sum(y1);
bar(x1,y1,CL(ic),'EdgeColor','none')

box off
subplot(clusterN+1,2,2+(clusterN)*2)
hold on
h=cdfplot(area10(Id1))

set(h,'Color',CL(ic))
 box off

% subplot(clusterN+1,3,3+(ic-1)*3)
% [y,x]=hist(areasi(Id1),xhist)
% histsiA{ic}=areasi(Id1);
% y=y./sum(y);
% bar(x,y,CL(ic),'EdgeColor','none')
% 
% box off
% subplot(clusterN+1,3,3+(clusterN)*3)
% hold on
% h=cdfplot(areasi(Id1))
% 
% set(h,'Color',CL(ic))
%  box off
end
hist60=cell(clusterN);
hist10=cell(clusterN);
histsi=cell(clusterN);
find_figure('hist_l56input')
clf
xhist1=0:40:1000;
xhist2=0:10:400;
for ic=1:clusterN
    Id1=find(Id==ic) 
subplot(clusterN+1,2,1+(ic-1)*2)
[y1,x1]=hist(chg60(Id1),xhist2);
hist60{ic}=chg60(Id1);
y1=y1./sum(y1);
bar(x1,y1,CL(ic),'EdgeColor','none')

box off
subplot(clusterN+1,2,1+(clusterN)*2)
hold on
h=cdfplot(chg60(Id1))

set(h,'Color',CL(ic))
 box off
 
subplot(clusterN+1,2,2+(ic-1)*2)
[y1,x1]=hist(chg10(Id1),xhist2)
hist10{ic}=chg10(Id1);
y1=y1./sum(y1);
bar(x1,y1,CL(ic),'EdgeColor','none')
 box off
subplot(clusterN+1,2,2+(clusterN)*2)
hold on
h=cdfplot(chg10(Id1))

set(h,'Color',CL(ic))
 box off
% subplot(clusterN+1,2,3+(ic-1)*2)
% [y,x]=hist(chgsi(Id1),xhist2)
% histsi{ic}=chgsi(Id1);
% y=y./sum(y);
% bar(x,y,CL(ic),'EdgeColor','none')
% 
%  box off
% subplot(clusterN+1,2,3+(clusterN)*2)
% hold on
% h=cdfplot(chgsi(Id1))
% 
% set(h,'Color',CL(ic))
%   box off
end

ks60=ones(clusterN,clusterN).*NaN;
ks10=ones(clusterN,clusterN).*NaN;
kssi=ones(clusterN,clusterN).*NaN;
P60=ones(clusterN,clusterN).*NaN;
P10=ones(clusterN,clusterN).*NaN;
Psi=ones(clusterN,clusterN).*NaN;
for i=1:size(hist60,1)-1
    for j=i+1:size(hist60,1)
   [ks60(i,j),P60(i,j)]=kstest2(hist60{i},hist60{j});
   [ks10(i,j),P10(i,j)]=kstest2(hist10{i},hist10{j});
%    [kssi(i,j),Psi(i,j)]=kstest2(histsi{i},histsi{j});
    end
end
ks60
ks10
% kssi
P60
P10
% Psi

ks60A=ones(clusterN,clusterN).*NaN;
ks10A=ones(clusterN,clusterN).*NaN;
% kssiA=ones(clusterN,clusterN).*NaN;
P60A=ones(clusterN,clusterN).*NaN;
P10A=ones(clusterN,clusterN).*NaN;
% PsiA=ones(clusterN,clusterN).*NaN;


RS60A=ones(clusterN,clusterN).*NaN;
RS10A=ones(clusterN,clusterN).*NaN;
% RSsiA=ones(clusterN,clusterN).*NaN;
Pr60A=ones(clusterN,clusterN).*NaN;
Pr10A=ones(clusterN,clusterN).*NaN;
% PrsiA=ones(clusterN,clusterN).*NaN;


RS60=ones(clusterN,clusterN).*NaN;
RS10=ones(clusterN,clusterN).*NaN;
% RSsi=ones(clusterN,clusterN).*NaN;
Pr60=ones(clusterN,clusterN).*NaN;
Pr10=ones(clusterN,clusterN).*NaN;
% Prsi=ones(clusterN,clusterN).*NaN;
for i=1:size(hist60A,1)-1
    for j=i+1:size(hist60A,1)
   [ks60A(i,j),P60A(i,j)]=kstest2(hist60A{i},hist60A{j});
   [ks10A(i,j),P10A(i,j)]=kstest2(hist10A{i},hist10A{j});
%    [kssiA(i,j),PsiA(i,j)]=kstest2(histsiA{i},histsiA{j});
    end
end


for i=1:size(hist60A,1)-1
    for j=i+1:size(hist60A,1)
   [RS60A(i,j),Pr60A(i,j)]=ranksum(hist60A{i},hist60A{j});
   [RS10A(i,j),Pr10A(i,j)]=ranksum(hist10A{i},hist10A{j});
%    [RSsiA(i,j),PrsiA(i,j)]=ranksum(histsiA{i},histsiA{j});
    end
end


for i=1:size(hist60A,1)-1
    for j=i+1:size(hist60A,1)
   [RS60(i,j),Pr60(i,j)]=ranksum(hist60{i},hist60{j});
   [RS10(i,j),Pr10(i,j)]=ranksum(hist10{i},hist10{j});
%    [RSsi(i,j),Prsi(i,j)]=ranksum(histsi{i},histsi{j});
    end
end
% figure(2)
%  scatterhist(xext, Ppk60,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
% figure(3)
%  scatterhist(xext2, Ppk10,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
% figure(5)
find_figure('Scatter_group')
clf;
scatterhist(D60y,D60x,'Group',Id,'Location','SouthEast',...
    'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
    'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
% % Y=pdist(X);
% % Z=linkage(Y);
% % dendrogram(Z);


%---------------------------------------------------------------------------
% accumation of columnar width for two groups

% Cw602=[];Cw601=[];Cw101=[];Cw102=[];Cwsi2=[];Cwsi1=[];
%
% xhist=0:5:1000;
%  for i=1:length(xhist)
%          t=y(Id2)<=xhist(i);
%          Cw602=[Cw602 sum(t)];
% end
% Cw602=Cw602./length(Id2);
% Cw601=[];
%
%  for i=1:length(xhist)
%          t=y(Id1)<=xhist(i);
%          Cw601=[Cw601 sum(t)];
% end
% Cw601=Cw601./length(Id1);
%
% Cw502=[];
%
%  for i=1:length(xhist)
%          t=y2(Id2)<=xhist(i);
%          Cw502=[Cw502 sum(t)];
% end
% Cw502=Cw502./length(Id2);
%
%
% Cw501=[];
%
%  for i=1:length(xhist)
%          t=y2(Id1)<=xhist(i);
%          Cw501=[Cw501 sum(t)];
% end
% Cw501=Cw501./length(Id1);
%
% Cwsi2=[];
%
%  for i=1:length(xhist)
%          t=y3(Id2)<=xhist(i);
%          Cwsi2=[Cwsi2 sum(t)];
% end
% Cwsi2=Cwsi2./length(Id2);
% Cwsi1=[];
%
%  for i=1:length(xhist)
%          t=y3(Id1)<=xhist(i);
%          Cwsi1=[Cwsi1 sum(t)];
% end
% Cwsi1=Cwsi1./length(Id1);
%
% figure(7)
% plot(xhist,Cw601,xhist,Cw501,xhist,Cwsi1,xhist,Cw602,xhist,Cw502,xhist,Cwsi2)
%----------------------------------------------------------------------------------

%----------------------------------------------------------------------------------
%plot the average map based on the group

minage=0;
maxage=1000;
Ntotal=length(find(Ag>=minage&Ag<=maxage));
density=0.1;
Ncellpergroup=[];
savepath=sprintf('%s/group%i',folderMeanpath,clusterN);
if ~isdir(savepath)
 mkdir(savepath);   
end
for c=1:clusterN
    
    Id1=find(Ag>=minage&Ag<=maxage&Id==c);
    N=length(find(Ag>=minage&Ag<=maxage&Id==c));
    Ncell=length(Id1);
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
    
%     Zpksi=[];
%     Zchgsi=[];
%     Zdensi=[];
    for i=1:length(Id1)
        Xavg=[Xavg Xaxis{Id1(i)}];
        Yavg=[Yavg Yaxis{Id1(i)}];
        Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
        Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
        Zden60=[Zden60 Zden60avg{Id1(i)}];
        Zpk10=[Zpk10 Zpk10avg{Id1(i)}];
        Zchg10=[Zchg10 Zchg10avg{Id1(i)}];
        Zden10=[Zden10 Zden10avg{Id1(i)}];
        Zlat10=[Zlat10 Zlat10avg{Id1(i)}];
        Zlat60=[Zlat60 Zlat60avg{Id1(i)}];
        
       %  ZpkD60=[ZpkD60 ZpkD60avg{Id1(i)}];
        ZchgD60=[ZchgD60 ZchgD60avg{Id1(i)}];
        ZdenD60=[ZdenD60 ZdenD60avg{Id1(i)}];
        
        ZlatD60=[ZlatD60 ZlatD60avg{Id1(i)}];
%         Zpksi=[Zpksi Zpksiavg{Id1(i)}];
%         Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
%         Zdensi=[Zdensi Zdensiavg{Id1(i)}];
    end
    
    
    
    
    xlim = [-300-dx-0.1 1000];
    ylim = [-500-dy-0.1 500];
    
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
    
    %__________________________________________________________________________________
   
%     avg = gridavg2(Xavg,Yavg,Zdensi,xtick,ytick,false,Ncell,0);
%     Indsi=find(avg<=density); avg(Indsi)=0;
%     currfig= find_figure(sprintf('Avgdensitysi,Id==%i',c));
%     clf;
%     avg=avg';
%     avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
%     axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
%         'DataAspectRatio',[1 1 1]);
%     imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colormap(winter)
%     colorbar
%     
%     title(sprintf('#cell %i',Ncell))
%     hold on
%     
%     
%     
%     plot(0,0,'wo','MarkerSize',15)
%     
%     hold(axes1,'all');
%     xmat=sum(avg);
%     ymat=sum(avg,2);
%     axes2 = axes('Parent',currfig,'XAxisLocation','top',...
%         'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
%     plot(yi,xmat,'k-','LineWidth',2)
%     set(gca,'XLim',[min(ytick) max(ytick)])
%     box off
%     axes3 = axes('Parent',currfig,'YDir','reverse',...
%         'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
%     plot(xi,ymat,'k-','LineWidth',2)
%     set(gca,'XLim',[min(xtick) max(xtick)])
%     
%     
%     view(90,90)
%     saveas(currfig,fullfile(savepath,sprintf('AvgdensitysiId%iAge%ito%i.fig',c,minage,maxage)),'fig');
%     avg= gridavg2(Xavg,Yavg,Zpksi,xtick,ytick,false,Ncell,0);
%     avg(Indsi)=0;
%     currfig= find_figure(sprintf('AvgPksi,Id==%i',c));
%     clf; avg=avg';
%     avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
%     axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
%         'DataAspectRatio',[1 1 1]);
%     imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colorbar
%     
%     title(sprintf('#cell %i',Ncell))
%     hold on
%     
%     
%     
%     plot(0,0,'wo','MarkerSize',15)
%     
%     hold(axes1,'all');
%     xmat=sum(avg);
%     ymat=sum(avg,2);
%     axes2 = axes('Parent',currfig,'XAxisLocation','top',...
%         'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
%     plot(yi,xmat,'k-','LineWidth',2)
%     set(gca,'XLim',[min(ytick) max(ytick)])
%     box off
%     axes3 = axes('Parent',currfig,'YDir','reverse',...
%         'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
%     plot(xi,ymat,'k-','LineWidth',2)
%     set(gca,'XLim',[min(xtick) max(xtick)])
%     
%     
%     view(90,90)
%     saveas(currfig,fullfile(savepath,sprintf('AvgpeaksiId%iAge%ito%i.fig',c,minage,maxage)),'fig');
%     avg =gridavg2(Xavg,Yavg,Zchgsi,xtick,ytick,false,Ncell,0);
%     avg(Indsi)=0;
%     currfig= find_figure(sprintf('Avgchgsi,Id==%i',c));
%     clf;
%     avg=avg';
%     avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
%     axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
%         'DataAspectRatio',[1 1 1]);
%     imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
%     colorbar
%     
%     title(sprintf('#cell %i',Ncell))
%     hold on
%     
%     
%     
%     plot(0,0,'wo','MarkerSize',15)
%     
% %    hold(axes1,'all');
%     xmat=sum(avg);
%     ymat=sum(avg,2);
%     axes2 = axes('Parent',currfig,'XAxisLocation','top',...
%         'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
%     plot(yi,xmat,'k-','LineWidth',2)
%     set(gca,'XLim',[min(ytick) max(ytick)])
%     box off
%     axes3=axes('Parent',currfig,'YDir','reverse',...
%         'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
%     plot(xi,ymat,'k-','LineWidth',2)
%     set(gca,'XLim',[min(xtick) max(xtick)])
%     
%     
%     view(90,90)
%     saveas(currfig,fullfile(savepath,sprintf('AvgchgsiId%iAge%ito%i.fig',c,minage,maxage)),'fig');
%     
    
    
    % combine two maps
    
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

%----------------------------------------------------------------------


figure(1)
%
%  cla;
% scatterhist(y,x,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
%
% figure(2)
%
% scatterhist(y2,x2,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
%
% figure(3)
% scatterhist(y3,x3,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
%
% figure(11)
%  cla;
% scatterhist(yy,xext,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
%
% figure(21)
% scatterhist(yy2,xext2,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
%
% figure(31)
% scatterhist(yy3,xextsi,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
% % x35=xext(Id35);
% y35=y(Id35);
% Id69=find(Ag>=6&Ag<=9);
% x69=xext(Id69);
% y69=y(Id69);
% Id10=find(Ag>=10&Ag<=15);
% x10=xext(Id10);
% y10=y(Id10);
% figure(4)
%  cla;
% ecdf(x35)
% hold on
% ecdf(x69)
% ecdf(x10)
% figure(10)
% scatterhist(y10,x10)
% title(sprintf('#Cell=%i',length(y10)))
% figure(69)
% scatterhist(y69,x69)
% title(sprintf('#Cell=%i',length(y69)))
% figure(35)
% scatterhist(y35,x35)
% title(sprintf('#Cell=%i',length(y35)))
%
% Id35=find(Ag>=3&Ag<=5);
% x35=xext2(Id35);
% y35=y2(Id35);
% Id69=find(Ag>=6&Ag<=9);
% x69=xext2(Id69);
% y69=y2(Id69);
% Id10=find(Ag>=10&Ag<=15);
% x10=xext2(Id10);
% y10=y2(Id10);
% figure(5)
%  cla;
% ecdf(x35)
% hold on
% ecdf(x69)
% ecdf(x10)
% x35=xextsi(Id35);
% y35=y3(Id35);
% Id69=find(Ag>=6&Ag<=9);
% x69=xextsi(Id69);
% y69=y3(Id69);
% Id10=find(Ag>=10&Ag<=15);
% x10=xextsi(Id10);
% y10=y3(Id10);
% figure(6)
%  cla;
%
%
%
% ecdf(x35);
%
% hold on
% ecdf(x69)
% ecdf(x10)