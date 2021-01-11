function [ output_args ] = COMPARE_CONTROLVsObject(Control_path, Object_path,savepath,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%type=1; %use ranksum 
folderMeanpath{1}=Control_path;
folderMeanpath{2}=Object_path;






for cmpind=1:2   
a=dir(fullfile(folderMeanpath{cmpind},'*.mat'));
Cellfile=size(a);

Pchglayer60=[];
Ppeaklayer60=[];
Pchglayer10=[];
Ppeaklayer10=[]; 

chglayer60=[];
peaklayer60=[];
chglayer10=[];
peaklayer10=[];


CorrHighMg60=[];
CorrHighMg10=[];

arealayer10=[];
arealayer60=[];
DirArea=[];
Ag=[];
D60Pk=[];
Bnd=[];
D60chg=[];
D10Pk=[];
D10chg=[];
n=0;
N=[];
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
                
%                 Soma= Cellrecord{fn}.SomaCoordinates;
%                 flipimg=Cellrecord{fn}.flipimg;
%                 pth=char(Cellrecord{fn}.Pth);
%                 stimcoordinates=Cellrecord{fn}.StimCoordinates;
%                 rotate_angle=Cellrecord{fn}.SpatialRotation;
%                 [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
%                 X=stimCoordinates(1,:);
%                 Y=stimCoordinates(2,:);
%                 Xaxis{N}=X;
%                 Yaxis{N}=Y;
%                 Zpk60=zeros(size(X)).*NaN;
%                 Zchg60=zeros(size(X)).*NaN;
%                 Zden60=zeros(size(X));
%                 Zpk10=zeros(size(X)).*NaN;
%                 Zden10=zeros(size(X));
%                 Zchg10=zeros(size(X)).*NaN;
%                 Zlat60=zeros(size(X)).*NaN;
%                 Zlat10=zeros(size(X)).*NaN;
%                 ZlatD60=zeros(size(X)).*NaN;
%                 ZchgD60=zeros(size(X)).*NaN;
%                 ZdenD60=zeros(size(X));
%                   
%                 ZlatD10=zeros(size(X)).*NaN;
%                 Zpksi=zeros(size(X)).*NaN;
%                 Zchgsi=zeros(size(X)).*NaN;
%                 Zdensi=zeros(size(X));
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
%                     k60=k60+1;
                    D60chg=[D60chg;Cellrecord{fn}.DistWidthChargeTp1Hp60];
                    D60Pk=[D60Pk;Cellrecord{fn}.DistWidthPeakTp1Hp60];
                    Pchglayer60=[Pchglayer60;Cellrecord{fn}.PchgTp1Hp60];
                    Ppeaklayer60=[Ppeaklayer60;Cellrecord{fn}.PpeakTp1Hp60];
                    chglayer60=[chglayer60;Cellrecord{fn}.chgTp1Hp60];
                    peaklayer60=[peaklayer60;Cellrecord{fn}.peakTp1Hp60];
                    arealayer60=[arealayer60;Cellrecord{fn}.AreaTp1Hp60];
                    LdistPk60{n}=Cellrecord{fn}.LaydistPkTp1Hp60;
                    LdistChg60{n}=Cellrecord{fn}.LaydistChgTp1Hp60;

                     Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
                    DirArea=[DirArea,length(Indden).*Cellrecord{fn}.PatternSpacing(1).^2];
%                     Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg60));
%                     Zpk60(Indpk)=Cellrecord{fn}.meanpeakHighMg60(Indpk);
%                     Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg60));
%                     Zchg60(Indchg)=Cellrecord{fn}.meanareaHighMg60(Indchg);
%                      Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg60));
%                     Zlat60(Indlat)=Cellrecord{fn}.meanlatencyHighMg60(Indlat);
%                     Indden=find((Cellrecord{fn}.meanflagHighMg60>0));
%                     Zden60(Indden)=Cellrecord{fn}.meanflagHighMg60(Indden);
%                     
%                   
%                     Indchg=find(~isnan(Cellrecord{fn}.meanDirareaHighMg60));
%                     ZchgD60(Indchg)=Cellrecord{fn}.meanDirareaHighMg60(Indchg);
%                      Indlat=find(~isnan(Cellrecord{fn}.meanDirlatencyHighMg60));
%                     ZlatD60(Indlat)=Cellrecord{fn}.meanDirlatencyHighMg60(Indlat);
%                     Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
%                     ZdenD60(Indden)=Cellrecord{fn}.meanflagDHighMg60(Indden);
                    
%                     
%                     StimcorrX60{k60}=X;
%                     StimcorrY60{k60}=Y;
%                     ChgCell60{k60}=Zchg60;
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
                     LdistPk60{n}=NaN;
                     LdistChg60{n}=NaN;
                end
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
%                     k10=k10+1;
                    D10chg=[D10chg;Cellrecord{fn}.DistWidthChargeTp1Hp10];
                    D10Pk=[D10Pk;Cellrecord{fn}.DistWidthPeakTp1Hp10];
                    Pchglayer10=[Pchglayer10;Cellrecord{fn}.PchgTp1Hp10];
                    Ppeaklayer10=[Ppeaklayer10;Cellrecord{fn}.PpeakTp1Hp10];
                     chglayer10=[chglayer10;Cellrecord{fn}.chgTp1Hp10];
                    peaklayer10=[peaklayer10;Cellrecord{fn}.peakTp1Hp10];
                    arealayer10=[arealayer10;Cellrecord{fn}.AreaTp1Hp10];
                     LdistPk10{n}=Cellrecord{fn}.LaydistPkTp1Hp10;
                    LdistChg10{n}=Cellrecord{fn}.LaydistChgTp1Hp10;
%                     Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg10));
%                     Zpk10(Indpk)=Cellrecord{fn}.meanpeakHighMg10(Indpk);
%                     Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg10));
%                     Zchg10(Indchg)=Cellrecord{fn}.meanareaHighMg10(Indchg);
%                      Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg10));
%                     Zlat10(Indlat)=Cellrecord{fn}.meanlatencyHighMg10(Indlat);
%                     Indden=find((Cellrecord{fn}.meanflagHighMg10>0));
%                     Zden10(Indden)=Cellrecord{fn}.meanflagHighMg10(Indden);
                     CorrHighMg10=[CorrHighMg10,Cellrecord{fn}.CorrflagHighMg10];
                
%                     StimcorrX10{k10}=X;
%                     StimcorrY10{k10}=Y;
%                     ChgCell10{k10}=Zchg10;
                    
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
                     LdistPk10{n}=NaN;
                     LdistChg10{n}=NaN;
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
%                 Zpk60avg{Nsi}=Zpk60;
%                 Zchg60avg{Nsi}=Zchg60;
%                 Zden60avg{Nsi}=Zden60;
%                 Zpk10avg{Nsi}=Zpk10;
%                 Zchg10avg{Nsi}=Zchg10;
%                 Zden10avg{Nsi}=Zden10;
%                 Zlat10avg{Nsi}=Zlat10;
%                 Zlat60avg{Nsi}=Zlat60;
%                   
%                 ZchgD60avg{Nsi}=ZchgD60;
%                 ZdenD60avg{Nsi}=ZdenD60;
%               
%                 ZlatD60avg{Nsi}=ZlatD60;
%                 Zpksiavg{Nsi}=Zpksi;
%                 Zchgsiavg{Nsi}=Zchgsi;
%                 Zdensiavg{Nsi}=Zdensi;
            end
            
            
            
        end
    end
end
Cs{cmpind}.Pchglayer60=Pchglayer60;
Cs{cmpind}.Ppeaklayer60=Ppeaklayer60;
Cs{cmpind}.Pchglayer10=Pchglayer10;
Cs{cmpind}.Ppeaklayer10=Ppeaklayer10; 

Cs{cmpind}.chglayer60=chglayer60;
Cs{cmpind}.peaklayer60=peaklayer60;
Cs{cmpind}.chglayer10=chglayer10;
Cs{cmpind}.peaklayer10=peaklayer10;


Cs{cmpind}.CorrHighMg60=CorrHighMg60;
Cs{cmpind}.CorrHighMg10=CorrHighMg10;

Cs{cmpind}.arealayer10=arealayer10;
Cs{cmpind}.arealayer60=arealayer60;
Cs{cmpind}.DirArea=DirArea;

Cs{cmpind}.D60Pk=D60Pk;
Cs{cmpind}.Bnd=Bnd;
Cs{cmpind}.D60chg=D60chg;
Cs{cmpind}.D10Pk=D10Pk;
Cs{cmpind}.D10chg=D10chg;
Cs{cmpind}.LdistPk60=LdistPk60;
Cs{cmpind}.LdistPk10=LdistPk10;
Cs{cmpind}.LdistChg60=LdistChg60;
Cs{cmpind}.LdistChg10=LdistChg10;
N=[N;n]
end



Pchg601=Cs{1}.Pchglayer60;
Pchg602=Cs{2}.Pchglayer60;
Ppk601=Cs{1}.Ppeaklayer60;
Ppk602=Cs{2}.Ppeaklayer60;

Pchg101=Cs{1}.Pchglayer10;
Pchg102=Cs{2}.Pchglayer10;
Ppk101=Cs{1}.Ppeaklayer10;
Ppk102=Cs{2}.Ppeaklayer10;

chg601=Cs{1}.chglayer60;
chg602=Cs{2}.chglayer60;
pk601=Cs{1}.peaklayer60;
pk602=Cs{2}.peaklayer60;

chg101=Cs{1}.chglayer10;
chg102=Cs{2}.chglayer10;
pk101=Cs{1}.peaklayer10;
pk102=Cs{2}.peaklayer10;

A601=Cs{1}.arealayer60;
A101=Cs{1}.arealayer10;

A602=Cs{2}.arealayer60;
A102=Cs{2}.arealayer10;


s=size(Pchg601,2);
if sum(nansum(Pchg601(:,end)))<5
    s=s-1;
end
%---------------------------------------------------------------------
for i=1:s
    [ h,p] =StAnalysis( Pchg601(:,i),Pchg602(:,i),type);
    currfig=find_figure('PercentageCharg_60');
    subplot(3,2,i)
    hold on
    h1=cdfplot(Pchg601(:,i));
    hold on
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    
    h2=cdfplot(Pchg602(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
    ii=i;
    if i==2
        ii=23;
    elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('CdfPchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPchg%i.eps',60)));
for i=1:s
    [ h,p] =StAnalysis( Pchg101(:,i),Pchg102(:,i),type);
     currfig=find_figure('PercentageCharg_10');
    subplot(3,2,i)
    hold on
    h1=cdfplot(Pchg101(:,i));
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    hold on
    h2=cdfplot(Pchg102(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
     ii=i;
    if i==2
        ii=23;
        elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('CdfPchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPchg%i.eps',10)));
for i=1:s
    [ h,p] =StAnalysis( Ppk601(:,i),Ppk602(:,i),type);
    currfig=find_figure('PercentagePeak_60');
    subplot(3,2,i)
    hold on
    h1=cdfplot(Ppk601(:,i));
    hold on
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    
    h2=cdfplot(Ppk602(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
    ii=i;
    if i==2
        ii=23;
    elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('CdfPpk%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPpk%i.eps',60)));
for i=1:s
    [ h,p] =StAnalysis( Ppk101(:,i),Ppk102(:,i),type);
     currfig=find_figure('PercentagePeak_10');
    subplot(3,2,i)
    hold on
    h1=cdfplot(Ppk101(:,i));
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    hold on
    h2=cdfplot(Ppk102(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
     ii=i;
    if i==2
        ii=23;
        elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('CdfPpk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPpk%i.eps',10)));
%----------------------------------------------------------------------------------

for i=1:s
    [ h,p] =StAnalysis( chg601(:,i),chg602(:,i),type);
    currfig=find_figure('Charg_60');
    subplot(3,2,i)
    hold on
    h1=cdfplot(chg601(:,i));
    hold on
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    
    h2=cdfplot(chg602(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
    ii=i;
    if i==2
        ii=23;
    elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('Cdfchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfchg%i.eps',60)));
for i=1:s
    [ h,p] =StAnalysis( chg101(:,i),chg102(:,i),type);
     currfig=find_figure('Charg_10');
    subplot(3,2,i)
    hold on
    h1=cdfplot(chg101(:,i));
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    hold on
    h2=cdfplot(chg102(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
     ii=i;
    if i==2
        ii=23;
        elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('Cdfchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfchg%i.eps',10)));
for i=1:s
    [ h,p] =StAnalysis( pk601(:,i),pk602(:,i),type);
    currfig=find_figure('Peak_60');
    subplot(3,2,i)
    hold on
    h1=cdfplot(pk601(:,i));
    hold on
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    
    h2=cdfplot(pk602(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
    ii=i;
    if i==2
        ii=23;
    elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('Cdfpk%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfpk%i.eps',60)));
for i=1:s
    [ h,p] =StAnalysis( pk101(:,i),pk102(:,i),type);
     currfig=find_figure('Peak_10');
    subplot(3,2,i)
    hold on
    h1=cdfplot(pk101(:,i));
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    hold on
    h2=cdfplot(pk102(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
     ii=i;
    if i==2
        ii=23;
        elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('Cdfpk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfpk%i.eps',10)));

%-------------------------------------------------------------------------

for i=1:s
    [ h,p] =StAnalysis( A601(:,i),A602(:,i),type);
    currfig=find_figure('Area_60');
    subplot(2,2,i)
    hold on
    h1=cdfplot(A601(:,i));
    hold on
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    
    h2=cdfplot(A602(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
    ii=i;
    if i==2
        ii=23;
    elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('CdfArea%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfArea%i.eps',60)));

for i=1:s
    [ h,p] =StAnalysis( A101(:,i),A102(:,i),type);
     currfig=find_figure('Area_10');
    subplot(2,2,i)
    hold on
    h1=cdfplot(A101(:,i));
    set(h1,'Color',[0.1,0.5,0.1],'LineWidth',1)
    hold on
    h2=cdfplot(A102(:,i));
    set(h2,'Color',[0,0,0],'LineWidth',1)
     ii=i;
    if i==2
        ii=23;
        elseif i==3
        ii=4;
    elseif i==4
        ii=56;
    end
    title(sprintf('L%d, h=%i, p=%i',ii,h,p))
    legend('NR','DE')
    box off
end
saveas(currfig,fullfile(savepath,sprintf('CdfArea%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfArea%i.eps',10)));

%___________________________________________________________________________
Bn1=Cs{1}.Bnd;
Bn2=Cs{2}.Bnd;

Rp1=deleteoutliers(abs(Bn1(:,3)./abs(Bn1(:,3)-Bn1(:,2))),0.05);
Rp2=deleteoutliers(abs(Bn2(:,3)./abs(Bn2(:,3)-Bn2(:,2))),0.05);
[h,p]=StAnalysis( Rp1,Rp2,type);

Xrp1=ones(size(Rp1));
Xrp2=ones(size(Rp2)).*2;
currfig=find_figure('relative positon of cell in layer2/3')
plot(Xrp1,Rp1,'.g',Xrp2,Rp2,'*b')
xlim([0,3])
box off
title(sprintf('relative position of cell in L23, h=%i, p=%i',h,p))
saveas(currfig,fullfile(savepath,'RelativePosition.fig'),'fig');
print(currfig,'-depsc2',fullfile(savepath,'RelativePosition.eps'));

