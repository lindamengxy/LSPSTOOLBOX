function [ output_args ] = COMPARE_CONTROLVsObject_MultiGroup_NoiseL23(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%type=1; %use ranksum 
%'/Users/lindameng/Dropbox (Kanoldlab)/silent_122313/meanValueeachcell_eventwindow50_direct8'
%'/Volumes/Samsung_T1/AnalysisNOise/meanValueeachcell_eventwindow50_direct8'
%COMPARE_CONTROLVsObject_MultiGroup_NoiseL23('/Users/xymeng/Dropbox (Kanoldlab)/NoiseBinghan/NoiseTillP32singlemapanalysis/NoiseRec',1,30,40,{'Control', 'Noise'},'/Users/xymeng/Dropbox (Kanoldlab)/NoiseBinghan/NoiseTillP32singlemapanalysis/P39NoiseRecControl/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Dropbox (Kanoldlab)/NoiseBinghan/NoiseTillP32singlemapanalysis/P39NoiseRec/meanValueeachcell_eventwindow50_direct8')
%COMPARE_CONTROLVsObject_MultiGroup_subplate('/Volumes/Samsung_T1/HI/cdfs_sp_P5T10',2,'/Volumes/Samsung_T1/HI/CTL_SHAM_SP_P5to10/meanValueeachcell_eventwindow50_direct8', '/Volumes/Samsung_T1/HI/HI_Mild_caut_SP_p5to10/sp_HI-cauterized/meanValueeachcell_eventwindow50_direct8','/Volumes/Samsung_T1/HI/Hypoxia_Sham_SP_P5to10/meanValueeachcell_eventwindow50_direct8','/Volumes/Samsung_T1/HI/HI_Ligated_Severe_SP_p5to10/meanValueeachcell_eventwindow50_direct8')
%COMPARE_CONTROLVsObject_MultiGroup_subplate_Noise('/Volumes/Samsung_T1 1/AnalysisDeaf',2,6,8, '/Users/xymeng/Dropbox (Kanoldlab)/silent_122313/meanValueeachcell_eventwindow50_direct8', '/Volumes/Samsung_T1 1/TMC/DKO_ALLData/meanValueeachcell_eventwindow50_direct8','/Volumes/Samsung_T1 1/AnalysisOTOF/meanValueeachcell_eventwindow50_direct8')
%COMPARE_CONTROLVsObject_Noise23('/Volumes/Samsung_T5/Noise_Binghan/AnalysisControl/meanValueeachcell_eventwindCOMPAREow50_direct8', '/Volumes/Samsung_T5/Noise_Binghan/AnalysisNoise/meanValueeachcell_eventwindow50_direct8',)
% run : close all; COMPARE_CONTROLVsObject_MultiGroup_subplate('/Users/xymeng/Documents/HI_Aminah/SP/AnalyzedData',2,'/Users/xymeng/Documents/HI_Aminah/SP/CtlSham/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Documents/HI_Aminah/SP/HL/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/HI_Aminah/SP/HS/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/HI_Aminah/SP/CtlUnmanip/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/HI_Aminah/SP/HICauterized/meanValueeachcell_eventwindow50_direct8')
% if nargin==6
% folderMeanpath{1}=Control_path;
% folderMeanpath{2}=Object_path1;
% Ninput=2;
% elseif nargin==7
%  folderMeanpath{1}=Control_path;
% folderMeanpath{2}=Object_path1;
% folderMeanpath{3}=Object_path2;
% Ninput=3;
% elseif nargin==8
%      folderMeanpath{1}=Control_path;
% folderMeanpath{2}=Object_path1;
% folderMeanpath{3}=Object_path2;
% folderMeanpath{4}=Object_path3;
% Ninput=4;
% elseif nargin==9
%      folderMeanpath{1}=Control_path;
% folderMeanpath{2}=Object_path1;
% folderMeanpath{3}=Object_path2;
% folderMeanpath{4}=Object_path3;
% folderMeanpath{5}=Object_path4;
% Ninput=5;
% end

%COMPARE_CONTROLVsObject_MultiGroup_L23Development('/Volumes/Samsung_T1/AnalysisDevelopementL23/CDFplot',1,'/Volumes/Samsung_T1/AnalysisDevelopementL23/meanValueeachcell_eventwindow50_direct8',5,6,7,9,10,15)




close all;
F=1;Agmin=[];Agmax=[];


    if length(varargin)<5
        F=0;
       
    end



Ninput=0;

if F
    savepath=varargin{1};
    type=varargin{2};
    Agmin=varargin{3};
    Agmax=varargin{4};
    GN=varargin{5};
    for i=6:length(varargin)
   folderMeanpath{i-5}=varargin{i};
   Ninput=Ninput+1;
    end
end

N=[]
savepath=fullfile(savepath,sprintf('P%dT%d',Agmin,Agmax));
if ~isdir(savepath)
    mkdir(savepath)
end
for cmpind=1:Ninput 
a=dir(fullfile(folderMeanpath{cmpind},'*.mat'));
Cellfile=size(a);
Parealayer60=[];
Pchglayer60=[];
Ppeaklayer60=[];
Pchglayer10=[];
Ppeaklayer10=[]; 
Parealayer10=[]; 
chglayer60=[];
peaklayer60=[];
chglayer10=[];

peaklayer10=[];

Mchglayer60=[];
Mpeaklayer60=[];
Mchglayer10=[];
Mpeaklayer10=[];

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

D60dens=[];
D10dens=[];
Dist8060L=[];
Dist8010L=[];

 ExpSum60=[];
    ChgExpSum60=[];
    PkExpSum60=[];
    
      ExpSum10=[];
    ChgExpSum10=[];
    PkExpSum10=[];
   
for i= 1:1:Cellfile
    if ~strncmp(a(i).name,'.',1)
    Cg= load(fullfile(folderMeanpath{cmpind},a(i).name));
    if isfield(Cg,'Cellrecord')
        Cellrecord=Cg.Cellrecord;
       for fn=1:Cg.Num_cell
           
            if (isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')||isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10'))&&Cellrecord{fn}.age>=Agmin&&Cellrecord{fn}.age<=Agmax
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
                    NL=length(Cellrecord{fn}.DistWidthPeakTp1Hp60);
                else NL=length(Cellrecord{fn}.DistWidthPeakTp1Hp10);
                end
                
                
                Ag=[Ag;Cellrecord{fn}.age];
                
                if length(Cellrecord{fn}.Boundry)>5
                        Bd(1:4)=Cellrecord{fn}.Boundry(1:4);
                        Bd(5)=Cellrecord{fn}.Boundry(6);
                    else Bd=Cellrecord{fn}.Boundry;
                    end
                
                Bnd=[Bnd;Bd];
               n=n+1;
               Boundry{n}=Cellrecord{fn}.Boundry;
                
               rb=abs(Cellrecord{fn}.Boundry(1,2)/abs(Cellrecord{fn}.Boundry(1,2)-Cellrecord{fn}.Boundry(1,1)));
               if rb>1
                   
                   Cellrecord{fn}.Pth
               end
%                
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
%                ZdenD60=zeros(size(X));
                  
%                 ZlatD10=zeros(size(X)).*NaN;
%                 Zpksi=zeros(size(X)).*NaN;
%                 Zchgsi=zeros(size(X)).*NaN;
%                 Zdensi=zeros(size(X));
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
%                     k60=k60+1;
                        D60chg=[D60chg;Cellrecord{fn}.DistWidthChargeTp1Hp60];
                        D60dens=[D60dens;Cellrecord{fn}.DistWidthdensTp1Hp60];
                        
                          
                    D60Pk=[D60Pk;Cellrecord{fn}.DistWidthPeakTp1Hp60];
                    Pchglayer60=[Pchglayer60;Cellrecord{fn}.PchgTp1Hp60];
                    Ppeaklayer60=[Ppeaklayer60;Cellrecord{fn}.PpeakTp1Hp60];
                      Parealayer60=[Parealayer60;Cellrecord{fn}.PareaTp1Hp60];
                    chglayer60=[chglayer60;Cellrecord{fn}.chgTp1Hp60];
                    peaklayer60=[peaklayer60;Cellrecord{fn}.peakTp1Hp60];
                    arealayer60=[arealayer60;Cellrecord{fn}.AreaTp1Hp60];
                     Mchglayer60=[Mchglayer60;Cellrecord{fn}.MchgTp1Hp60];
                    Mpeaklayer60=[Mpeaklayer60;Cellrecord{fn}.MpeakTp1Hp60];
                    LdistPk60{n}=Cellrecord{fn}.LaydistPkTp1Hp60;
                    LdistChg60{n}=Cellrecord{fn}.LaydistChgTp1Hp60;
Dist8060L=[Dist8060L;Cellrecord{fn}.Dist8060LayerTp1Hp60];
                     Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
                    DirArea=[DirArea;length(Indden).*Cellrecord{fn}.PatternSpacing(1).^2];
%                     Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg60));
%                     Zpk60(Indpk)=Cellrecord{fn}.meanpeakHighMg60(Indpk);
%                     Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg60));
%                     Zchg60(Indchg)=Cellrecord{fn}.meanareaHighMg60(Indchg);
%                      Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg60));
%                     Zlat60(Indlat)=Cellrecord{fn}.meanlatencyHighMg60(Indlat);
%                       Indden=find((Cellrecord{fn}.meanflagHighMg60>0));
%                     Zden60(Indden)=Cellrecord{fn}.meanflagHighMg60(Indden);
%                     
%                   
%                     Indchg=find(~isnan(Cellrecord{fn}.meanDirareaHighMg60));
%                     ZchgD60(Indchg)=Cellrecord{fn}.meanDirareaHighMg60(Indchg);
%                      Indlat=find(~isnan(Cellrecord{fn}.meanDirlatencyHighMg60));
%                     ZlatD60(Indlat)=Cellrecord{fn}.meanDirlatencyHighMg60(Indlat);
%                      DirArea=[DirArea;sum((Cellrecord{fn}.meanflagDHighMg60>0))];
                  
                     ExpSum60=[ExpSum60;Cellrecord{fn}.expSum60];
    ChgExpSum60=[ChgExpSum60;Cellrecord{fn}.ChgexpSum60];
    PkExpSum60=[PkExpSum60;Cellrecord{fn}.PkexpSum60];
                        
%                     
%                     StimcorrX60{k60}=X;
%                     StimcorrY60{k60}=Y;
%                     ChgCell60{k60}=Zchg60;
                   CorrHighMg60=[CorrHighMg60,Cellrecord{fn}.CorrflagHighMg60];
                     LaminaArea60{n}=Cellrecord{fn}.LaminaAreaTp1Hp60;
                     LaminaPeak60{n}=Cellrecord{fn}.LaminaPeakTp1Hp60;
                      LaminaChg60{n}=Cellrecord{fn}.LaminaChgTp1Hp60;
                          ColumaArea60{n}=Cellrecord{fn}.ColumaAreaTp1Hp60;
                     ColumaPeak60{n}=Cellrecord{fn}.ColumaPeakTp1Hp60;
                      ColumaChg60{n}=Cellrecord{fn}.ColumaChgTp1Hp60;
                else  Nanmaxtric= ones(1,NL).*NaN;
                    D60chg=[D60chg; Nanmaxtric];
                    D60Pk=[D60Pk; Nanmaxtric];
                     D60dens=[D60dens; Nanmaxtric];
                    Nanmaxtric2=ones(1,3).*NaN;
                    ExpSum60=[ExpSum60;Nanmaxtric2];
                        ChgExpSum60=[ChgExpSum60;Nanmaxtric2];
                        PkExpSum60=[PkExpSum60;Nanmaxtric2];
                        
                    Pchglayer60=[Pchglayer60;Nanmaxtric2];
                    Parealayer60=[Parealayer60;Nanmaxtric2];
                    Ppeaklayer60=[Ppeaklayer60;Nanmaxtric2];
                    chglayer60=[chglayer60;Nanmaxtric2];
                    peaklayer60=[peaklayer60;Nanmaxtric2];
                    Mchglayer60=[Mchglayer60;Nanmaxtric2];
                    Mpeaklayer60=[Mpeaklayer60;Nanmaxtric2];
                    arealayer60=[arealayer60;Nanmaxtric2];
                     CorrHighMg60=[CorrHighMg60,NaN];
                     LdistPk60{n}=NaN;
                     LdistChg60{n}=NaN;
                      Dist8060L=[Dist8060L;Nanmaxtric2];
                         LaminaArea60{n}=NaN;
                     LaminaPeak60{n}=NaN;
                      LaminaChg60{n}=NaN;
                         ColumaArea60{n}=NaN;
                     ColumaPeak60{n}=NaN;
                      ColumaChg60{n}=NaN;
                end
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
%                     k10=k10+1;
                    D10chg=[D10chg;Cellrecord{fn}.DistWidthChargeTp1Hp10];
                    D10Pk=[D10Pk;Cellrecord{fn}.DistWidthPeakTp1Hp10];
                      D10dens=[D10dens;Cellrecord{fn}.DistWidthdensTp1Hp10];
                    Pchglayer10=[Pchglayer10;Cellrecord{fn}.PchgTp1Hp10];
                    Ppeaklayer10=[Ppeaklayer10;Cellrecord{fn}.PpeakTp1Hp10];
                    Parealayer10=[Parealayer10;Cellrecord{fn}.PareaTp1Hp10];
                     chglayer10=[chglayer10;Cellrecord{fn}.chgTp1Hp10];
                    peaklayer10=[peaklayer10;Cellrecord{fn}.peakTp1Hp10];
                    arealayer10=[arealayer10;Cellrecord{fn}.AreaTp1Hp10];
                     Mchglayer10=[Mchglayer10;Cellrecord{fn}.MchgTp1Hp10];
                    Mpeaklayer10=[Mpeaklayer10;Cellrecord{fn}.MpeakTp1Hp10];
                     LdistPk10{n}=Cellrecord{fn}.LaydistPkTp1Hp10;
                    LdistChg10{n}=Cellrecord{fn}.LaydistChgTp1Hp10;
                    Dist8010L=[Dist8010L;Cellrecord{fn}.Dist8010LayerTp1Hp10];
%                     Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg10));
%                     Zpk10(Indpk)=Cellrecord{fn}.meanpeakHighMg10(Indpk);
%                     Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg10));
%                     Zchg10(Indchg)=Cellrecord{fn}.meanareaHighMg10(Indchg);
%                      Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg10));
%                     Zlat10(Indlat)=Cellrecord{fn}.meanlatencyHighMg10(Indlat);
%                     Indden=find((Cellrecord{fn}.meanflagHighMg10>0));
%                     Zden10(Indden)=Cellrecord{fn}.meanflagHighMg10(Indden);
                     CorrHighMg10=[CorrHighMg10,Cellrecord{fn}.CorrflagHighMg10];
                 ExpSum10=[ExpSum10;Cellrecord{fn}.expSum10];
                        ChgExpSum10=[ChgExpSum10;Cellrecord{fn}.ChgexpSum10];
                        PkExpSum10=[PkExpSum10;Cellrecord{fn}.PkexpSum10];
%                     StimcorrX10{k10}=X;
%                     StimcorrY10{k10}=Y;
%                     ChgCell10{k10}=Zchg10;
 LaminaArea10{n}=Cellrecord{fn}.LaminaAreaTp1Hp10;
                     LaminaPeak10{n}=Cellrecord{fn}.LaminaPeakTp1Hp10;
                      LaminaChg10{n}=Cellrecord{fn}.LaminaChgTp1Hp10;
                       ColumaArea10{n}=Cellrecord{fn}.ColumaAreaTp1Hp10;
                     ColumaPeak10{n}=Cellrecord{fn}.ColumaPeakTp1Hp10;
                      ColumaChg10{n}=Cellrecord{fn}.ColumaChgTp1Hp10;
                    
                else  Nanmaxtric= ones(1,NL).*NaN;
                    D10chg=[D10chg; Nanmaxtric];
                    D10Pk=[D10Pk; Nanmaxtric];
                    D10dens=[D10dens; Nanmaxtric];
                    Nanmaxtric2=ones(1,3).*NaN;
                    Pchglayer10=[Pchglayer10;Nanmaxtric2];
                     Parealayer10=[Parealayer10;Nanmaxtric2];
                    Ppeaklayer10=[Ppeaklayer10;Nanmaxtric2];
                     arealayer10=[arealayer10;Nanmaxtric2];
                        chglayer10=[chglayer10;Nanmaxtric2];
                    peaklayer10=[peaklayer10;Nanmaxtric2];
                     Mchglayer10=[Mchglayer10;Nanmaxtric2];
                    Mpeaklayer10=[Mpeaklayer10;Nanmaxtric2];
                    CorrHighMg10=[CorrHighMg10,NaN];
                     LdistPk10{n}=NaN;
                     LdistChg10{n}=NaN;
                     Dist8010L=[Dist8010L;Nanmaxtric2];
                        ExpSum10=[ExpSum10;Nanmaxtric2];
                        ChgExpSum10=[ChgExpSum10;Nanmaxtric2];
                        PkExpSum10=[PkExpSum10;Nanmaxtric2];
                       LaminaArea10{n}=NaN;
                     LaminaPeak10{n}=NaN;
                      LaminaChg10{n}=NaN;
                      ColumaArea10{n}=NaN;
                     ColumaPeak10{n}=NaN;
                      ColumaChg10{n}=NaN;
                      
                      
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
end
Cs{cmpind}.Pchglayer60=Pchglayer60;
Cs{cmpind}.Ppeaklayer60=Ppeaklayer60;
Cs{cmpind}.Pchglayer10=Pchglayer10;
Cs{cmpind}.Ppeaklayer10=Ppeaklayer10; 
Cs{cmpind}.Parealayer60=Parealayer60;
Cs{cmpind}.Parealayer10=Parealayer10;


Cs{cmpind}.chglayer60=chglayer60;
Cs{cmpind}.peaklayer60=peaklayer60;                
Cs{cmpind}.chglayer10=chglayer10;
Cs{cmpind}.peaklayer10=peaklayer10;

Cs{cmpind}.Mchglayer60=Mchglayer60;
Cs{cmpind}.Mpeaklayer60=Mpeaklayer60;                        
Cs{cmpind}.Mchglayer10=Mchglayer10;
Cs{cmpind}.Mpeaklayer10=Mpeaklayer10;

Cs{cmpind}.CorrHighMg60=CorrHighMg60;
Cs{cmpind}.CorrHighMg10=CorrHighMg10;

Cs{cmpind}.arealayer10=arealayer10;
Cs{cmpind}.arealayer60=arealayer60;
Cs{cmpind}.DirArea=DirArea;

Cs{cmpind}.D60Pk=D60Pk;
Cs{cmpind}.D60dens=D60dens;
Cs{cmpind}.Bnd=Bnd;
Cs{cmpind}.Boundry=Boundry;
Cs{cmpind}.D60chg=D60chg;
Cs{cmpind}.D10Pk=D10Pk;
Cs{cmpind}.D10chg=D10chg;
Cs{cmpind}.D10dens=D10dens;
Cs{cmpind}.LdistPk60=LdistPk60;
Cs{cmpind}.LdistPk10=LdistPk10;
Cs{cmpind}.LdistChg60=LdistChg60;
Cs{cmpind}.LdistChg10=LdistChg10;
Cs{cmpind}.Dist8060L=Dist8060L;
Cs{cmpind}.Dist8010L=Dist8010L;

Cs{cmpind}.ExpSum60=ExpSum60;
Cs{cmpind}.ChgExpSum60=ChgExpSum60;
Cs{cmpind}.PkExpSum60=PkExpSum60;
Cs{cmpind}.Ag=Ag;
Cs{cmpind}.ExpSum10=ExpSum10;
Cs{cmpind}.ChgExpSum10=ChgExpSum10;
Cs{cmpind}.PkExpSum10=PkExpSum10;
% GN{cmpind}=sprintf('P%d-%d',Agmin(cmpind),Agmax(cmpind));
Cs{cmpind}.LaminaPeak60= LaminaPeak60;
Cs{cmpind}.LaminaArea60= LaminaArea60;
Cs{cmpind}.LaminaChg60= LaminaChg60;
Cs{cmpind}.LaminaPeak10= LaminaPeak10;
Cs{cmpind}.LaminaArea10= LaminaArea10;
Cs{cmpind}.LaminaChg10= LaminaChg10;
Cs{cmpind}.ColumaPeak60= ColumaPeak60;
Cs{cmpind}.ColumaArea60= ColumaArea60;
Cs{cmpind}.ColumaChg60= ColumaChg60;
Cs{cmpind}.ColumaPeak10= ColumaPeak10;
Cs{cmpind}.ColumaArea10= ColumaArea10;
Cs{cmpind}.ColumaChg10= ColumaChg10;


N=[N;n]
end



%%
LaLength60=[];AgM=[]; LaLength10=[];
for cmpind=1:Ninput
    
    Lamina60=Cs{cmpind}.LaminaArea60;
    Lamina10=Cs{cmpind}.LaminaArea10;
    for ii=1:N(cmpind)
        LaLength60=[LaLength60,length(Lamina60{ii})];
        LaLength10=[LaLength10,length(Lamina10{ii})];
    end
end

mx=max(max(LaLength60),max(LaLength10));
my60=length(find(~isnan(LaLength60)));
my10=length(find(~isnan(LaLength10)));

MatrixArea60=NaN(mx,my60);MatrixPeak60=NaN(mx,my60);MatrixChg60=NaN(mx,my60);
MatrixArea10=NaN(mx,my60);MatrixPeak10=NaN(mx,my60);MatrixChg10=NaN(mx,my60);

MatrixPercArea60=NaN(mx,my60);MatrixPercPeak60=NaN(mx,my60);MatrixPercChg60=NaN(mx,my60);
MatrixPercArea10=NaN(mx,my60);MatrixPercPeak10=NaN(mx,my60);MatrixPercChg10=NaN(mx,my60);
ic=1;ic10=1;
soma=[];
x=20:40:40*mx;soma10=[];
for cmpind=1:Ninput
    Ag=Cs{cmpind}.Ag;
    
    [Ag1 I]=sort(Ag);
    AgM=[AgM;Ag1];
    
    LaminaArea60=Cs{cmpind}.LaminaArea60;
    LaminaArea10=Cs{cmpind}.LaminaArea10;
    
    LaminaPeak60=Cs{cmpind}.LaminaPeak60;
    LaminaPeak10=Cs{cmpind}.LaminaPeak10;
    LaminaChg60=Cs{cmpind}.LaminaChg60;
    LaminaChg10=Cs{cmpind}.LaminaChg10;
    
    
    for ii=1:N(cmpind)
        if length(LaminaArea60{I(ii)})
            MatrixArea60(1:length(LaminaArea60{I(ii)}),ic)=LaminaArea60{I(ii)};
            MatrixPercArea60(1:length(LaminaArea60{I(ii)}),ic)=LaminaArea60{I(ii)}./sum(LaminaArea60{I(ii)});
            MatrixPeak60(1:length(LaminaPeak60{I(ii)}),ic)=LaminaPeak60{I(ii)};
            MatrixPercPeak60(1:length(LaminaPeak60{I(ii)}),ic)=LaminaPeak60{I(ii)}./sum(LaminaPeak60{I(ii)});
            MatrixChg60(1:length(LaminaChg60{I(ii)}),ic)=LaminaChg60{I(ii)};
            MatrixPercChg60(1:length(LaminaChg60{I(ii)}),ic)=LaminaChg60{I(ii)}./sum(LaminaChg60{I(ii)});
            ic=ic+1;
            soma=[soma,-Cs{cmpind}.Boundry{I(ii)}(1)];
            Bd1{ic}=Cs{cmpind}.Boundry{I(ii)}-Cs{cmpind}.Boundry{I(ii)}(1);
            
        end
        
         if length(LaminaArea10{I(ii)})
            MatrixArea10(1:length(LaminaArea10{I(ii)}),ic)=LaminaArea10{I(ii)};
            MatrixPercArea10(1:length(LaminaArea10{I(ii)}),ic)=LaminaArea10{I(ii)}./sum(LaminaArea10{I(ii)});
            MatrixPeak10(1:length(LaminaPeak10{I(ii)}),ic)=LaminaPeak10{I(ii)};
            MatrixPercPeak10(1:length(LaminaPeak10{I(ii)}),ic)=LaminaPeak10{I(ii)}./sum(LaminaPeak10{I(ii)});
            MatrixChg10(1:length(LaminaChg10{I(ii)}),ic)=LaminaChg10{I(ii)};
            MatrixPercChg10(1:length(LaminaChg10{I(ii)}),ic)=LaminaChg10{I(ii)}./sum(LaminaChg10{I(ii)});
            ic10=ic10+1;
            soma10=[soma10,-Cs{cmpind}.Boundry{I(ii)}(1)];
            Bd10{ic}=Cs{cmpind}.Boundry{I(ii)}-Cs{cmpind}.Boundry{I(ii)}(1);
            
        end
        
        
        
        
    end
end


currfig=find_figure('ColumnPercArea60')
imagesc(1:length(soma),x,MatrixPercArea60)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixPercArea60(:)))*4]);
hold on
plot(1:length(soma),soma,'wo')
for ii=1:ic-1
    for jj=1:length(Bd1{ii})
        line([ii-0.5 ii+0.5],[Bd1{ii}(jj) Bd1{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnPercArea60.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnPercArea60.eps'));

   
   currfig=find_figure('ColumnArea60')
imagesc(1:length(soma),x,MatrixArea60)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixArea60(:)))*4]);
hold on
plot(1:length(soma),soma,'wo')
for ii=1:ic-1
    for jj=1:length(Bd1{ii})
        line([ii-0.5 ii+0.5],[Bd1{ii}(jj) Bd1{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnArea60.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnArea60.eps'));


   
   
   
   
 currfig=find_figure('ColumnPercChg60')
imagesc(1:length(soma),x,MatrixPercChg60)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixPercChg60(:)))*4]);
hold on
plot(1:length(soma),soma,'wo')
for ii=1:ic-1
    for jj=1:length(Bd1{ii})
        line([ii-0.5 ii+0.5],[Bd1{ii}(jj) Bd1{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnPercChg60.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnPercChg60.eps'));

   
   currfig=find_figure('ColumnChg60')
imagesc(1:length(soma),x,MatrixChg60)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixChg60(:)))*4]);
hold on
plot(1:length(soma),soma,'wo')
for ii=1:ic-1
    for jj=1:length(Bd1{ii})
        line([ii-0.5 ii+0.5],[Bd1{ii}(jj) Bd1{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnChg60.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnChg60.eps'));  
   
   
   
   
   
 currfig=find_figure('ColumnPercArea10')
imagesc(1:length(soma10),x,MatrixPercArea10)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixPercArea10(:)))*4]);
hold on
plot(1:length(soma10),soma10,'wo')
for ii=1:ic-1
    for jj=1:length(Bd10{ii})
        line([ii-0.5 ii+0.5],[Bd10{ii}(jj) Bd10{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnPercArea10.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnPercArea10.eps'));

   
   currfig=find_figure('ColumnArea10')
imagesc(1:length(soma10),x,MatrixArea10)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixArea10(:)))*4]);
hold on
plot(1:length(soma10),soma10,'wo')
for ii=1:ic-1
    for jj=1:length(Bd10{ii})
        line([ii-0.5 ii+0.5],[Bd10{ii}(jj) Bd10{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnArea10.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnArea10.eps'));


   
   
   
   
 currfig=find_figure('ColumnPercChg10')
imagesc(1:length(soma10),x,MatrixPercChg10)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixPercChg10(:)))*4]);
hold on
plot(1:length(soma10),soma10,'wo')
for ii=1:ic-1
    for jj=1:length(Bd10{ii})
        line([ii-0.5 ii+0.5],[Bd10{ii}(jj) Bd10{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnPercChg10.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnPercChg10.eps'));

   
   currfig=find_figure('ColumnChg10')
imagesc(1:length(soma10),x,MatrixChg10)
colormap('hot')
caxis([0,nanmean(nanmean(MatrixChg10(:)))*4]);
hold on
plot(1:length(soma10),soma10,'wo')
for ii=1:ic-1
    for jj=1:length(Bd10{ii})
        line([ii-0.5 ii+0.5],[Bd10{ii}(jj) Bd10{ii}(jj)],'Color','w')
    end
end
nn=0;
for ii=1:length(N)
    nn=nn+N(ii);
    line([nn+0.5 nn+0.5],[0,40*mx])
end
hold off

saveas(currfig,fullfile(savepath,'ColumnChg10.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'ColumnChg10.eps'));  
  
   
   
   
   


%%
AGE=[];
DensWidth60=[];
ChgWidth60=[];
PkWidth60=[];
DensRatio60=[];
ChgRatio60=[];
PkRatio60=[];
ElementN=size(Cs{1}.arealayer60,2)-1;
for j=1:Ninput
     AGE=[AGE;Cs{j}.Ag];
 DensWidth60=[DensWidth60;Cs{j}.D60dens(:,end-ElementN:end)];
 PkWidth60=[PkWidth60;Cs{j}.D60Pk(:,end-ElementN:end)];
 ChgWidth60=[ChgWidth60;Cs{j}.D60chg(:,end-ElementN:end)];
 DensRatio60=[DensRatio60;Cs{j}.arealayer60(:,2)./Cs{j}.arealayer60(:,1)];
 ChgRatio60=[ChgRatio60;Cs{j}.chglayer60(:,2)./Cs{j}.chglayer60(:,1)];
 PkRatio60=[PkRatio60;Cs{j}.peaklayer60(:,2)./Cs{j}.peaklayer60(:,1)];
end
AGE(AGE>30)=30;

A=unique(AGE);medianDens60=[];medianChg60=[];medianPk60=[];medianDenRatio60=[];medianChgRatio60=[];medianPkRatio60=[];
for ia=1:length(A);
    indA=find(AGE==A(ia));
    medianDens60=[medianDens60;nanmedian(DensWidth60(indA,:))];
    medianChg60=[medianChg60;nanmedian(ChgWidth60(indA,:))];
    medianPk60=[medianPk60;nanmedian(PkWidth60(indA,:))];
    medianDenRatio60=[medianDenRatio60;nanmedian(DensRatio60(indA))];
    medianChgRatio60=[medianChgRatio60;nanmedian(ChgRatio60(indA))];
    medianPkRatio60=[medianPkRatio60;nanmedian(PkRatio60(indA))];
end


find_figure('denswidth60')
for jj=1:3
subplot(1,3,jj)
plot(AGE,DensWidth60(:,jj),'ko')
hold on
plot(A,medianDens60(:,jj),'r*-')
hold off
end
find_figure('Chargewidth60')
for jj=1:3
subplot(1,3,jj)
plot(AGE,ChgWidth60(:,jj),'ko')
hold on
plot(A,medianChg60(:,jj),'r*-')
hold off
end
find_figure('Pkwidth60')
for jj=1:3
subplot(1,3,jj)
plot(AGE,PkWidth60(:,jj),'ko')
hold on
plot(A,medianPk60(:,jj),'r*-')
hold off
end






find_figure('Ratio60')

subplot(1,3,1)
plot(AGE,DensRatio60,'ko')
hold on
plot(A,medianDenRatio60,'r*-')
hold off
title('Area, L4:L23')
subplot(1,3,2)
plot(AGE,ChgRatio60,'ko')
hold on
plot(A,medianChgRatio60,'r*-')
hold off
title('Charge, L4:L23')

subplot(1,3,3)
plot(AGE,PkRatio60,'ko')
hold on
plot(A,medianPkRatio60,'r*-')
hold off
title('Peak, L4:L23')






DensWidth10=[];
ChgWidth10=[];
PkWidth10=[];AGE=[];
for j=1:Ninput
     AGE=[AGE;Cs{j}.Ag];
 DensWidth10=[DensWidth10;Cs{j}.D10dens(:,end-ElementN:end)];
 PkWidth10=[PkWidth10;Cs{j}.D10Pk(:,end-ElementN:end)];
 ChgWidth10=[ChgWidth10;Cs{j}.D10chg(:,end-ElementN:end)];
end


A=unique(AGE);medianDens10=[];medianChg10=[];medianPk10=[];
for ia=1:length(A);
    indA=find(AGE==A(ia));
    medianDens10=[medianDens10;nanmedian(DensWidth10(indA,:))];
    medianChg10=[medianChg10;nanmedian(ChgWidth10(indA,:))];
    medianPk10=[medianPk10;nanmedian(PkWidth10(indA,:))];
end


find_figure('denswidth10')
for jj=1:3
subplot(1,3,jj)
plot(AGE,DensWidth10(:,jj),'ko')
hold on
plot(A,medianDens10(:,jj),'r*-')
hold off
end
find_figure('Chargewidth10')
for jj=1:3
subplot(1,3,jj)
plot(AGE,ChgWidth10(:,jj),'ko')
hold on
plot(A,medianChg10(:,jj),'r*-')
hold off
end
find_figure('Pkwidth10')
for jj=1:3
subplot(1,3,jj)
plot(AGE,PkWidth10(:,jj),'ko')
hold on
plot(A,medianPk10(:,jj),'r*-')
hold off
end















s=size(Cs{1}.Pchglayer60,2);
NGroup=Ninput;

Cl=colormap('hsv');
    nl=floor(size(Cl,1)/NGroup);
%-------------------------------------------------------------       
    MExpsum=[];VarExpsum=[];MExpsum10=[];VarExpsum10=[];
    for jd=1:NGroup
        ExpSum60=Cs{jd}.ExpSum60;
         ExpSum10=Cs{jd}.ExpSum10;
        ExpSum60(isnan(ExpSum60))=0;
        ExpSum10(isnan(ExpSum10))=0;
         mdist=[];vdist=[];mdist10=[];vdist10=[];
         for jlayer=1:size(Dist8060L,2)
             mdist=[mdist nanmean(ExpSum60(:,jlayer))];
             vdist=[vdist nanstd(ExpSum60(:,jlayer))./sqrt(size(ExpSum60,1))];
                  mdist10=[mdist10 nanmean(ExpSum10(:,jlayer))];
             vdist10=[vdist10 nanstd(ExpSum10(:,jlayer))./sqrt(size(ExpSum10,1))];
         end
       MExpsum=[ MExpsum;mdist];
     VarExpsum=[VarExpsum;vdist]; 
       MExpsum10=[ MExpsum10;mdist10];
     VarExpsum10=[VarExpsum10;vdist10]; 
      msfile=sprintf('MeanStd_Expsum60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mdist','vdist','-ascii') 
    
    msfile=sprintf('MeanStd_Expsum10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mdist10','vdist10','-ascii') 
     
    end
      
    
    
    
  
  MExpsum=MExpsum';
  VarExpsum=VarExpsum';
  MExpsum10=MExpsum10';
  VarExpsum10=VarExpsum10';
  
%boundline

currfig=find_figure('Expsum60Boundline')
xa=1:length(MExpsum(1,:));
[l,p] = boundedline(xa,MExpsum(1,:), VarExpsum(1,:), '-g*', xa,MExpsum(2,:), VarExpsum(2,:), '-r*',xa,MExpsum(3,:), VarExpsum(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Expsum60Boundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Expsum60Boundline.eps'));

  

currfig=find_figure(sprintf('Expsum60Layer',i))

  h=bar(MExpsum);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MExpsum,1)
  numbars=size(MExpsum,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MExpsum(:,jj),VarExpsum(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MExpsum,1)
      Ds={}; noNaNLength=[];
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.ExpSum60(:,jj);
          noNaNLength=[noNaNLength;length(find(~isnan(Ds{ii})))];
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if ~length(find(noNaNLength<4))
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
      
       H=sigstar(groups,P);
      else P=NaN;p=NaN;
      end
     pfile=sprintf('statisticP_ExpSum60_G%d_Layer%d.txt',NGroup,jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarExpSum60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarExpSum60_G%d.eps',NGroup)));
    
    currfig=find_figure(sprintf('Expsum10Layer',i))
 h=bar(MExpsum10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MExpsum10,1)
  numbars=size(MExpsum10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MExpsum10(:,jj),VarExpsum10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
for jj=1:size(MExpsum10,1)
      Ds={};
      noNaNLength=[];
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.ExpSum10(:,jj);
          noNaNLength=[noNaNLength;length(find(~isnan(Ds{ii})))];
      end
      
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if ~length(find(noNaNLength<4))
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
      else P=NaN;p=NaN;
      end
     pfile=sprintf('statisticP_Expsum10_G%d_Layer%d.txt',NGroup,jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii')  
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarExpsum10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarExpsum10_G%d.eps',NGroup)));
    
   currfig=find_figure('Expsum60Boundline')
xa=1:length(MExpsum10(1,:));
[l,p] = boundedline(xa,MExpsum10(1,:), VarExpsum10(1,:), '-g*', xa,MExpsum10(2,:), VarExpsum10(2,:), '-r*',xa,MExpsum10(3,:), VarExpsum10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Expsum10Boundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Expsum10Boundline.eps')); 
   
   
%-----------------------------------------------------------------------------    
MDist=[];VarDist=[];MDist10=[];VarDist10=[];
    for jd=1:NGroup
        Dist8060L=Cs{jd}.Dist8060L;Dist8010L=Cs{jd}.Dist8010L;
        Dist8060L(isnan(Dist8060L))=0;
        Dist8010L(isnan(Dist8010L))=0;
         mdist=[];vdist=[];mdist10=[];vdist10=[];
         for jlayer=1:size(Dist8060L,2)
             mdist=[mdist nanmean(Dist8060L(:,jlayer))];
             vdist=[vdist nanstd(Dist8060L(:,jlayer))./sqrt(size(Dist8060L,1))];
                  mdist10=[mdist10 nanmean(Dist8010L(:,jlayer))];
             vdist10=[vdist10 nanstd(Dist8010L(:,jlayer))./sqrt(size(Dist8010L,1))];
         end
      MDist=[MDist;mdist];
     VarDist=[VarDist;vdist]; 
      MDist10=[MDist10;mdist10];
     VarDist10=[VarDist10;vdist10]; 
      msfile=sprintf('MeanStd_Dist8060_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mdist','vdist','-ascii') 
    
    msfile=sprintf('MeanStd_Dist8010_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mdist10','vdist10','-ascii') 
     
    end
      
  
  MDist=MDist';
  VarDist=VarDist';
  MDist10=MDist10';
  VarDist10=VarDist10';
  currfig=find_figure('Dist8060LayerBoundline')
xa=1:length(MDist(1,:));
[l,p] = boundedline(xa,MDist(1,:), VarDist(1,:), '-g*', xa,MDist(2,:), VarDist(2,:), '-r*',xa,MDist(3,:), VarDist(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Dist8060LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Dist8060LayerBoundline.eps'));
  
  currfig=find_figure('Dist8010LayerBoundline')
xa=1:length(MDist(1,:));
[l,p] = boundedline(xa,MDist10(1,:), VarDist10(1,:), '-g*', xa,MDist10(2,:), VarDist10(2,:), '-r*',xa,MDist10(3,:), VarDist10(3,:));
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Dist8010LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Dist8010LayerBoundline.eps'));
  

currfig=find_figure(sprintf('Dist8060Layer',i))

  h=bar(MDist);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MDist,1)
  numbars=size(MDist,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MDist(:,jj),VarDist(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Dist8060L(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Dist8060_G%d_Layer%d.txt',NGroup,jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarDist8060_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarDist8060_G%d.eps',NGroup)));

%---------------------------------------------------------------------
currfig=find_figure(sprintf('Dist8010Layer',i))
 h=bar(MDist10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MDist10,1)
  numbars=size(MDist10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MDist10(:,jj),VarDist10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
for jj=1:size(MDist10,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Dist8010L(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Dist8010_G%d_Layer%d.txt',NGroup,jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii')  
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarDist8010_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarDist8010_G%d.eps',NGroup)));








%---------------------------------------------------------------------

if Ninput==2
        [hh,c]=StAnalysis( Cs{1}.DirArea,Cs{2}.DirArea,type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.DirArea,Cs{2}.DirArea,Cs{3}.DirArea)
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.DirArea,Cs{2}.DirArea,Cs{3}.DirArea,Cs{4}.DirArea)
         elseif Ninput==5
        c =StAnalysis5( Cs{1}.DirArea,Cs{2}.DirArea,Cs{3}.DirArea,Cs{4}.DirArea,Cs{5}.DirArea)
end

MDirArea=[];VarDirArea=[];
    for jd=1:NGroup
        DirArea60=Cs{jd}.DirArea;
         mDirArea=[];vDirArea=[];
         for jlayer=1:size(DirArea60,2)
             mDirArea=[mDirArea nanmean(DirArea60(:,jlayer))];
             vDirArea=[vDirArea nanstd(DirArea60(:,jlayer))./sqrt(length(DirArea60))];
               
         end
      MDirArea=[MDirArea;mDirArea];
     VarDirArea=[VarDirArea;vDirArea]; 
      
      msfile=sprintf('MeanStd_DirectArea60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mDirArea','vDirArea','-ascii') 
    
  
     
    end
      
  
  MDirArea=MDirArea;
  VarDirArea=VarDirArea;
 currfig=find_figure('DirectAreaBoundline')
xa=1:length(MDirArea);
[l,p] = boundedline(xa,MDirArea',VarDirArea', '-k*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'DirectAreaBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'DirectAreaBoundline.eps'));


 currfig=find_figure('BarDirArea_60');
 
  h=bar(MDirArea);
%   for ih=1:NGroup
%       h(ih).FaceColor=Cl(ih.*nl,:,:);
%   end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MDirArea,1)
  numbars=size(MDirArea,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MDirArea(:,jj),VarDirArea(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(DirArea,2)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.DirArea(:,jj);
      end
      XXaxis=XXaxis(:);
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_DirArea_G%d.txt',NGroup);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarDirArea_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarDirArea_G%d.eps',NGroup)));





 currfig=find_figure('DirArea_60');
 
 hold on
 for i=1:Ninput
     subplot(1,2,1)
     h1=cdfplot(Cs{i}.DirArea);
     hold on
      set(h1,'LineWidth',2,'Color',Cl(i.*nl,:,:))
     
     
 end
 legend(GN)
 subplot(1,2,2)
 %# A 5-by-5 matrix of random values from 0 to 1
 sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
 sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
  sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else sgm(1,2)=c(ind,end);
  sgm(2,1)=c(ind,end);
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
title('DirectArea60')
saveas(currfig,fullfile(savepath,sprintf('CdfDirectArea%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfDirectArea%i.eps',60)));




MPchglayer60=[];VarPchglayer60=[];MPchglayer10=[];VarPchglayer10=[];
    for jd=1:NGroup
        Pchglayer60=Cs{jd}.Pchglayer60;Pchglayer10=Cs{jd}.Pchglayer10;
         mPchglayer60=[];vPchglayer60=[];mPchglayer10=[];vPchglayer10=[];
         for jlayer=1:size(Dist8060L,2)
             mPchglayer60=[mPchglayer60 nanmean(Pchglayer60(:,jlayer))];
             vPchglayer60=[vPchglayer60 nanstd(Pchglayer60(:,jlayer))./sqrt(size(Pchglayer60,1))];
                  mPchglayer10=[mPchglayer10 nanmean(Pchglayer10(:,jlayer))];
             vPchglayer10=[vPchglayer10 nanstd(Pchglayer10(:,jlayer))./sqrt(size(Pchglayer10,1))];
         end
      MPchglayer60=[MPchglayer60;mPchglayer60];
     VarPchglayer60=[VarPchglayer60;vPchglayer60]; 
      MPchglayer10=[MPchglayer10;mPchglayer10];
     VarPchglayer10=[VarPchglayer10;vPchglayer10]; 
      msfile=sprintf('MeanStd_Pchglayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mPchglayer60','vPchglayer60','-ascii') 
    
    msfile=sprintf('MeanStd_Pchglayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mPchglayer10','vPchglayer10','-ascii') 
     
    end
      
  
  MPchglayer60=MPchglayer60';
  VarPchglayer60=VarPchglayer60';
  MPchglayer10=MPchglayer10';
  VarPchglayer10=VarPchglayer10';
  
  
  currfig=find_figure('Pchglayer60LayerBoundline')
xa=1:length(MDist(1,:));
[l,p] = boundedline(xa,MPchglayer60(1,:), VarPchglayer60(1,:), '-g*', xa,MPchglayer60(2,:), VarPchglayer60(2,:), '-r*',xa,MPchglayer60(3,:), VarPchglayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Pchglayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Pchglayer60LayerBoundline.eps'));
   
   currfig=find_figure('Pchglayer10LayerBoundline')
xa=1:length(MDist(1,:));
[l,p] = boundedline(xa,MPchglayer10(1,:), VarPchglayer10(1,:), '-g*', xa,MPchglayer10(2,:), VarPchglayer10(2,:), '-r*',xa,MPchglayer10(3,:), VarPchglayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Pchglayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Pchglayer10LayerBoundline.eps'));
   
currfig=find_figure('Pchglayer60Layer_stacked')  
h=bar(MPchglayer60','stacked');

set(gca, 'XTickLabel',GN)

   saveas(currfig,fullfile(savepath,sprintf('BarPchglayer60_G%d_stacked.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarPchglayer60_G%d_stacked.eps',NGroup)));
currfig=find_figure(sprintf('Pchglayer60Layer',i))

  h=bar(MPchglayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MPchglayer60,1)
  numbars=size(MPchglayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MPchglayer60(:,jj),VarPchglayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Pchglayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Pchglayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarPchglayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarPchglayer60_G%d.eps',NGroup)));

  currfig=find_figure('Pchglayer10Layer_stacked')  
h=bar(MPchglayer10','stacked')

set(gca, 'XTickLabel',GN)

   saveas(currfig,fullfile(savepath,sprintf('BarPchglayer10_G%d_stacked.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarPchglayer10_G%d_stacked.eps',NGroup))); 
   
   

currfig=find_figure(sprintf('Pchglayer10Layer',i))

  h=bar(MPchglayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MPchglayer10,1)
  numbars=size(MPchglayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MPchglayer10(:,jj),VarPchglayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Pchglayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Pchglayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarPchglayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarPchglayer10_G%d.eps',NGroup)));








L=ElementN+1;
for i=1:L
if Ninput==2
        [hh,c]=StAnalysis( Cs{1}.Pchglayer60(:,i),Cs{2}.Pchglayer60(:,i),type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.Pchglayer60(:,i),Cs{2}.Pchglayer60(:,i),Cs{3}.Pchglayer60(:,i))
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.Pchglayer60(:,i),Cs{2}.Pchglayer60(:,i),Cs{3}.Pchglayer60(:,i),Cs{4}.Pchglayer60(:,i))
         elseif Ninput==5
        c =StAnalysis5(Cs{1}.Pchglayer60(:,i),Cs{2}.Pchglayer60(:,i),Cs{3}.Pchglayer60(:,i),Cs{4}.Pchglayer60(:,i),Cs{5}.Pchglayer60(:,i))
end


 currfig=find_figure('PercentageCharg_60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Pchglayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('PercentageCharge60')

saveas(currfig,fullfile(savepath,sprintf('CdfPchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPchg%i.eps',60)));


for i=1:L
if Ninput==2
       [hh,c]=StAnalysis( Cs{1}.Pchglayer10(:,i),Cs{2}.Pchglayer10(:,i),type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.Pchglayer10(:,i),Cs{2}.Pchglayer10(:,i),Cs{3}.Pchglayer10(:,i))
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.Pchglayer10(:,i),Cs{2}.Pchglayer10(:,i),Cs{3}.Pchglayer10(:,i),Cs{4}.Pchglayer10(:,i))
         elseif Ninput==5
        c =StAnalysis5(Cs{1}.Pchglayer10(:,i),Cs{2}.Pchglayer10(:,i),Cs{3}.Pchglayer10(:,i),Cs{4}.Pchglayer10(:,i),Cs{5}.Pchglayer10(:,i))
end


 currfig=find_figure('PercentageCharg_10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Pchglayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('PercentageCharge10')
saveas(currfig,fullfile(savepath,sprintf('CdfPchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPchg%i.eps',10)));





MParealayer60=[];VarParealayer60=[];MParealayer10=[];VarParealayer10=[];
    for jd=1:NGroup
        Parealayer60=Cs{jd}.Parealayer60;Parealayer10=Cs{jd}.Parealayer10;
         mParealayer60=[];vParealayer60=[];mParealayer10=[];vParealayer10=[];
         for jlayer=1:size(Dist8060L,2)
             mParealayer60=[mParealayer60 nanmean(Parealayer60(:,jlayer))];
             vParealayer60=[vParealayer60 nanstd(Parealayer60(:,jlayer))./sqrt(size(Parealayer60,1))];
                  mParealayer10=[mParealayer10 nanmean(Parealayer10(:,jlayer))];
             vParealayer10=[vParealayer10 nanstd(Parealayer10(:,jlayer))./sqrt(size(Parealayer10,1))];
         end
      MParealayer60=[MParealayer60;mParealayer60];
     VarParealayer60=[VarParealayer60;vParealayer60]; 
      MParealayer10=[MParealayer10;mParealayer10];
     VarParealayer10=[VarParealayer10;vParealayer10]; 
      msfile=sprintf('MeanStd_Parealayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mParealayer60','vParealayer60','-ascii') 
    
    msfile=sprintf('MeanStd_Parealayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mParealayer10','vParealayer10','-ascii') 
     
    end
      
  
  MParealayer60=MParealayer60';
  VarParealayer60=VarParealayer60';
  MParealayer10=MParealayer10';
  VarParealayer10=VarParealayer10';
  currfig=find_figure('Parealayer60LayerBoundline')
xa=1:length(MDist(1,:));
[l,p] = boundedline(xa,MParealayer60(1,:), VarParealayer60(1,:), '-g*', xa,MParealayer60(2,:), VarParealayer60(2,:), '-r*',xa,MParealayer60(3,:), VarParealayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Parealayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Parealayer60rBoundline.eps'));
  
   
   currfig=find_figure('Parealayer10LayerBoundline')
xa=1:length(MDist(1,:));
[l,p] = boundedline(xa,MParealayer10(1,:), VarParealayer10(1,:), '-g*', xa,MParealayer10(2,:), VarParealayer10(2,:), '-r*',xa,MParealayer10(3,:), VarParealayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Parealayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Parealayer10rBoundline.eps'));
  

currfig=find_figure(sprintf('Parealayer60Layer',i))

  h=bar(MParealayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MParealayer60,1)
  numbars=size(MParealayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MParealayer60(:,jj),VarParealayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Parealayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Parealayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarParealayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarParealayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('Parealayer10Layer',i))

  h=bar(MParealayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MParealayer10,1)
  numbars=size(MParealayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MParealayer10(:,jj),VarParealayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Parealayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Parealayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarParealayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarParealayer10_G%d.eps',NGroup)));















for i=1:L
if Ninput==2
        [hh,c]=StAnalysis( Cs{1}.Parealayer60(:,i),Cs{2}.Parealayer60(:,i),type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.Parealayer60(:,i),Cs{2}.Parealayer60(:,i),Cs{3}.Parealayer60(:,i))
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.Parealayer60(:,i),Cs{2}.Parealayer60(:,i),Cs{3}.Parealayer60(:,i),Cs{4}.Parealayer60(:,i))
         elseif Ninput==5
        c =StAnalysis5(Cs{1}.Parealayer60(:,i),Cs{2}.Parealayer60(:,i),Cs{3}.Parealayer60(:,i),Cs{4}.Parealayer60(:,i),Cs{5}.Parealayer60(:,i))
end


 currfig=find_figure('PercentageArea_60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Parealayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('PercentageArea60')
saveas(currfig,fullfile(savepath,sprintf('CdfPAREA%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPAREA%i.eps',60)));


for i=1:L
if Ninput==2
        [hh,c]=StAnalysis( Cs{1}.Parealayer10(:,i),Cs{2}.Parealayer10(:,i),type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.Parealayer10(:,i),Cs{2}.Parealayer10(:,i),Cs{3}.Parealayer10(:,i))
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.Parealayer10(:,i),Cs{2}.Parealayer10(:,i),Cs{3}.Parealayer10(:,i),Cs{4}.Parealayer10(:,i))
         elseif Ninput==5
        c =StAnalysis5(Cs{1}.Parealayer10(:,i),Cs{2}.Parealayer10(:,i),Cs{3}.Parealayer10(:,i),Cs{4}.Parealayer10(:,i),Cs{5}.Parealayer10(:,i))
end


 currfig=find_figure('PercentageArea10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Parealayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('PercentageArea10')
saveas(currfig,fullfile(savepath,sprintf('CdfPAREA%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPAREA%i.eps',10)));







MPpeaklayer60=[];VarPpeaklayer60=[];MPpeaklayer10=[];VarPpeaklayer10=[];
    for jd=1:NGroup
        Ppeaklayer60=Cs{jd}.Ppeaklayer60;Ppeaklayer10=Cs{jd}.Ppeaklayer10;
         mPpeaklayer60=[];vPpeaklayer60=[];mPpeaklayer10=[];vPpeaklayer10=[];
         for jlayer=1:size(Dist8060L,2)
             mPpeaklayer60=[mPpeaklayer60 nanmean(Ppeaklayer60(:,jlayer))];
             vPpeaklayer60=[vPpeaklayer60 nanstd(Ppeaklayer60(:,jlayer))./sqrt(size(Ppeaklayer60,1))];
                  mPpeaklayer10=[mPpeaklayer10 nanmean(Ppeaklayer10(:,jlayer))];
             vPpeaklayer10=[vPpeaklayer10 nanstd(Ppeaklayer10(:,jlayer))./sqrt(size(Ppeaklayer10,1))];
         end
      MPpeaklayer60=[MPpeaklayer60;mPpeaklayer60];
     VarPpeaklayer60=[VarPpeaklayer60;vPpeaklayer60]; 
      MPpeaklayer10=[MPpeaklayer10;mPpeaklayer10];
     VarPpeaklayer10=[VarPpeaklayer10;vPpeaklayer10]; 
      msfile=sprintf('MeanStd_Ppeaklayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mPpeaklayer60','vPpeaklayer60','-ascii') 
    
    msfile=sprintf('MeanStd_Ppeaklayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mPpeaklayer10','vPpeaklayer10','-ascii') 
     
    end
      
  
  MPpeaklayer60=MPpeaklayer60';
  VarPpeaklayer60=VarPpeaklayer60';
  MPpeaklayer10=MPpeaklayer10';
  VarPpeaklayer10=VarPpeaklayer10';
  currfig=find_figure('Ppeaklayer60LayerBoundline')
xa=1:length( MPpeaklayer60(1,:));
[l,p] = boundedline(xa,MPpeaklayer60(1,:), VarPpeaklayer60(1,:), '-g*', xa,MPpeaklayer60(2,:), VarPpeaklayer60(2,:), '-r*',xa,MPpeaklayer60(3,:), VarPpeaklayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Ppeaklayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Ppeaklayer60rBoundline.eps'));
  
  currfig=find_figure('Ppeaklayer10LayerBoundline')
xa=1:length(MPpeaklayer10(1,:));
[l,p] = boundedline(xa,MPpeaklayer10(1,:), VarPpeaklayer10(1,:), '-g*', xa,MPpeaklayer10(2,:), VarPpeaklayer10(2,:), '-r*',xa,MPpeaklayer10(3,:), VarPpeaklayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Ppeaklayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Ppeaklayer10rBoundline.eps'));
  
  

currfig=find_figure(sprintf('Ppeaklayer60Layer',i))

  h=bar(MPpeaklayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MPpeaklayer60,1)
  numbars=size(MPpeaklayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MPpeaklayer60(:,jj),VarPpeaklayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Ppeaklayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Ppeaklayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarPpeaklayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarPpeaklayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('Ppeaklayer10Layer',i))

  h=bar(MPpeaklayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MPpeaklayer10,1)
  numbars=size(MPpeaklayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MPpeaklayer10(:,jj),VarPpeaklayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Ppeaklayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Ppeaklayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarPpeaklayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarPpeaklayer10_G%d.eps',NGroup)));












for i=1:L
if Ninput==2
        [hh,c]=StAnalysis( Cs{1}.Ppeaklayer60(:,i),Cs{2}.Ppeaklayer60(:,i),type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.Ppeaklayer60(:,i),Cs{2}.Ppeaklayer60(:,i),Cs{3}.Ppeaklayer60(:,i))
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.Ppeaklayer60(:,i),Cs{2}.Ppeaklayer60(:,i),Cs{3}.Ppeaklayer60(:,i),Cs{4}.Ppeaklayer60(:,i))
         elseif Ninput==5
        c =StAnalysis5(Cs{1}.Ppeaklayer60(:,i),Cs{2}.Ppeaklayer60(:,i),Cs{3}.Ppeaklayer60(:,i),Cs{4}.Ppeaklayer60(:,i),Cs{5}.Ppeaklayer60(:,i))
end


 currfig=find_figure('PercentagePeak60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Ppeaklayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('PercentagePeak60')
saveas(currfig,fullfile(savepath,sprintf('CdfPpk%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPpk%i.eps',60)));

for i=1:L
if Ninput==2
        [hh,c]=StAnalysis( Cs{1}.Ppeaklayer10(:,i),Cs{2}.Ppeaklayer10(:,i),type)
    elseif Ninput==3
        c=StAnalysis3(Cs{1}.Ppeaklayer10(:,i),Cs{2}.Ppeaklayer10(:,i),Cs{3}.Ppeaklayer10(:,i))
    elseif Ninput==4
        c =StAnalysis4(Cs{1}.Ppeaklayer10(:,i),Cs{2}.Ppeaklayer10(:,i),Cs{3}.Ppeaklayer10(:,i),Cs{4}.Ppeaklayer10(:,i))
         elseif Ninput==5
        c =StAnalysis5(Cs{1}.Ppeaklayer10(:,i),Cs{2}.Ppeaklayer10(:,i),Cs{3}.Ppeaklayer10(:,i),Cs{4}.Ppeaklayer10(:,i),Cs{5}.Ppeaklayer10(:,i))
end


 currfig=find_figure('PercentagePeak10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Ppeaklayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('PercentagePeak10')
saveas(currfig,fullfile(savepath,sprintf('CdfPpk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfPpk%i.eps',10)));
%----------------------------------------------------------------------------------






Mchglayer60=[];Varchglayer60=[];Mchglayer10=[];Varchglayer10=[];
    for jd=1:NGroup
        chglayer60=Cs{jd}.chglayer60;chglayer10=Cs{jd}.chglayer10;
         mchglayer60=[];vchglayer60=[];mchglayer10=[];vchglayer10=[];
         for jlayer=1:size(Dist8060L,2)
             mchglayer60=[mchglayer60 nanmean(chglayer60(:,jlayer))];
             vchglayer60=[vchglayer60 nanstd(chglayer60(:,jlayer))./sqrt(size(chglayer60,1))];
                  mchglayer10=[mchglayer10 nanmean(chglayer10(:,jlayer))];
             vchglayer10=[vchglayer10 nanstd(chglayer10(:,jlayer))./sqrt(size(chglayer10,1))];
         end
      Mchglayer60=[Mchglayer60;mchglayer60];
     Varchglayer60=[Varchglayer60;vchglayer60]; 
      Mchglayer10=[Mchglayer10;mchglayer10];
     Varchglayer10=[Varchglayer10;vchglayer10]; 
      msfile=sprintf('MeanStd_chglayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mchglayer60','vchglayer60','-ascii') 
    
    msfile=sprintf('MeanStd_chglayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mchglayer10','vchglayer10','-ascii') 
     
    end
      
  
  Mchglayer60=Mchglayer60';
  Varchglayer60=Varchglayer60';
  Mchglayer10=Mchglayer10';
  Varchglayer10=Varchglayer10';
  
  
  currfig=find_figure('chargelayer60LayerBoundline')
xa=1:length( Mchglayer60(1,:));
[l,p] = boundedline(xa,Mchglayer60(1,:), Varchglayer60(1,:), '-g*', xa,Mchglayer60(2,:), Varchglayer60(2,:), '-r*',xa,Mchglayer60(3,:), Varchglayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'chargelayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'chargelayer60rBoundline.eps'));
  
  currfig=find_figure('chargelayer10LayerBoundline')
xa=1:length(Mchglayer10(1,:));
[l,p] = boundedline(xa,Mchglayer10(1,:), Varchglayer10(1,:), '-g*', xa,Mchglayer10(2,:), Varchglayer10(2,:), '-r*',xa,Mchglayer10(3,:), Varchglayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'chargelayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'chargelayer10rBoundline.eps'));
  

currfig=find_figure(sprintf('chglayer60Layer',i))

  h=bar(Mchglayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mchglayer60,1)
  numbars=size(Mchglayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mchglayer60(:,jj),Varchglayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          
          Ds{ii}=Cs{ii}.chglayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_chglayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barchglayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barchglayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('chglayer10Layer',i))

  h=bar(Mchglayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mchglayer10,1)
  numbars=size(Mchglayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mchglayer10(:,jj),Varchglayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
        Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.chglayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_chglayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barchglayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barchglayer10_G%d.eps',NGroup)));






for i=1:L
    
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.chglayer60(:,i);
      end
     
      if NGroup>2
          c=StAnalysisMulti(Ds);
      else
           [hh,c]=StAnalysis( Cs{1}.chglayer60(:,i),Cs{2}.chglayer60(:,i),type);
      end


 currfig=find_figure('Charge60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.chglayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Charge60')
saveas(currfig,fullfile(savepath,sprintf('Cdfchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfchg%i.eps',60)));

%------------------------------total charge 60----------------------------------


if Ninput==2
        [hh,c]=StAnalysis(sum(Cs{1}.chglayer60,2),sum(Cs{2}.chglayer60,2),type)
    elseif Ninput==3
        c=StAnalysis3(sum(Cs{1}.chglayer60,2),sum(Cs{2}.chglayer60,2),sum(Cs{3}.chglayer60,2))
    elseif Ninput==4
        c =StAnalysis4(sum(Cs{1}.chglayer60,2),sum(Cs{2}.chglayer60,2),sum(Cs{3}.chglayer60,2),sum(Cs{4}.chglayer60,2))
         elseif Ninput==5
        c =StAnalysis5(sum(Cs{1}.chglayer60,2),sum(Cs{2}.chglayer60,2),sum(Cs{3}.chglayer60,2),sum(Cs{4}.chglayer60,2),sum(Cs{5}.chglayer60,2))
end


 currfig=find_figure('TotalCharge60');
subplot(1,2,1)

hold on
for j=1:Ninput
h1=cdfplot(sum(Cs{j}.chglayer60,2));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off

 title('TotalCharge60')
saveas(currfig,fullfile(savepath,sprintf('CdfTotalchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfTotalchg%i.eps',60)));


%--------------------------------------------------------------------------



for i=1:L
         Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.chglayer10(:,i);
      end
     
      if NGroup>2
          c=StAnalysisMulti(Ds);
      else
          [hh,c]=StAnalysis( Cs{1}.chglayer10(:,i),Cs{2}.chglayer10(:,i),type);
      end



 currfig=find_figure('Charge10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.chglayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Charge10')
saveas(currfig,fullfile(savepath,sprintf('Cdfchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfchg%i.eps',10)));


%------------------------------total charge 10----------------------------------


if Ninput==2
        [hh,c]=StAnalysis(sum(Cs{1}.chglayer10,2),sum(Cs{2}.chglayer10,2),type)
    elseif Ninput==3
        c=StAnalysis3(sum(Cs{1}.chglayer10,2),sum(Cs{2}.chglayer10,2),sum(Cs{3}.chglayer10,2))
    elseif Ninput==4
        c =StAnalysis4(sum(Cs{1}.chglayer10,2),sum(Cs{2}.chglayer10,2),sum(Cs{3}.chglayer10,2),sum(Cs{4}.chglayer10,2))
         elseif Ninput==5
        c =StAnalysis5(sum(Cs{1}.chglayer10,2),sum(Cs{2}.chglayer10,2),sum(Cs{3}.chglayer10,2),sum(Cs{4}.chglayer10,2),sum(Cs{5}.chglayer10,2))
end


 currfig=find_figure('TotalCharge10');

subplot(1,2,1)
hold on
for j=1:Ninput
h1=cdfplot(sum(Cs{j}.chglayer10,2));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off

 title('TotalCharge10')
saveas(currfig,fullfile(savepath,sprintf('CdfTotalchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfTotalchg%i.eps',10)));


%--------------------------------------------------------------------------


Mpeaklayer60=[];Varpeaklayer60=[];Mpeaklayer10=[];Varpeaklayer10=[];
    for jd=1:NGroup
        peaklayer60=Cs{jd}.peaklayer60;peaklayer10=Cs{jd}.peaklayer10;
         mpeaklayer60=[];vpeaklayer60=[];mpeaklayer10=[];vpeaklayer10=[];
         for jlayer=1:size(Dist8060L,2)
             mpeaklayer60=[mpeaklayer60 nanmean(peaklayer60(:,jlayer))];
             vpeaklayer60=[vpeaklayer60 nanstd(peaklayer60(:,jlayer))./length(peaklayer60(:,jlayer))];
                  mpeaklayer10=[mpeaklayer10 nanmean(peaklayer10(:,jlayer))];
             vpeaklayer10=[vpeaklayer10 nanstd(peaklayer10(:,jlayer))./length(peaklayer60(:,jlayer))];
         end
      Mpeaklayer60=[Mpeaklayer60;mpeaklayer60];
     Varpeaklayer60=[Varpeaklayer60;vpeaklayer60]; 
      Mpeaklayer10=[Mpeaklayer10;mpeaklayer10];
     Varpeaklayer10=[Varpeaklayer10;vpeaklayer10]; 
      msfile=sprintf('MeanStd_peaklayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mpeaklayer60','vpeaklayer60','-ascii') 
    
    msfile=sprintf('MeanStd_peaklayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mpeaklayer10','vpeaklayer10','-ascii') 
     
    end
      
  
  Mpeaklayer60=Mpeaklayer60';
  Varpeaklayer60=Varpeaklayer60';
  Mpeaklayer10=Mpeaklayer10';
  Varpeaklayer10=Varpeaklayer10';
  currfig=find_figure('Peaklayer60LayerBoundline')
xa=1:length( Mpeaklayer60(1,:));
[l,p] = boundedline(xa,Mpeaklayer60(1,:), Varpeaklayer60(1,:), '-g*', xa,Mpeaklayer60(2,:), Varpeaklayer60(2,:), '-r*',xa,Mpeaklayer60(3,:), Varpeaklayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Peaklayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Peaklayer60rBoundline.eps'));
  
  currfig=find_figure('Peaklayer10LayerBoundline')
xa=1:length(Mpeaklayer10(1,:));
[l,p] = boundedline(xa,Mpeaklayer10(1,:), Varpeaklayer10(1,:), '-g*', xa,Mpeaklayer10(2,:), Varpeaklayer10(2,:), '-r*',xa,Mpeaklayer10(3,:), Varpeaklayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Peaklayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Peaklayer10rBoundline.eps'));
  
  

currfig=find_figure(sprintf('peaklayer60Layer',i))

  h=bar(Mpeaklayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mpeaklayer60,1)
  numbars=size(Mpeaklayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mpeaklayer60(:,jj),Varpeaklayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.peaklayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_peaklayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barpeaklayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barpeaklayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('peaklayer10Layer',i))

  h=bar(Mpeaklayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mpeaklayer10,1)
  numbars=size(Mpeaklayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mpeaklayer10(:,jj),Varpeaklayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.peaklayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_peaklayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barpeaklayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barpeaklayer10_G%d.eps',NGroup)));













for i=1:L
    
         Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.peaklayer60(:,i);
      end
     
      if NGroup>2
          c=StAnalysisMulti(Ds);
      else
          [hh,c]=StAnalysis( Cs{1}.peaklayer60(:,i),Cs{2}.peaklayer60(:,i),type);
      end


 currfig=find_figure('Peak60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.peaklayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Peak60')
saveas(currfig,fullfile(savepath,sprintf('Cdfpk%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfpk%i.eps',60)));

%------------------------------total Peak 60----------------------------------


if Ninput==2
        [hh,c]=StAnalysis(sum(Cs{1}.peaklayer60,2),sum(Cs{2}.peaklayer60,2),type)
    elseif Ninput==3
        c=StAnalysis3(sum(Cs{1}.peaklayer60,2),sum(Cs{2}.peaklayer60,2),sum(Cs{3}.peaklayer60,2))
    elseif Ninput==4
        c =StAnalysis4(sum(Cs{1}.peaklayer60,2),sum(Cs{2}.peaklayer60,2),sum(Cs{3}.peaklayer60,2),sum(Cs{4}.peaklayer60,2))
         elseif Ninput==5
        c =StAnalysis5(sum(Cs{1}.peaklayer60,2),sum(Cs{2}.peaklayer60,2),sum(Cs{3}.peaklayer60,2),sum(Cs{4}.peaklayer60,2),sum(Cs{5}.peaklayer60,2))
end


 currfig=find_figure('TotalPeak60');

subplot(1,2,1)
hold on
for j=1:Ninput
h1=cdfplot(sum(Cs{j}.peaklayer60,2));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off

 title('TotalPeak60')
saveas(currfig,fullfile(savepath,sprintf('CdfTotalPk%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfTotalPk%i.eps',60)));


%--------------------------------------------------------------------------




for i=1:L
           Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.peaklayer10(:,i);
      end
     
      if NGroup>2
          c=StAnalysisMulti(Ds);
      else
          [hh,c]=StAnalysis( Cs{1}.peaklayer10(:,i),Cs{2}.peaklayer10(:,i),type);
      end


 currfig=find_figure('Peak10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.peaklayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Peak10')
saveas(currfig,fullfile(savepath,sprintf('Cdfpk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('Cdfpk%i.eps',10)));
%-------------------------------------------------------------------------

%------------------------------total Peak 10----------------------------------


if Ninput==2
        [hh,c]=StAnalysis(sum(Cs{1}.peaklayer10,2),sum(Cs{2}.peaklayer10,2),type)
    elseif Ninput==3
        c=StAnalysis3(sum(Cs{1}.peaklayer10,2),sum(Cs{2}.peaklayer10,2),sum(Cs{3}.peaklayer10,2))
    elseif Ninput==4
        c =StAnalysis4(sum(Cs{1}.peaklayer10,2),sum(Cs{2}.peaklayer10,2),sum(Cs{3}.peaklayer10,2),sum(Cs{4}.peaklayer10,2))
         elseif Ninput==5
        c =StAnalysis5(sum(Cs{1}.peaklayer10,2),sum(Cs{2}.peaklayer10,2),sum(Cs{3}.peaklayer10,2),sum(Cs{4}.peaklayer10,2),sum(Cs{5}.peaklayer10,2))
end


 currfig=find_figure('TotalPeak10');

subplot(1,2,1)
hold on
for j=1:Ninput
h1=cdfplot(sum(Cs{j}.peaklayer10,2));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off

 title('TotalPeak10')
saveas(currfig,fullfile(savepath,sprintf('CdfTotalPk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfTotalPk%i.eps',10)));


%--------------------------------------------------------------------------




MMchglayer60=[];VarMchglayer60=[];MMchglayer10=[];VarMchglayer10=[];
    for jd=1:NGroup
        Mchglayer60=Cs{jd}.Mchglayer60;
        Mchglayer10=Cs{jd}.Mchglayer10;
        Mchglayer60(isnan(Mchglayer60))=0;
            Mchglayer10(isnan(Mchglayer10))=0;
         mMchglayer60=[];vMchglayer60=[];mMchglayer10=[];vMchglayer10=[];
         for jlayer=1:size(Dist8060L,2)
             
             mMchglayer60=[mMchglayer60 nanmean(Mchglayer60(:,jlayer))];
             vMchglayer60=[vMchglayer60 nanstd(Mchglayer60(:,jlayer))./sqrt(size(Mchglayer60,1))];
                  mMchglayer10=[mMchglayer10 nanmean(Mchglayer10(:,jlayer))];
             vMchglayer10=[vMchglayer10 nanstd(Mchglayer10(:,jlayer))./sqrt(size(Mchglayer10,1))];
         end
      MMchglayer60=[MMchglayer60;mMchglayer60];
     VarMchglayer60=[VarMchglayer60;vMchglayer60]; 
      MMchglayer10=[MMchglayer10;mMchglayer10];
     VarMchglayer10=[VarMchglayer10;vMchglayer10]; 
      msfile=sprintf('MeanStd_Meanchglayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mMchglayer60','vMchglayer60','-ascii') 
    
    msfile=sprintf('MeanStd_Meanchglayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mMchglayer10','vMchglayer10','-ascii') 
     
    end
      
  
  MMchglayer60=MMchglayer60';
  VarMchglayer60=VarMchglayer60';
  MMchglayer10=MMchglayer10';
  VarMchglayer10=VarMchglayer10';
  currfig=find_figure('MeanChargelayer60LayerBoundline')
xa=1:length(MMchglayer60(1,:));
[l,p] = boundedline(xa,MMchglayer60(1,:), VarMchglayer60(1,:), '-g*', xa,MMchglayer60(2,:), VarMchglayer60(2,:), '-r*',xa,MMchglayer60(3,:), VarMchglayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MeanChargelayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MeanChargelayer60rBoundline.eps'));
  
  currfig=find_figure('MeanChargelayer10LayerBoundline')
xa=1:length(MMchglayer10(1,:));
[l,p] = boundedline(xa,MMchglayer10(1,:), VarMchglayer10(1,:), '-g*', xa,MMchglayer10(2,:), VarMchglayer10(2,:), '-r*',xa,MMchglayer10(3,:), VarMchglayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MeanChargelayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MeanChargelayer10rBoundline.eps'));
  
  

currfig=find_figure(sprintf('Meanchglayer60Layer',i))

  h=bar(MMchglayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MMchglayer60,1)
  numbars=size(MMchglayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MMchglayer60(:,jj),VarMchglayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Mchglayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Meanchglayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarMeanchglayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarMeanchglayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('Mchglayer10Layer',i))

  h=bar(MMchglayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MMchglayer10,1)
  numbars=size(MMchglayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MMchglayer10(:,jj),VarMchglayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      noNaNLength=[];
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Mchglayer10(:,jj);
          noNaNLength=[noNaNLength;length(find(~isnan(Ds{ii})))];
      end
      
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if ~length(find(noNaNLength<4))
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
      else P=NaN; p=NaN;
      end
     pfile=sprintf('statisticP_Meanchglayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarMeanchglayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarMeanchglayer10_G%d.eps',NGroup)));











for i=1:L
           Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Mchglayer60(:,i);
      end
     
      if NGroup>2
          c=StAnalysisMulti(Ds);
      else
       [hh,c]=StAnalysis( Cs{1}.Mchglayer60(:,i),Cs{2}.Mchglayer60(:,i),type);
      end



 currfig=find_figure('MeanCharge60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Mchglayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MeanCharge60')
saveas(currfig,fullfile(savepath,sprintf('CdfMeanchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMeanchg%i.eps',60)));
for i=1:L
    Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.Mchglayer10(:,i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds)
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end



 currfig=find_figure('MeanCharge10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Mchglayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MeanCharge10')
saveas(currfig,fullfile(savepath,sprintf('CdfMeanchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMeanchg%i.eps',10)));


MMpeaklayer60=[];VarMpeaklayer60=[];MMpeaklayer10=[];VarMpeaklayer10=[];
    for jd=1:NGroup
        Mpeaklayer60=Cs{jd}.Mpeaklayer60;Mpeaklayer10=Cs{jd}.Mpeaklayer10;
            Mpeaklayer60(isnan(Mpeaklayer60))=0;
             Mpeaklayer10(isnan(Mpeaklayer10))=0;
         mMpeaklayer60=[];vMpeaklayer60=[];mMpeaklayer10=[];vMpeaklayer10=[];
         for jlayer=1:size(Dist8060L,2)
             mMpeaklayer60=[mMpeaklayer60 nanmean(Mpeaklayer60(:,jlayer))];
             vMpeaklayer60=[vMpeaklayer60 nanstd(Mpeaklayer60(:,jlayer))./sqrt(size(Mpeaklayer60,1))];
                  mMpeaklayer10=[mMpeaklayer10 nanmean(Mpeaklayer10(:,jlayer))];
             vMpeaklayer10=[vMpeaklayer10 nanstd(Mpeaklayer10(:,jlayer))./sqrt(size(Mpeaklayer10,1))];
         end
      MMpeaklayer60=[MMpeaklayer60;mMpeaklayer60];
     VarMpeaklayer60=[VarMpeaklayer60;vMpeaklayer60]; 
      MMpeaklayer10=[MMpeaklayer10;mMpeaklayer10];
     VarMpeaklayer10=[VarMpeaklayer10;vMpeaklayer10]; 
      msfile=sprintf('MeanStd_Meanpeaklayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mMpeaklayer60','vMpeaklayer60','-ascii') 
    
    msfile=sprintf('MeanStd_Meanpeaklayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mMpeaklayer10','vMpeaklayer10','-ascii') 
     
    end
      
  
  MMpeaklayer60=MMpeaklayer60';
  VarMpeaklayer60=VarMpeaklayer60';
  MMpeaklayer10=MMpeaklayer10';
  VarMpeaklayer10=VarMpeaklayer10';
  
   currfig=find_figure('MeanPeaklayer60LayerBoundline')
xa=1:length(MMpeaklayer60(1,:));
[l,p] = boundedline(xa,MMpeaklayer60(1,:), VarMpeaklayer60(1,:), '-g*', xa,MMpeaklayer60(2,:), VarMpeaklayer60(2,:), '-r*',xa,MMpeaklayer60(3,:), VarMpeaklayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MeanPeaklayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MeanPeaklayer60rBoundline.eps'));
  
  currfig=find_figure('MeanPeaklayer10LayerBoundline')
xa=1:length(MMpeaklayer10(1,:));
[l,p] = boundedline(xa,MMpeaklayer10(1,:), VarMpeaklayer10(1,:), '-g*', xa,MMpeaklayer10(2,:), VarMpeaklayer10(2,:), '-r*',xa,MMpeaklayer10(3,:), VarMpeaklayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MeanPeaklayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MeanPeaklayer10rBoundline.eps'));
  

currfig=find_figure(sprintf('Mpeaklayer60Layer',i))

  h=bar(MMpeaklayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MMpeaklayer60,1)
  numbars=size(MMpeaklayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MMpeaklayer60(:,jj),VarMpeaklayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Mpeaklayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Meanpeaklayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarMeanpeaklayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarMeanpeaklayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('Mpeaklayer10Layer',i))

  h=bar(MMpeaklayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MMpeaklayer10,1)
  numbars=size(MMpeaklayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MMpeaklayer10(:,jj),VarMpeaklayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Mpeaklayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Meanpeaklayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarMeanpeaklayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarMeanpeaklayer10_G%d.eps',NGroup)));
















for i=1:L
  Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.Mpeaklayer60(:,i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   
    


 currfig=find_figure('MeanPeak60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Mpeaklayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MeanPeak60')
saveas(currfig,fullfile(savepath,sprintf('CdfMeanpk%i.fig',60)),'fig'); 
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMeanpk%i.eps',60)));
for i=1:L
   
     Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.Mpeaklayer10(:,jj);
      end
     
      if NGroup>2
          c=StAnalysisMulti(Ds);
      else
          [hh,c]=StAnalysis( Cs{1}.Mpeaklayer10(:,i),Cs{2}.Mpeaklayer10(:,i),type)
      end



 currfig=find_figure('Mpeaklayer10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.Mpeaklayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Mpeaklayer10')
saveas(currfig,fullfile(savepath,sprintf('CdfMeanpk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMeanpk%i.eps',10)));

%-------------------------------------------------------------------------


Marealayer60=[];Vararealayer60=[];Marealayer10=[];Vararealayer10=[];
    for jd=1:NGroup
        arealayer60=Cs{jd}.arealayer60;arealayer10=Cs{jd}.arealayer10;
         marealayer60=[];varealayer60=[];marealayer10=[];varealayer10=[];
         for jlayer=1:size(Dist8060L,2)
             marealayer60=[marealayer60 nanmean(arealayer60(:,jlayer))];
             varealayer60=[varealayer60 nanstd(arealayer60(:,jlayer))./sqrt(size(arealayer60,1))];
                  marealayer10=[marealayer10 nanmean(arealayer10(:,jlayer))];
             varealayer10=[varealayer10 nanstd(arealayer10(:,jlayer))./sqrt(size(arealayer10,1))];
         end
      Marealayer60=[Marealayer60;marealayer60];
     Vararealayer60=[Vararealayer60;varealayer60]; 
      Marealayer10=[Marealayer10;marealayer10];
     Vararealayer10=[Vararealayer10;varealayer10]; 
      msfile=sprintf('MeanStd_Arealayer60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'marealayer60','varealayer60','-ascii') 
    
    msfile=sprintf('MeanStd_Arealayer10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'marealayer10','varealayer10','-ascii') 
     
    end
      
  
  Marealayer60=Marealayer60';
  Vararealayer60=Vararealayer60';
  Marealayer10=Marealayer10';
  Vararealayer10=Vararealayer10';
  
  currfig=find_figure('Arealayer60LayerBoundline')
xa=1:length(MMpeaklayer60(1,:));
[l,p] = boundedline(xa,Marealayer60(1,:), Vararealayer60(1,:), '-g*', xa,Marealayer60(2,:), Vararealayer60(2,:), '-r*',xa,Marealayer60(3,:), Vararealayer60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Arealayer60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Arealayer60rBoundline.eps'));
  
  currfig=find_figure('Arealayer10LayerBoundline')
xa=1:length(Marealayer10(1,:));
[l,p] = boundedline(xa,Marealayer10(1,:), Vararealayer10(1,:), '-g*', xa,Marealayer10(2,:), Vararealayer10(2,:), '-r*',xa,Marealayer10(3,:), Vararealayer10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'Arealayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'Arealayer10rBoundline.eps'));
  

currfig=find_figure(sprintf('arealayer60Layer',i))

  h=bar(Marealayer60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Marealayer60,1)
  numbars=size(Marealayer60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Marealayer60(:,jj),Vararealayer60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.arealayer60(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Arealayer60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Bararealayer60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Bararealayer60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('arealayer10Layer',i))

  h=bar(Marealayer10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Marealayer10,1)
  numbars=size(Marealayer10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Marealayer10(:,jj),Vararealayer10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MDist,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.arealayer10(:,jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_Arealayer10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarArealayer10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarArealayer10_G%d.eps',NGroup)));














for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.arealayer60(:,i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('Area60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.arealayer60(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Area60')
saveas(currfig,fullfile(savepath,sprintf('CdfArea%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfArea%i.eps',60)));

%------------------------------total Area 60----------------------------------


















if Ninput==2
        [hh,c]=StAnalysis(sum(Cs{1}.arealayer60,2),sum(Cs{2}.arealayer60,2),type)
    elseif Ninput==3
        c=StAnalysis3(sum(Cs{1}.arealayer60,2),sum(Cs{2}.arealayer60,2),sum(Cs{3}.arealayer60,2))
    elseif Ninput==4
        c =StAnalysis4(sum(Cs{1}.arealayer60,2),sum(Cs{2}.arealayer60,2),sum(Cs{3}.arealayer60,2),sum(Cs{4}.arealayer60,2))
         elseif Ninput==5
        c =StAnalysis5(sum(Cs{1}.arealayer60,2),sum(Cs{2}.arealayer60,2),sum(Cs{3}.arealayer60,2),sum(Cs{4}.arealayer60,2),sum(Cs{5}.arealayer60,2))
end


 currfig=find_figure('TotalArea60');
subplot(1,2,1)

hold on
for j=1:Ninput
h1=cdfplot(sum(Cs{j}.arealayer60,2));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off

 title('TotalArea60')
saveas(currfig,fullfile(savepath,sprintf('CdfTotalArea%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfTotalArea%i.eps',60)));


%--------------------------------------------------------------------------


for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.arealayer10(:,i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('Area10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.arealayer10(:,i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Area10')
saveas(currfig,fullfile(savepath,sprintf('CdfArea%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfArea%i.eps',10)));

%___________________________________________________________________________

%------------------------------total Area 10----------------------------------


if Ninput==2
        [hh,c]=StAnalysis(sum(Cs{1}.arealayer10,2),sum(Cs{2}.arealayer10,2),type)
    elseif Ninput==3
        c=StAnalysis3(sum(Cs{1}.arealayer10,2),sum(Cs{2}.arealayer10,2),sum(Cs{3}.arealayer10,2))
    elseif Ninput==4
        c =StAnalysis4(sum(Cs{1}.arealayer10,2),sum(Cs{2}.arealayer10,2),sum(Cs{3}.arealayer10,2),sum(Cs{4}.arealayer10,2))
         elseif Ninput==5
        c =StAnalysis5(sum(Cs{1}.arealayer10,2),sum(Cs{2}.arealayer10,2),sum(Cs{3}.arealayer10,2),sum(Cs{4}.arealayer10,2),sum(Cs{5}.arealayer10,2))
end


 currfig=find_figure('TotalArea10');

subplot(1,2,1)
hold on
for j=1:Ninput
h1=cdfplot(sum(Cs{j}.arealayer10,2));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off

 title('TotalArea10')
saveas(currfig,fullfile(savepath,sprintf('CdfTotalArea%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfTotalArea%i.eps',10)));


%--------------------------------------------------------------------------


Mmarginalchg60=[];Varmarginalchg60=[];Mmarginalchg10=[];Varmarginalchg10=[];
    for jd=1:NGroup
        marginalchg60=Cs{jd}.D60chg(:,end-ElementN:end);marginalchg10=Cs{jd}.D10chg(:,end-ElementN:end);
        marginalchg60(isnan(marginalchg60))=0;
         marginalchg10(isnan(marginalchg10))=0;
         mmarginalchg60=[];vmarginalchg60=[];mmarginalchg10=[];vmarginalchg10=[];
         
         for jlayer=1:size(Dist8060L,2)
             mmarginalchg60=[mmarginalchg60 nanmean(marginalchg60(:,jlayer))];
             vmarginalchg60=[vmarginalchg60 nanstd(marginalchg60(:,jlayer))./sqrt(size(marginalchg60,1))];
                  mmarginalchg10=[mmarginalchg10 nanmean(marginalchg10(:,jlayer))];
             vmarginalchg10=[vmarginalchg10 nanstd(marginalchg10(:,jlayer))./sqrt(size(marginalchg10,1))];
         end
      Mmarginalchg60=[Mmarginalchg60;mmarginalchg60];
     Varmarginalchg60=[Varmarginalchg60;vmarginalchg60]; 
      Mmarginalchg10=[Mmarginalchg10;mmarginalchg10];
     Varmarginalchg10=[Varmarginalchg10;vmarginalchg10]; 
      msfile=sprintf('MeanStd_marginalchg60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mmarginalchg60','vmarginalchg60','-ascii') 
    
    msfile=sprintf('MeanStd_marginalchg10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mmarginalchg10','vmarginalchg10','-ascii') 
     
    end
      
  
  Mmarginalchg60=Mmarginalchg60';
  Varmarginalchg60=Varmarginalchg60';
  Mmarginalchg10=Mmarginalchg10';
  Varmarginalchg10=Varmarginalchg10';
  currfig=find_figure('MarginalChargelayer60LayerBoundline')
xa=1:length(MMpeaklayer60(1,:));
[l,p] = boundedline(xa,Mmarginalchg60(1,:), Varmarginalchg60(1,:), '-g*', xa,Mmarginalchg60(2,:), Varmarginalchg60(2,:), '-r*',xa,Mmarginalchg60(3,:), Varmarginalchg60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MarginalCharge60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MarginalCharge60LayerBoundline.eps'));
  
  currfig=find_figure('MarginalChargelayer10LayerBoundline')
xa=1:length(Mmarginalchg10(1,:));
[l,p] = boundedline(xa,Mmarginalchg10(1,:), Varmarginalchg10(1,:), '-g*', xa,Mmarginalchg10(2,:), Varmarginalchg10(2,:), '-r*',xa,Mmarginalchg10(3,:), Varmarginalchg10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MarginalChargelayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MarginalChargelayer10rBoundline.eps'));

currfig=find_figure(sprintf('marginalchg60',i))

  h=bar(Mmarginalchg60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mmarginalchg60,1)
  numbars=size(Mmarginalchg60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mmarginalchg60(:,jj),Varmarginalchg60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(Mmarginalchg60,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.D60chg(:,end-ElementN+jj-1);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_marginalchg60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barmarginalchg60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barmarginalchg60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('marginalchg10',i))

  h=bar(Mmarginalchg10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mmarginalchg10,1)
  numbars=size(Mmarginalchg10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mmarginalchg10(:,jj),Varmarginalchg10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(Mmarginalchg10,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.D10chg(:,end-ElementN+jj-1);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_marginalchg10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barmarginalchg10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barmarginalchg10_G%d.eps',NGroup)));







for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.D60chg(:,end-ElementN+i-1);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('Marginalchg60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.D60chg(:,end-ElementN+i-1));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Marginalchg60')
saveas(currfig,fullfile(savepath,sprintf('CdfMarginalchg%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMarginalchg%i.eps',60)));

%___________________________________________________________________________
for i=1:L
    
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.D10chg(:,end-ElementN+i-1);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('Marginalchg10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.D10chg(:,end-ElementN+i-1));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('Marginalchg10')
saveas(currfig,fullfile(savepath,sprintf('CdfMarginalchg%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMarginalchg%i.eps',10)));

%___________________________________________________________________________


MmarginalPk60=[];VarmarginalPk60=[];MmarginalPk10=[];VarmarginalPk10=[];
    for jd=1:NGroup
        marginalPk60=Cs{jd}.D60Pk(:,end-ElementN:end);marginalPk10=Cs{jd}.D10Pk(:,end-ElementN:end);
        marginalPk60(isnan(marginalPk60))=0;
         marginalPk10(isnan(marginalPk10))=0;
         mmarginalPk60=[];vmarginalPk60=[];mmarginalPk10=[];vmarginalPk10=[];
         for jlayer=1:size(Dist8060L,2)
             mmarginalPk60=[mmarginalPk60 nanmean(marginalPk60(:,jlayer))];
             vmarginalPk60=[vmarginalPk60 nanstd(marginalPk60(:,jlayer))./sqrt(size(marginalPk60,1))];
                  mmarginalPk10=[mmarginalPk10 nanmean(marginalPk10(:,jlayer))];
             vmarginalPk10=[vmarginalPk10 nanstd(marginalPk10(:,jlayer))./sqrt(size(marginalPk10,1))];
         end
      MmarginalPk60=[MmarginalPk60;mmarginalPk60];
     VarmarginalPk60=[VarmarginalPk60;vmarginalPk60]; 
      MmarginalPk10=[MmarginalPk10;mmarginalPk10];
     VarmarginalPk10=[VarmarginalPk10;vmarginalPk10]; 
      msfile=sprintf('MeanStd_marginalPk60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mmarginalPk60','vmarginalPk60','-ascii') 
    
    msfile=sprintf('MeanStd_marginalPk10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mmarginalPk10','vmarginalPk10','-ascii') 
     
    end
      
  
  MmarginalPk60=MmarginalPk60';
  VarmarginalPk60=VarmarginalPk60';
  MmarginalPk10=MmarginalPk10';
  VarmarginalPk10=VarmarginalPk10';
  
  currfig=find_figure('MarginalPeaklayer60LayerBoundline')
xa=1:length(MMpeaklayer60(1,:));
[l,p] = boundedline(xa,MmarginalPk60(1,:), VarmarginalPk60(1,:), '-g*', xa,MmarginalPk60(2,:), VarmarginalPk60(2,:), '-r*',xa,MmarginalPk60(3,:), VarmarginalPk60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MarginalPeak60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MarginalPeak60LayerBoundline.eps'));
  
  currfig=find_figure('MarginalPeaklayer10LayerBoundline')
xa=1:length(MmarginalPk10(1,:));
[l,p] = boundedline(xa,MmarginalPk10(1,:), VarmarginalPk10(1,:), '-g*', xa,MmarginalPk10(2,:), VarmarginalPk10(2,:), '-r*',xa,MmarginalPk10(3,:), VarmarginalPk10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MarginalPeaklayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MarginalPeaklayer10rBoundline.eps'));
  
  

currfig=find_figure(sprintf('marginalPk60',i))

  h=bar(MmarginalPk60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MmarginalPk60,1)
  numbars=size(MmarginalPk60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MmarginalPk60(:,jj),VarmarginalPk60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MmarginalPk60,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.D60Pk(:,end-4+jj);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_marginalPk60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarmarginalPk60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarmarginalPk60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('marginalPk10',i))

  h=bar(MmarginalPk10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(MmarginalPk10,1)
  numbars=size(MmarginalPk10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,MmarginalPk10(:,jj),VarmarginalPk10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(MmarginalPk10,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.D10Pk(:,end-ElementN+jj-1);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_marginalPk10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('BarmarginalPk10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('BarmarginalPk10_G%d.eps',NGroup)));





for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.D60Pk(:,end-L+i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('MarginalPeak60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.D60Pk(:,end-L+i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MarginalPeak60')
saveas(currfig,fullfile(savepath,sprintf('CdfMarginalPk%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMarginalPk%i.eps',60)));

%___________________________________________________________________________
for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.D10Pk(:,end-L+i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('MarginalPk10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.D10Pk(:,end-L+i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MarginalPk10')
saveas(currfig,fullfile(savepath,sprintf('CdfMarginalPk%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMarginalPk%i.eps',10)));


%---------------------------------------------------------------------------
%___________________________________________________________________________





Mmarginaldens60=[];Varmarginaldens60=[];Mmarginaldens10=[];Varmarginaldens10=[];
    for jd=1:NGroup
        marginaldens60=Cs{jd}.D60dens(:,end-ElementN:end);marginaldens10=Cs{jd}.D10dens(:,end-ElementN:end);
        
        marginaldens60(isnan(marginaldens60))=0;
         marginaldens10(isnan(marginaldens10))=0;
         mmarginaldens60=[];vmarginaldens60=[];mmarginaldens10=[];vmarginaldens10=[];
         for jlayer=1:size(Dist8060L,2)
             mmarginaldens60=[mmarginaldens60 nanmean(marginaldens60(:,jlayer))];
             vmarginaldens60=[vmarginaldens60 nanstd(marginaldens60(:,jlayer))./sqrt(size(marginaldens60,1))];
                  mmarginaldens10=[mmarginaldens10 nanmean(marginaldens10(:,jlayer))];
             vmarginaldens10=[vmarginaldens10 nanstd(marginaldens10(:,jlayer))./sqrt(size(marginaldens10,1))];
         end
      Mmarginaldens60=[Mmarginaldens60;mmarginaldens60];
     Varmarginaldens60=[Varmarginaldens60;vmarginaldens60]; 
      Mmarginaldens10=[Mmarginaldens10;mmarginaldens10];
     Varmarginaldens10=[Varmarginaldens10;vmarginaldens10]; 
      msfile=sprintf('MeanStd_marginaldens60_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mmarginaldens60','vmarginaldens60','-ascii') 
    
    msfile=sprintf('MeanStd_marginaldens10_%s.txt',GN{jd});
   msfile=fullfile(savepath,msfile);
    save(msfile,'mmarginaldens10','vmarginaldens10','-ascii') 
     
    end
      
  
  Mmarginaldens60=Mmarginaldens60';
  Varmarginaldens60=Varmarginaldens60';
  Mmarginaldens10=Mmarginaldens10';
  Varmarginaldens10=Varmarginaldens10';
  currfig=find_figure('MarginalDenslayer60LayerBoundline')
xa=1:length(MMpeaklayer60(1,:));
[l,p] = boundedline(xa,Mmarginaldens60(1,:), Varmarginaldens60(1,:), '-g*', xa,Mmarginaldens60(2,:), Varmarginaldens60(2,:), '-r*',xa,Mmarginaldens60(3,:), Varmarginaldens60(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MarginalDens60LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MarginalDens60LayerBoundline.eps'));
  
  currfig=find_figure('MarginalDenslayer10LayerBoundline')
xa=1:length(Mmarginaldens10(1,:));
[l,p] = boundedline(xa,Mmarginaldens10(1,:), Varmarginaldens10(1,:), '-g*', xa,Mmarginaldens10(2,:), Varmarginaldens10(2,:), '-r*',xa,Mmarginaldens10(3,:), Varmarginaldens10(3,:), '-b*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
saveas(currfig,fullfile(savepath,'MarginalDenslayer10LayerBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'MarginalDenslayer10rBoundline.eps'));
  
  

currfig=find_figure(sprintf('marginaldens60',i))

  h=bar(Mmarginaldens60);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mmarginaldens60,1)
  numbars=size(Mmarginaldens60,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mmarginaldens60(:,jj),Varmarginaldens60(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(Mmarginaldens60,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.D60dens(:,end-ElementN+jj-1);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_marginaldens60_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barmarginaldens60_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barmarginaldens60_G%d.eps',NGroup)));


currfig=find_figure(sprintf('marginaldens10',i))

  h=bar(Mmarginaldens10);
  for ih=1:NGroup
      h(ih).FaceColor=Cl(ih.*nl,:,:);
  end
   hold on
  set(h,'BarWidth',1);
  set(gca,'YGrid','on');
  set(gca,'GridLineStyle','-');
 set(gca, 'XTickLabel',{'L2/3','L4','L5/6'})
  numgroups=size(Mmarginaldens10,1)
  numbars=size(Mmarginaldens10,2);
  groupwidth=min(0.8,numbars/(numbars+1.5));
  XXaxis=[];
  for jj=1:numbars
      xxx=(1:numgroups)-groupwidth/2+(2*jj-1)*groupwidth/(2*numbars);
      XXaxis=[XXaxis; xxx];
      errorbar(xxx,Mmarginaldens10(:,jj),Varmarginaldens10(:,jj),'*','Color',Cl(jj.*nl,:,:))
  end
  
  for jj=1:size(Mmarginaldens10,1)
      Ds={};
      for ii=1:NGroup
          Ds{ii}=Cs{ii}.D10dens(:,end-ElementN+jj-1);
      end
      X1=XXaxis(:,jj);
      groups={};
      P=[];
      if NGroup>2
          p=StAnalysisMulti(Ds)
          for ii=1:size(p,1)
              groups{ii}=[X1(p(ii,1)) X1(p(ii,2))];
              P=[P,p(ii,end)];
          end
          
         
      else 
          [h,p]=StAnalysis(Ds{1},Ds{2},type);
          groups{1}=[X1(1),X1(2)];
          P=[P p];
      end
       H=sigstar(groups,P);
     pfile=sprintf('statisticP_marginaldens10_Layer%d.txt',jj);
   pfile=fullfile(savepath,pfile);
    save(pfile,'p','-ascii') 
    
     
  end
  
  
   
    
   saveas(currfig,fullfile(savepath,sprintf('Barmarginaldens10_G%d.fig',NGroup)),'fig');
   print(currfig,'-depsc2',fullfile(savepath,sprintf('Barmarginaldens10_G%d.eps',NGroup)));
















for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.D60dens(:,end-L+i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('MarginalArea60');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.D60dens(:,end-L+i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
     if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MarginalArea60')
saveas(currfig,fullfile(savepath,sprintf('CdfMarginalArea%i.fig',60)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMarginalArea%i.eps',60)));

%___________________________________________________________________________
for i=1:L
     Ds={};
    for ii=1:NGroup
        Ds{ii}=Cs{ii}.D10dens(:,end-L+i);
    end
    
    if NGroup>2
        c=StAnalysisMulti(Ds);
        
        
        
    else
        [hh,c]=StAnalysis(Ds{1},Ds{2},type);
        
    end   



 currfig=find_figure('MarginalArea10');

subplot(L,2,i*2-1)
hold on
for j=1:Ninput
h1=cdfplot(Cs{j}.D10dens(:,end-L+i));
  set(h1,'LineWidth',2,'Color',Cl(j.*nl,:,:))
end
legend(GN)
subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
end
 title('MarginalArea10')
saveas(currfig,fullfile(savepath,sprintf('CdfMarginalArea%i.fig',10)),'fig');
print(currfig,'-depsc2',fullfile(savepath,sprintf('CdfMarginalArea%i.eps',10)));




%-------------------------------------------------------




for i=1:Ninput
    Bnn=Cs{i}.Bnd;
 Rp{i}=deleteoutliers(abs(Bnn(:,3)./abs(Bnn(:,2)-Bnn(:,3))),0.05);
end

  
if Ninput==2
        [hh,c]=StAnalysis(Rp{1},Rp{2},type);
else c=StAnalysisMulti(Rp);
   
end

% Bn1=Cs{1}.Bnd;
% Bn2=Cs{2}.Bnd;
% Bn3=Cs{3}.Bnd;
% Bn4=Cs{4}.Bnd;
% Bn5=Cs{5}.Bnd;
% Rp1=deleteoutliers(abs(Bn1(:,4)./abs(Bn1(:,5)-Bn1(:,4))),0.05);
% Rp2=deleteoutliers(abs(Bn2(:,4)./abs(Bn2(:,5)-Bn2(:,4))),0.05);
% Rp3=deleteoutliers(abs(Bn3(:,4)./abs(Bn3(:,5)-Bn3(:,4))),0.05);    
% Rp4=deleteoutliers(abs(Bn4(:,4)./abs(Bn4(:,5)-Bn4(:,4))),0.05); 
% Rp5=deleteoutliers(abs(Bn5(:,4)./abs(Bn5(:,5)-Bn5(:,4))),0.05);  
% c=StAnalysis5( Rp1,Rp2,Rp3,Rp4,Rp5);
GG=[];XX=[];
for i=1:Ninput
Rp{i}=Rp{i}(Rp{i}<=1&Rp{i}>=0);Xrp=ones(size(Rp{i})).*i;XX=[XX;Xrp]
currfig=find_figure('position of cell in subplate')
subplot(1,2,1)
hold on; plot(Xrp,Rp{i},'o'); GG=[GG; Rp{i}]
end
hold on
boxplot(GG,XX,'Label',GN)
box off
subplot(1,2,2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end
 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off


saveas(currfig,fullfile(savepath,'Position.fig'),'fig');
print(currfig,'-depsc2',fullfile(savepath,'Position.eps'));

for i=1:L
   Rchg={};Rpk={};Rdens={};
for ii=1:Ninput
      
        c60total=Cs{ii}.chglayer60(:,i);
        c10total=Cs{ii}.chglayer10(:,i);
        d60total=Cs{ii}.arealayer60(:,i);
        d10total=Cs{ii}.arealayer10(:,i);
        ind10=find(c10total);
        rchg=c60total(ind10)./c10total(ind10);
        Rchg{ii}=rchg;
        p60total=Cs{ii}.peaklayer60(:,i);
        p10total=Cs{ii}.peaklayer10(:,i);
        ind10=find(p10total);
        rpk=p60total(ind10)./p10total(ind10);
        Rpk{ii}=rpk;
        
        ind10=find(d10total);
        rdens=d60total(ind10)./d10total(ind10);
        Rdens{ii}=rdens;
    end
    for ii=1:Ninput
        Rchg{ii}=deleteoutliers(Rchg{ii});
    end
    
    for ii=1:Ninput
        Rpk{ii}=deleteoutliers(Rpk{ii});
    end
    for ii=1:Ninput
        Rdens{ii}=deleteoutliers(Rdens{ii});
    end
    
    if NGroup>2
        c=StAnalysisMulti(Rchg);
        
        
        
    else
        [hh,c]=StAnalysis(Rchg{1},Rchg{2},type);
        
    end 
    
    
    currfig= find_figure('cdfEIRatioCharge'); 
   subplot(L,2,i*2-1)
    hold on
    for ii=1:Ninput
        
       h1=cdfplot(Rchg{ii})
        set(h1,'LineWidth',2,'Color',Cl(ii.*nl,:,:))
    end
     legend(GN)
    subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
 title('cdfEIRatioCharge')
%     if Ninput==2
%         title(sprintf(' p12=%i',c))
%         legend('CTRL','HI')
%         box off
%     elseif Ninput==3
%         title(sprintf(' p12=%i,p13=%i,p23=%i',c(1,end),c(2,end),c(3,end)))
%         legend('CTRL','HI','HYPX')
%         box off
%     elseif Ninput==4
%         title(sprintf(' p12=%i,p13=%i,p14=%i,p23=%i,p24=%i,p34=%i',c(1,end),c(2,end),c(3,end),c(4,end),c(5,end),c(6,end)))
%         legend('CTRL','HI','HYPXSham','HIligated')
%         box off
%     end
  
    saveas(currfig,fullfile(savepath,'CdfEIRatioCharge.fig'),'fig');
    print(currfig,'-depsc2',fullfile(savepath,'CdfEIRatioCharge.eps'));
    
    
    
   if NGroup>2
        c=StAnalysisMulti(Rpk);
        
        
        
    else
        [hh,c]=StAnalysis(Rpk{1},Rpk{2},type);
        
    end  
    
    
    
    
     currfig= find_figure('cdfEIRatioPeak'); 
   subplot(L,2,i*2-1)
    hold on
    for ii=1:Ninput
        
       h1=cdfplot(Rpk{ii})
        set(h1,'LineWidth',2,'Color',Cl(ii.*nl,:,:))
    end
     legend(GN)
    subplot(L,2,i*2)
sgm=eye(Ninput);
  for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
 title('cdfEIRatioPeak')
%     if Ninput==2
%         title(sprintf(' p12=%i',c))
%         legend('CTRL','HI')
%         box off
%     elseif Ninput==3
%         title(sprintf(' p12=%i,p13=%i,p23=%i',c(1,end),c(2,end),c(3,end)))
%         legend('CTRL','HI','HYPX')
%         box off
%     elseif Ninput==4
%         title(sprintf(' p12=%i,p13=%i,p14=%i,p23=%i,p24=%i,p34=%i',c(1,end),c(2,end),c(3,end),c(4,end),c(5,end),c(6,end)))
%         legend('CTRL','HI','HYPXSham','HIligated')
%         box off
%     end
  
    saveas(currfig,fullfile(savepath,'CdfEIRatioPeak.fig'),'fig');
    print(currfig,'-depsc2',fullfile(savepath,'CdfEIRatioPeak.eps'));
    
 
    
    if NGroup>2
        c=StAnalysisMulti(Rdens);
        
        
        
    else
        [hh,c]=StAnalysis(Rdens{1},Rdens{2},type);
        
    end  
    
    
    
    
     currfig= find_figure('cdfEIRatioDense'); 
   subplot(L,2,i*2-1)
    hold on
    for ii=1:Ninput
        
       h1=cdfplot(Rdens{ii})
        set(h1,'LineWidth',2,'Color',Cl(ii.*nl,:,:))
    end
     legend(GN)
    subplot(L,2,i*2)
sgm=eye(Ninput);
 for ind=1:size(c,1)
    if size(c,1)>1
         sgm(int8(c(ind,1)),int8(c(ind,2)))=c(ind,end);
         sgm(int8(c(ind,2)),int8(c(ind,1)))=c(ind,end);
     else
         sgm(1,2)=c;
     end 
 end
 imagesc(sgm);
 colormap('gray');
 textStrings = num2str(sgm(:));
 textStrings = strtrim(cellstr(textStrings));
 [x,y] = meshgrid(1:Ninput);   %# Create x and y coordinates for the strings
 hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
     'HorizontalAlignment','center');
 
 midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
 textColors = repmat(sgm(:)<midValue,1,3);  %# Choose white or black for the
 %#   text color of the strings so
 %#   they can be easily seen over
 %#   the background color
 set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
 
 set(gca,'XTick',1:Ninput,...                         %# Change the axes tick marks
     'XTickLabel',GN,...  %#   and tick labels
     'YTick',1:Ninput,...
     'YTickLabel',GN,...
     'TickLength',[0 0]);
 box off
 title('cdfEIRatioPeak')
%     if Ninput==2
%         title(sprintf(' p12=%i',c))
%         legend('CTRL','HI')
%         box off
%     elseif Ninput==3
%         title(sprintf(' p12=%i,p13=%i,p23=%i',c(1,end),c(2,end),c(3,end)))
%         legend('CTRL','HI','HYPX')
%         box off
%     elseif Ninput==4
%         title(sprintf(' p12=%i,p13=%i,p14=%i,p23=%i,p24=%i,p34=%i',c(1,end),c(2,end),c(3,end),c(4,end),c(5,end),c(6,end)))
%         legend('CTRL','HI','HYPXSham','HIligated')
%         box off
%     end
  
    saveas(currfig,fullfile(savepath,'CdfEIRatioDens.fig'),'fig');
    print(currfig,'-depsc2',fullfile(savepath,'CdfEIRatioDens.eps'));
    
    
end






colorLine={'k','g','b','r','y'}
currfig=find_figure('summery60')
hold on

for i=1:L
    subplot(4,5,1+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.arealayer60(:,i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end
for i=1:L
    subplot(4,5,2+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.chglayer60(:,i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end
for i=1:L
    subplot(4,5,3+5*(i-1))
    hold on
 for j=1:Ninput
   h1=cdfplot(Cs{j}.Mchglayer60(:,i));

  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end





for i=1:L
    subplot(4,5,4+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.D60dens(:,end-L+i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end


for i=1:L
    subplot(4,5,5+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.Pchglayer60(:,i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end


hold off

saveas(currfig,fullfile(savepath,'summery60.fig'),'fig');
print(currfig,'-depsc2',fullfile(savepath,'summery60.eps'));



currfig=find_figure('summery10')
hold on

for i=1:L
    subplot(4,5,1+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.arealayer10(:,i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end
for i=1:L
    subplot(4,5,2+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.chglayer10(:,i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end
for i=1:L
    subplot(4,5,3+5*(i-1))
    hold on
 for j=1:Ninput
   h1=cdfplot(Cs{j}.Mchglayer10(:,i));

  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end





for i=1:L
    subplot(L,5,4+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.D10dens(:,end-L+i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end

for i=1:L
    subplot(4,5,5+5*(i-1))
    hold on
 for j=1:Ninput
h1=cdfplot(Cs{j}.Pchglayer10(:,i));
  set(h1,'LineWidth',1,'Color',colorLine{j})
 end
end
hold off

saveas(currfig,fullfile(savepath,'summery10.fig'),'fig');
print(currfig,'-depsc2',fullfile(savepath,'summery10.eps'));


for i=1:L
   Rchg={};Rpk={};Rdens={};
for ii=1:Ninput
      
        c60total=Cs{ii}.chglayer60(:,i);
        c10total=Cs{ii}.chglayer10(:,i);
        d60total=Cs{ii}.arealayer60(:,i);
        d10total=Cs{ii}.arealayer10(:,i);
        ind10=find(c10total);
        rchg=c60total(ind10)./c10total(ind10);
        Rchg{ii}=rchg;
        p60total=Cs{ii}.peaklayer60(:,i);
        p10total=Cs{ii}.peaklayer10(:,i);
        ind10=find(p10total);
        rpk=p60total(ind10)./p10total(ind10);
        Rpk{ii}=rpk;
        
        ind10=find(d10total);
        rdens=d60total(ind10)./d10total(ind10);
        Rdens{ii}=rdens;
end

MRchg=[];
VarRchg=[];
MRpk=[];
VarRpk=[];
MRdens=[];
VarRdens=[];

 for ii=1:Ninput
        Rchg{ii}=deleteoutliers(Rchg{ii});
        MRchg=[MRchg;nanmean(Rchg{ii})];
VarRchg=[VarRchg;nanstd(Rchg{ii})./sqrt(sum(~isnan(Rchg{ii})))];
    end
    
    for ii=1:Ninput
        Rpk{ii}=deleteoutliers(Rpk{ii});
          MRpk=[MRpk;nanmean(Rpk{ii})];
VarRpk=[VarRpk;nanstd(Rchg{ii})./sqrt(sum(~isnan(Rpk{ii})))];
    end
    for ii=1:Ninput
        Rdens{ii}=deleteoutliers(Rdens{ii});
          MRdens=[MRdens;nanmean(Rdens{ii})];
VarRdens=[VarRdens;nanstd(Rdens{ii})./sqrt(sum(~isnan(Rdens{ii})))];
    end
    currfig=find_figure('EIratioBoundline')
   
    subplot(1,3,1)
xa=1:length(MRchg);
[l,p] = boundedline(xa,MRchg, VarRchg, '-g*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
subplot(1,3,2)
xa=1:length(MRpk);
[l,p] = boundedline(xa,MRpk, VarRpk, '-g*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
subplot(1,3,3)
xa=1:length(MRdens);
[l,p] = boundedline(xa,MRdens, VarRdens, '-g*');
outlinebounds(l,p);set(gca, 'XTick',1:Ninput,'XTickLabel',GN);set(gca, 'XTick',1:Ninput,'XTickLabel',GN)
    
saveas(currfig,fullfile(savepath,'EIratioBoundline.fig'),'fig');
   print(currfig,'-depsc2',fullfile(savepath,'EIratioBoundline.eps'));
  
  
    
   currfig= find_figure('cdfEIRatio');       
   subplot(L,4,4*(i-1)+1)
   hold on
   for ii=1:Ninput
   h1=cdfplot(Rchg{ii});
  set(h1,'LineWidth',1,'Color',colorLine{ii})
   end
  
   subplot(L,4,4*(i-1)+2)
   hold on
   for ii=1:Ninput
   h1=cdfplot(Rpk{ii});
  set(h1,'LineWidth',1,'Color',colorLine{ii})
   end
  
    subplot(L,4,4*(i-1)+3)
    hold on
  for ii=1:Ninput
   h1=cdfplot(Rdens{ii});
  set(h1,'LineWidth',1,'Color',colorLine{ii})
   end
end
saveas(currfig,fullfile(savepath,'cdfEIRatio.fig'),'fig');
print(currfig,'-depsc2',fullfile(savepath,'cdfEIRatio.eps'));
