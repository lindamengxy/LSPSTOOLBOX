function  handles= mapave_gui2(datafolder,handles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

folderMeanname=sprintf('meanValueeachcell_eventwindow%i_direct%i',handles.eventwindow,handles.direct_t);
folderMeanpath=fullfile(datafolder,folderMeanname);
if (exist(folderMeanpath) == 0)
    mkdir (folderMeanpath);
end
if ~isfield(handles,'eventwindowValue')||isempty(handles.eventwindowValue)
        handles.eventwindow=50;
end
foldername=sprintf('%s%i','eventwindow',handles.eventwindow);
datafolder=fullfile(datafolder,foldername);
a=dir(fullfile(datafolder,'*.mat'));
Cellfile=size(a);

C=cell(1,1);

% xaxis=[];
% yaxis=[];
% Pk=[];
% Chg=[];
% Den=[];
% Lat=[];
% dX=[];
% AGE=[];
% N=0;
% if strcmp(handles.experimenttype,'Whole cell')



Pth=cell(Cellfile);
Tp=[];
Hpotential=[];
flag=[];
for fn= 1:1:Cellfile
    C=load(fullfile(datafolder,a(fn).name));
    if isfield(C,'cells')
        cells=C.cells;
        Pth{fn}=cells.cell_folder;
        Tp=[Tp,cells.Tp];
        Hpotential=[Hpotential,cells.header.holdingPotential];
        flag=[flag,1];
    else flag=[flag,0];Tp=[Tp,NaN];Hpotential=[Hpotential,NaN];Pth{fn}=NaN;
    end
end

NC=0;
ff=0;
Ind=find(flag);
folderCellNum=5;
totalCellNum=0;
while length(Ind)
    
    ff=ff+1;
    Num_cell1=0;
    
    Cellrecord=cell(folderCellNum,1);
    
    while Num_cell1<folderCellNum &&length(Ind)
        
        fn=Ind(1);
        pt=Pth{fn};
        bs=find(strcmp(pt,Pth));
        HighMgInd60=[];
        HighMgInd10=[];
        PtxInd60=[];
        PtxInd50=[];
        PtxInd10=[];
        apvInd60=[];
        apvInd10=[];
        apvInd50=[];
        ttxInd60=[];
        ttxInd10=[];
        ttxInd50=[];
        otherInd60=[];
        otherInd10=[];
        BUD=[];
        if length(bs)
            Num_cell1=Num_cell1+1;
            totalCellNum=totalCellNum+1
            
            C1=load(fullfile(datafolder,a(bs(1)).name));
            for j=1:length(bs)
                
                
                CC=load(fullfile(datafolder,a(bs(j)).name));
                if isfield(CC.cells,'Boundary')
                    
                BUD=[BUD;CC.cells.Boundary(1,:)];
                else
                    Boundry=boundary(CC.cells.header,CC.cells.flipimg);
                     BUD=[BUD;Boundry(1,:)];
                end
                if CC.cells.header.holdingPotential<=-50&CC.cells.Tp==1
                    
                    HighMgInd60=[HighMgInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==1
                    HighMgInd10=[HighMgInd10,bs(j)];
                elseif CC.cells.header.holdingPotential<=-50&CC.cells.Tp==2
                    PtxInd60=[PtxInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==2
                    PtxInd10=[PtxInd10,bs(j)];
                elseif CC.cells.header.holdingPotential>30&CC.cells.Tp==2
                    PtxInd50=[PtxInd50,bs(j)];
                elseif CC.cells.header.holdingPotential<=-50&CC.cells.Tp==3
                    apvInd60=[apvInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==3
                    apvInd10=[apvInd10,bs(j)];
                elseif CC.cells.header.holdingPotential>30&CC.cells.Tp==3
                    apvInd50=[apvInd50,bs(j)];
                elseif CC.cells.header.holdingPotential<=-50&CC.cells.Tp==4
                    ttxInd60=[ttxInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==4
                    ttxInd10=[ttxInd10,bs(j)];
                elseif CC.cells.header.holdingPotential>30&CC.cells.Tp==4
                    ttxInd50=[ttxInd50,bs(j)];
                elseif CC.cells.header.holdingPotential<=-50&CC.cells.Tp==5
                    ttxInd60=[otherInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==5
                    otherInd10=[otherInd10,bs(j)];
                    
                end
            end
             if length(HighMgInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] =  Meanvalue_Gui_2(a,datafolder,HighMgInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyHighMg60',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurHighMg60',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakHighMg60',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauHighMg60',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaHighMg60',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagHighMg60',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanrateHighMg60',Num_cell1) '=meanrate' ';']);
            eval([sprintf('Cellrecord{%i}.meanflagDHighMg60',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyHighMg60',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaHighMg60',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagHighMg60',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinateHighMgInd10=HighMgInd10(1);
            end
            if length(HighMgInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] =  Meanvalue_Gui_2(a,datafolder,HighMgInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyHighMg10',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurHighMg10',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakHighMg10',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauHighMg10',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaHighMg10',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagHighMg10',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanrateHighMg10',Num_cell1) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDHighMg10',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyHighMg10',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaHighMg10',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagHighMg10',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinateHighMgInd10=HighMgInd10(1);
            end
            if length(PtxInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] =  Meanvalue_Gui_2(a,datafolder,PtxInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyPtxInd60',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurPtxInd60',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakPtxInd60',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauPtxInd60',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaPtxInd60',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagPtxInd60',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanratePtxInd60',Num_cell1) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDPtxInd60',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyPtxInd60',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaPtxInd60',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagPtxInd60',Num_cell1) '=Corrflag' ';']);
                %                      Cellrecord{Num_cell}.stimcoordinatePtxInd60=PtxInd60(1);
            end
            if length(PtxInd50)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] =  Meanvalue_Gui_2(a,datafolder,PtxInd50,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyPtxInd50',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurPtxInd50',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakPtxInd50',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauPtxInd50',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaPtxInd50',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagPtxInd50',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanratePtxInd50',Num_cell1) '=meanrate' ';']);
                 eval([sprintf('Cellrecord{%i}.meanflagDPtxInd50',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyPtxInd50',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaPtxInd50',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagPtxInd50',Num_cell1) '=Corrflag' ';']);
                %                    Cellrecord{Num_cell}.stimcoordinatePtxInd50=PtxInd50(1);
            end
            if length(PtxInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,PtxInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyPtxInd10',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurPtxInd10',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakPtxInd10',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauPtxInd10',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaPtxInd10',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagPtxInd10',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDPtxInd10',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyPtxInd10',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaPtxInd10',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagPtxInd10',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(apvInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,apvInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyapvInd60',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurapvInd60',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakapvInd60',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauapvInd60',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaapvInd60',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagapvInd60',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDapvInd60',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyapvInd60',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaapvInd60',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagapvInd60',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(apvInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,apvInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyapvInd10',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurapvInd10',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakapvInd10',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauapvInd10',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaapvInd10',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagapvInd10',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDapvInd10',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyapvInd10',Num_cell1) '=mDirlatency' ';']);
                eval([sprintf('Cellrecord{%i}.CorrflagapvInd10',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(apvInd50)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,apvInd50,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyapvInd50',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurapvInd50',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakapvInd50',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauapvInd50',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaapvInd50',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagapvInd50',Num_cell1) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanflagDapvInd50',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyapvInd50',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaapvInd50',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagapvInd50',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(ttxInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,ttxInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
  
                eval([sprintf('Cellrecord{%i}.meanlatencyttxInd60',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurttxInd60',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakttxInd60',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauttxInd60',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareattxInd60',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagttxInd60',Num_cell1) '=meanflag' ';']);
                  eval([sprintf('Cellrecord{%i}.meanflagDttxInd60',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyttxInd60',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareattxInd60',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagttxInd60',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
                end
            if length(ttxInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,ttxInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                
                eval([sprintf('Cellrecord{%i}.meanlatencyttxInd10',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurttxInd10',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakttxInd10',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauttxInd10',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareattxInd10',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagttxInd10',Num_cell1) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanflagDttxInd10',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyttxInd10',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareattxInd10',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagttxInd10',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(ttxInd50)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,ttxInd50,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                
                eval([sprintf('Cellrecord{%i}.meanlatencyttxInd50',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurttxInd50',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakttxInd50',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauttxInd50',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareattxInd50',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagttxInd50',Num_cell1) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanflagDttxInd50',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyttxInd50',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareattxInd50',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagttxInd50',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(otherInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,otherInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                
                eval([sprintf('Cellrecord{%i}.meanlatencyotherInd60',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurotherInd60',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakotherInd60',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauotherInd60',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaotherInd60',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagotherInd60',Num_cell1) '=meanflag' ';']);
                   eval([sprintf('Cellrecord{%i}.meanflagDotherInd60',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyotherInd60',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaotherInd60',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagotherInd60',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(otherInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea,Corrflag] = Meanvalue_Gui_2(a,datafolder,otherInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                
                eval([sprintf('Cellrecord{%i}.meanlatencyotherInd10',Num_cell1) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurotherInd10',Num_cell1) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakotherInd10',Num_cell1) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauotherInd10',Num_cell1) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaotherInd10',Num_cell1) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagotherInd10',Num_cell1) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDotherInd10',Num_cell1) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyotherInd10',Num_cell1) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaotherInd10',Num_cell1) '=mDirarea' ';']);
                 eval([sprintf('Cellrecord{%i}.CorrflagotherInd10',Num_cell1) '=Corrflag' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            Cellrecord{Num_cell1}.SomaCoordinates=C1.cells.header.Soma1Coordinates;
            
            
            Cellrecord{Num_cell1}.StimCoordinates=C1.cells.header.StimCoordinates;
            Cellrecord{Num_cell1}.SpatialRotation=C1.cells.header.SpatialRotation;
            Cellrecord{Num_cell1}.PatternSpacing=C1.cells.header.PatternSpacing;
            Cellrecord{Num_cell1}.age=C1.cells.age;
            Cellrecord{Num_cell1}.darkexp=C1.cells.darkexp;
            Cellrecord{Num_cell1}.flipimg=C1.cells.flipimg;
            Cellrecord{Num_cell1}.Pth=pt;
           Cellrecord{Num_cell1}.Boundry=mean(BUD,1);
           Cellrecord{Num_cell1}.ImageData=C1.cells.header.ImageData;
           Cellrecord{Num_cell1}.stim_start=find(C1.cells.header.stimulatorOutput(4,:));
           Cellrecord{Num_cell1}.ImageScale=C1.cells.header.ImageScale;
           
           
            if ~isfield(C1.cells,'experimenttypeValue')||isempty(C1.cells.experimenttypeValue)
                Cellrecord{Num_cell1}.experimenttype='Whole cell';
            else
                Cellrecord{Num_cell1}.experimenttype=C1.cells.experimenttype;
                
            end
        end
        flag(bs)=0;
        Ind=find(flag);
        
    end
    Num_cell=Num_cell1-1;
   
    filename=sprintf('MeanValue%i',ff);
    filepath=fullfile(folderMeanpath,filename)
    save(filepath,'Cellrecord','Num_cell','filepath');
end

handles.folderMeanpath=folderMeanpath;
end
