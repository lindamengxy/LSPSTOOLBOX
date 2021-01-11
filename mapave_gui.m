function  handles= mapave_gui(datafolder,handles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


a=dir(fullfile(datafolder,'*.mat'));
Cellfile=size(a);
C=cell(1,1);  

xaxis=[];
yaxis=[];
Pk=[];
Chg=[];
Den=[];
Lat=[];
dX=[];
AGE=[];
N=0;
% if strcmp(handles.experimenttype,'Whole cell')
    
  if ~isfield(handles,'avgTp')||isempty(handles.avgTp)||ishandle(handles.avgTp)
        handles.avgTp=1;
        set(handles.error_run,'string','Forget to choose the drugs type. The default version is High Mg. If not, please choose the right one')
        else
        set(handles.error_run,'string','')
  end
  if ~isfield(handles,'minage')||isempty(handles.minage)||ishandle(handles.minage)
        handles.minage=0;
        
  end
  if ~isfield(handles,'cutoffdensity')||isempty(handles.cutoffdensity)||ishandle(handles.cutoffdensity)
        handles.cutoffdensity=0;
        
  end
    if ~isfield(handles,'maxage')||isempty(handles.maxage)||ishandle(handles.maxage)
        handles.maxage=Inf;
  end
  if ~isfield(handles,'avgdarkexp')||isempty(handles.avgdarkexp)||ishandle(handles.avgdarkexp)
        handles.avgdarkexp=0;
        
  end
  if ~isfield(handles,'avgHoldingPotential')||isempty(handles.avgHoldingPotential)||ishandle(handles.avgHoldingPotential)
        handles.avgHoldingPotential=-60;
        set(handles.error_run,'string','Forget to fill in the holding potential. The default value is -60. If not, please choose the right one')
        else
        set(handles.error_run,'string','')
  end
  
Pth=cell(Cellfile);
Tp=[];
Hpotential=[];
flag=[];
for fn= 1:1:Cellfile
     C=load(fullfile(datapath,a(fn).name));
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
    Num_cell=0;
    
    Cellrecord=cells(folderCellNum,1);
    
    while Num_cell<=folderCellNum
        
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
        
        if length(bs)
            Num_cell=Num_cell+1;
            totalCellNum=totalCellNum+1
            C1=load(fullfile(datapath,a(bs(1)).name));
            for j=1:length(bs)
                
                CC=load(fullfile(datapath,a(bs(j)).name));
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
                    PtxInd60=[apvInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==3
                    PtxInd10=[apvInd10,bs(j)];
                elseif CC.cells.header.holdingPotential>30&CC.cells.Tp==3
                    PtxInd50=[apvInd50,bs(j)];
                elseif CC.cells.header.holdingPotential<=-50&CC.cells.Tp==4
                    PtxInd60=[ttxInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==4
                    PtxInd10=[ttxInd10,bs(j)];
                elseif CC.cells.header.holdingPotential>30&CC.cells.Tp==4
                    PtxInd50=[ttxInd50,bs(j)];
                elseif CC.cells.header.holdingPotential<=-50&CC.cells.Tp==5
                    PtxInd60=[otherInd60,bs(j)];
                elseif CC.cells.header.holdingPotential>=-20&CC.cells.header.holdingPotential<30&CC.cells.Tp==5
                    PtxInd10=[otherInd10,bs(j)];
                    
                end
            end
             if length(HighMgInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] =  Meanvalue_Gui_2(a,datapath,HighMgInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyHighMg60',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurHighMg60',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakHighMg60',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauHighMg60',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaHighMg60',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagHighMg60',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanrateHighMg60',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDHighMg60',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyHighMg60',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaHighMg60',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinateHighMgInd10=HighMgInd10(1);
            end
            if length(HighMgInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] =  Meanvalue_Gui_2(a,datapath,HighMgInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyHighMg10',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurHighMg10',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakHighMg10',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauHighMg10',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaHighMg10',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagHighMg10',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanrateHighMg10',Num_cell) '=meanrate' ';']);
                 eval([sprintf('Cellrecord{%i}.meanflagDHighMg10',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyHighMg10',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaHighMg10',Num_cell) '=mDirarea' ';']);
                
                %                     Cellrecord{Num_cell}.stimcoordinateHighMgInd10=HighMgInd10(1);
            end
            if length(PtxInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] =  Meanvalue_Gui_2(a,datapath,PtxInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyPtxInd60',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurPtxInd60',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakPtxInd60',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauPtxInd60',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaPtxInd60',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagPtxInd60',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanratePtxInd60',Num_cell) '=meanrate' ';']);
                 eval([sprintf('Cellrecord{%i}.meanflagDPtxInd60',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyPtxInd60',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaPtxInd60',Num_cell) '=mDirarea' ';']);
                %                      Cellrecord{Num_cell}.stimcoordinatePtxInd60=PtxInd60(1);
            end
            if length(PtxInd50)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] =  Meanvalue_Gui_2(a,datapath,PtxInd50,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyPtxInd50',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurPtxInd50',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakPtxInd50',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauPtxInd50',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaPtxInd50',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagPtxInd50',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanratePtxInd50',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDPtxInd50',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyPtxInd50',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaPtxInd50',Num_cell) '=mDirarea' ';']);
                %                    Cellrecord{Num_cell}.stimcoordinatePtxInd50=PtxInd50(1);
            end
            if length(PtxInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(a,datapath,PtxInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyPtxInd10',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurPtxInd10',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakPtxInd10',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauPtxInd10',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaPtxInd10',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagPtxInd10',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanratePtxInd10',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDPtxInd10',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyPtxInd10',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaPtxInd10',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(apvInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(a,datapath,apvInd60,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyapvInd60',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurapvInd60',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakapvInd60',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauapvInd60',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaapvInd60',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagapvInd60',Num_cell) '=meanflag' ';']);
                  eval([sprintf('Cellrecord{%i}.meanrateapvInd60',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDapvInd60',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyapvInd60',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaapvInd60',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(apvInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(a,datapath,apvInd10,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyapvInd10',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurapvInd10',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakapvInd10',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauapvInd10',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaapvInd10',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagapvInd10',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanrateapvInd10',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDapvInd10',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyapvInd10',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaapvInd10',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(apvInd50)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(a,datapath,apvInd50,C1.cells.direct_t,C1.cells.direct_t,C1.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyapvInd50',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurapvInd50',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakapvInd50',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauapvInd50',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaapvInd50',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagapvInd50',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanrateapvInd50',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDapvInd50',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyapvInd50',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaapvInd50',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(ttxInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(C,ttxInd60,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyttxInd60',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurttxInd60',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakttxInd60',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauttxInd60',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareattxInd60',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagttxInd60',Num_cell) '=meanflag' ';']);
                eval([sprintf('Cellrecord{%i}.meanratettxInd60',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDttxInd60',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyttxInd60',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareattxInd60',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(ttxInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(C,ttxInd10,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyttxInd10',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurttxInd10',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakttxInd10',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauttxInd10',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareattxInd10',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagttxInd10',Num_cell) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanratettxInd10',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDttxInd10',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyttxInd10',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareattxInd10',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(ttxInd50)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(C,ttxInd50,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyttxInd50',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurttxInd50',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakttxInd50',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauttxInd50',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareattxInd50',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagttxInd50',Num_cell) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanratettxInd50',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDttxInd50',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyttxInd50',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareattxInd50',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(otherInd60)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(C,otherInd60,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyotherInd60',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurotherInd60',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakotherInd60',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauotherInd60',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaotherInd60',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagotherInd60',Num_cell) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanrateotherInd60',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDotherInd60',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyotherInd60',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaotherInd60',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            if length(otherInd10)
                [meanlatency,meantau,meandur,meanarea,meanpeak,meanflag,meanrate,meanflagD,mDirlatency,mDirarea] = Meanvalue_Gui_2(C,otherInd10,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.direct_t,C{bs(1)}.cells.eventwindow,bs(1));
                
                eval([sprintf('Cellrecord{%i}.meanlatencyotherInd10',Num_cell) '=meanlatency' ';']);
                eval([sprintf('Cellrecord{%i}.meandurotherInd10',Num_cell) '=meandur' ';']);
                eval([sprintf('Cellrecord{%i}.meanpeakotherInd10',Num_cell) '=meanpeak' ';']);
                eval([sprintf('Cellrecord{%i}.meantauotherInd10',Num_cell) '=meantau' ';']);
                eval([sprintf('Cellrecord{%i}.meanareaotherInd10',Num_cell) '=meanarea' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagotherInd10',Num_cell) '=meanflag' ';']);
                 eval([sprintf('Cellrecord{%i}.meanrateotherInd10',Num_cell) '=meanrate' ';']);
                eval([sprintf('Cellrecord{%i}.meanflagDotherInd10',Num_cell) '=meanflagD' ';']);
                eval([sprintf('Cellrecord{%i}.meanDirlatencyotherInd10',Num_cell) '=mDirlatency' ';']);
                 eval([sprintf('Cellrecord{%i}.meanDirareaotherInd10',Num_cell) '=mDirarea' ';']);
                %                     Cellrecord{Num_cell}.stimcoordinatePtxInd10=PtxInd10(1);
            end
            Cellrecord{Num_cell}.SomaCoordinates=C1.cells.header.Soma1Coordinates;
            
            
            Cellrecord{Num_cell}.StimCoordinates=C1.cells.header.StimCoordinates;
            Cellrecord{Num_cell}.SpatialRotation=C1.cells.header.SpatialRotation;
            Cellrecord{Num_cell}.PatternSpacing=C1.cells.header.PatternSpacing;
            Cellrecord{Num_cell}.age=C1.cells.age;
            Cellrecord{Num_cell}.darkexp=C1.cells.darkexp;
            Cellrecord{Num_cell}.flipimg=C1.cells.flipimg;
            Cellrecord{Num_cell}.Pth=pt;
            if ~isfield(C1.cells,'experimenttype')|isempty(C1.cells.experimenttype)
                Cellrecord{Num_cell}.experimenttype='Whole cell';
            else
                Cellrecord{Num_cell}.experimenttype;
                
            end
        end
        flag(bs)=0;
        Ind=find(flag);
        
    end
    folderMeanname=sprintf('meanValueeachcell_eventwindow%i_direct%i',handles.eventwindow,handles.direct_t);
    filepath=fullfile(datafolder,folderMeanname);
    if (exist(filepath) == 0)
        mkdir (filepath);
    end
    filename=sprintf('MeanValue%i',ff);
    filepath=fullfile(filepath,filename)
    save(filepath,'Cellrecord','Num_cell','filepath');
end


     
   direct_t=handles.direct_t;
   eventwindow=handles.eventwindow;
  for fn=1:1:Num_cell
      
      
%       if isfield(C{fn},'cells')
%       cells=C{fn}.cells;
%       N=N+1;
      h=Cellrecord{fn}.experimenttype;
     
%       header=cells.header;
%       events=cells.events;
%       stim_start=find(header.stimulatorOutput(4,:));
      if Cellrecord{fn}.age>=handles.minage&Cellrecord{fn}.age<handles.maxage&abs(Cellrecord{fn}.darkexp-handles.avgdarkexp)<0.1&strcmp(h,'Whole cell')
          Soma= Cellrecord{fn}.SomaCoordinates;
        flipimg=Cellrecord{fn}.flipimg;
        pth=char(Cellrecord{fn}.Pth);
        stimcoordinates=Cellrecord{fn}.StimCoordinates;
        rotate_angle=Cellrecord{fn}.SpatialRotation;
        [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
        X=stimCoordinates(1,:);
        Y=stimCoordinates(2,:);
          xaxis=[xaxis,X];
          yaxis=[yaxis,Y];
          AGE=[AGE,Cellrecord{fn}.age.*ones(1,header.nPts)];
          dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
          
          if 
          
              Zpk=events.peakAmp;
         
              Zchg=events.area;
         
              Zlat=events.startSamp;
        
          z1=[];z2=[];z3=[];z4=[];
          for k=1:1:header.nPts
              fevent=find(events.flag{k}(:)>0);
              
              if length(fevent)>0
                  if events.flag{k}(fevent(1))==1
                      
                     
                      z1=[z1 NaN];
                      
                      z2=[z2 NaN];
                      
                      z3=[z3 NaN];
                      z4=[z4 0];
                  elseif events.flag{k}(fevent(1))==2
                      z1=[z1 Zpk{k}(fevent(1))]; 
                      z2=[z2 Zchg{k}(fevent(1))]; 
                      z3=[z3 Zlat{k}(fevent(1))]; 
                      z4=[z4 1];
                  end
              else z1=[z1 NaN];
                   z2=[z2 NaN];
                   z3=[z3 NaN];
                   z4=[z4 0];
              end
          end
          
       
              Pk=[Pk,z1];
              
          
              Chg=[Chg,z2];
          
              Lat=[Lat,z3.*1000./header.sampleRate-stim_start(1).*1000./header.sampleRate];
              Den=[Den z4];
      end
      end
      
  end
  handles.avgPk=Pk;
  handles.avgChg=Chg;
  handles.Lat=Lat;
  handles.dX=dX;
  handles.AvgAge=AGE;
  handles.avgxaxis=xaxis;
  handles.avgyaxis=yaxis;
  handles.avgDensity=Den;
  if N<1
     set(handles.error_run,'string','No maps meets your standards')
        else
        set(handles.error_run,'string','')
  end
% end
end
