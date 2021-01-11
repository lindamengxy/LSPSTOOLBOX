
La=[];DU=[];
    AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];xaxis60=[];yaxis60=[];
  load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1634.mat'  
        Soma= cells.header.Soma1Coordinates;
        flipimg=cells.flipimg;
        
      header60=cells.header;
        stimcoordinates=cells.header.StimCoordinates;
        rotate_angle=cells.header.SpatialRotation;
        [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
        xaxis60=stimCoordinates(1,:);
        yaxis60=stimCoordinates(2,:);
        
    for k=1:cells.header.nPts
    
    %for dd=1:1:length(Indx)
        %CL=CC{dd};
        header1=cells.header;
        stim_start=find(header1.stimulatorOutput(4,:));
        
        %xxx=INMatrix(k,dd);
       
        fevent=find(cells.events.flag{k}(:)>0);
     if length(fevent)
            if cells.events.flag{k}(fevent(1))==1
                ff=NaN;
            elseif cells.events.flag{k}(fevent(1))==2
                ff=fevent(1);
                
            end
        else ff=NaN;
        end
        
        
        if ~isnan(ff)
            latency=cells.events.startSamp{k}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
            dur=cells.events.duration{k}(ff).*1000./cells.header.sampleRate;
            area=(cells.events.area{k}(ff));
            peak=(cells.events.peakAmp{k}(ff));
            tau=(cells.events.tau{k}(ff));
            flag_latency=1;
            
            %            eventcnt=eventcnt+1;
            
        else
            
            latency=NaN;
            dur=NaN;
            area=NaN;
            peak=NaN;
            tau=NaN;
            flag_latency=0;
        end
        
        
        La=[La,latency];
        DU=[DU dur];
        AREA=[AREA area];
        PA=[PA peak];
        TAU=[TAU tau];
        Flg=[Flg flag_latency];
    end
nanmean(PA)



Chg60=AREA;
Pk60=PA;
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1706.mat'
La=[];DU=[];
    AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
   
       header1=cells.header;
        X=sameind(header60,header1);
    for k=1:cells.header.nPts
    
    %for dd=1:1:length(Indx)
        %CL=CC{dd};
     
        stim_start=find(header1.stimulatorOutput(4,:));
        
        %xxx=INMatrix(k,dd);
       
        fevent=find(cells.events.flag{X(k)}(:)>0);
        
      if length(fevent)
            if cells.events.flag{X(k)}(fevent(1))==1
                ff=NaN;
            elseif cells.events.flag{X(k)}(fevent(1))==2
                ff=fevent(1);
                
            end
        else ff=NaN;
        end
        
        
        
        if ~isnan(ff)
            latency=cells.events.startSamp{X(k)}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
            dur=cells.events.duration{X(k)}(ff).*1000./cells.header.sampleRate;
            area=(cells.events.area{X(k)}(ff));
            peak=(cells.events.peakAmp{X(k)}(ff));
            tau=(cells.events.tau{X(k)}(ff));
            flag_latency=1;
            
            %            eventcnt=eventcnt+1;
            
        else
            
            latency=NaN;
            dur=NaN;
            area=NaN;
            peak=NaN;
            tau=NaN;
            flag_latency=0;
        end
        
        
        La=[La,latency];
        DU=[DU dur];
        AREA=[AREA area];
        PA=[PA peak];
        TAU=[TAU tau];
        Flg=[Flg flag_latency];
    end
PKAPV=PA;
ChgAPV=AREA;

PkAPV60=PKAPV;
ChgAPV60=ChgAPV;
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1639.mat'
La=[];DU=[];
    AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
       header1=cells.header;
        X=sameind(header60,header1);
    for k=1:cells.header.nPts
    
    %for dd=1:1:length(Indx)
        %CL=CC{dd};
       
        stim_start=find(header1.stimulatorOutput(4,:));
        
        %xxx=INMatrix(k,dd);
       
        fevent=find(cells.events.flag{X(k)}(:)>0);
        
        if length(fevent)
            if cells.events.flag{X(k)}(fevent(1))==1
                ff=NaN;
            elseif cells.events.flag{X(k)}(fevent(1))==2
                ff=fevent(1);
                
            end
        else ff=NaN;
        end
        
        
        if ~isnan(ff)
            latency=cells.events.startSamp{X(k)}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
            dur=cells.events.duration{X(k)}(ff).*1000./cells.header.sampleRate;
            area=(cells.events.area{X(k)}(ff));
            peak=(cells.events.peakAmp{X(k)}(ff));
            tau=(cells.events.tau{X(k)}(ff));
            flag_latency=1;
            
            %            eventcnt=eventcnt+1;
            
        else
            
            latency=NaN;
            dur=NaN;
            area=NaN;
            peak=NaN;
            tau=NaN;
            flag_latency=0;
        end
        
        
        La=[La,latency];
        DU=[DU dur];
        AREA=[AREA area];
        PA=[PA peak];
        TAU=[TAU tau];
        Flg=[Flg flag_latency];
    end
Pk40=PA;
Chg40=AREA;
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1711.mat'
La=[];DU=[];
    AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
       header1=cells.header;
        X=sameind(header60,header1);
    for k=1:cells.header.nPts
    
    %for dd=1:1:length(Indx)
        %CL=CC{dd};
        
        stim_start=find(header1.stimulatorOutput(4,:));
        
        %xxx=INMatrix(k,dd);
       
        fevent=find(cells.events.flag{X(k)}(:)>0);
        
       if length(fevent)
            if cells.events.flag{X(k)}(fevent(1))==1
                ff=NaN;
            elseif cells.events.flag{X(k)}(fevent(1))==2
                ff=fevent(1);
                
            end
        else ff=NaN;
        end
        
        
        if ~isnan(ff)
            latency=cells.events.startSamp{X(k)}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
            dur=cells.events.duration{X(k)}(ff).*1000./cells.header.sampleRate;
            area=(cells.events.area{X(k)}(ff));
            peak=(cells.events.peakAmp{X(k)}(ff));
            tau=(cells.events.tau{X(k)}(ff));
            flag_latency=1;
            
            %            eventcnt=eventcnt+1;
            
        else
            
            latency=NaN;
            dur=NaN;
            area=NaN;
            peak=NaN;
            tau=NaN;
            flag_latency=0;
        end
        
        
        La=[La,latency];
        DU=[DU dur];
        AREA=[AREA area];
        PA=[PA peak];
        TAU=[TAU tau];
        Flg=[Flg flag_latency];
    end
Pk40APV=PA;
Chg40APV=AREA;
mPK60=nanmean(Pk60);
mChg60=nanmean(Chg60);
mPK40=nanmean(Pk40);
mChg40=nanmean(Chg40);

mPK60APV=nanmean(PkAPV60);
mChg60APV=nanmean(ChgAPV60);
mPK40APV=nanmean(Pk40APV);
mChg40APV=nanmean(Chg40APV);
vPK60=nanstd(Pk60);
vChg60=nanstd(Chg60);
vPK40=nanstd(Pk40);
vChg40=nanstd(Chg40);
vPK60APV=nanstd(PkAPV60);
vChg60APV=nanstd(ChgAPV60);
vPK40APV=nanstd(Pk40APV);
vChg40APV=nanstd(Chg40APV);
figure(1)
 x=[0.86 1.14 1.86 2.14];
errorbar(x,[mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV],'b*')
hold on
bar([mPK60 mPK60APV ; mPK40 mPK40APV])

figure(3)
x=[0.86 1.14 1.86 2.14];
errorbar(x,[mChg60 mChg60APV mChg40 mChg40APV],[vChg60 vChg60APV vChg40 vChg40APV],'b*')
hold on

bar([mChg60 mChg60APV ; mChg40 mChg40APV])
PkAPV60(find(isnan(PkAPV60)))=0;
Pk60(find(isnan(Pk60)))=0;
Pk40APV(find(isnan(Pk40APV)))=0;
Pk40(find(isnan(Pk40)))=0;
DifPK60=(Pk60-PkAPV60)./Pk60;
DifPK40=(Pk40-Pk40APV)./Pk40;

ChgAPV60(find(isnan(ChgAPV60)))=0;
Chg60(find(isnan(Chg60)))=0;
Chg40APV(find(isnan(Chg40APV)))=0;
Chg40(find(isnan(Chg40)))=0;
DifChg60=(Chg60-ChgAPV60)./Chg60;
DifChg40=(Chg40-Chg40APV)./Chg40;


dx = cells.header.PatternSpacing(2); dy = cells.header.PatternSpacing(1);N=1;density=0;
        xlim = [ceil(min(xaxis60))-dx-0.1 ceil(max(xaxis60))]; ylim = [ceil(min(yaxis60))-dy-0.1 ceil(max(yaxis60))];
        
        xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
        xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
        pos = [xaxis60;yaxis60];
        xii=xi(1):1:xi(end);
    yii=yi(1):1:yi(end);
    [yi1,xi1]=meshgrid(yi,xi);
    [yii,xii]=meshgrid(yii,xii);
         avg = gridavg2(pos(1,:),pos(2,:),DifPK60,xtick,ytick,false,N,density);
         
         avg=avg';
    currfig=find_figure('different_map-60')

    avg2=interp2(yi1,xi1,avg,yii,xii);
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yi,xi,avg); axis image;axis ij;
    colormap(jet)
    colorbar
    currfig=find_figure('different_map40')
    avg = gridavg2(pos(1,:),pos(2,:),DifPK40,xtick,ytick,false,N,density);
    avg=avg';     
    avg2=interp2(yi1,xi1,avg,yii,xii);
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yi,xi,avg); axis image;axis ij;
    colormap(jet)
    colorbar
    
    
     avg = gridavg2(pos(1,:),pos(2,:),DifChg60,xtick,ytick,false,N,density);
         
         avg=avg';
    currfig=find_figure('different_mapChg-60')

    avg2=interp2(yi1,xi1,avg,yii,xii);
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yi,xi,avg); axis image;axis ij;
    colormap(jet)
    colorbar
    currfig=find_figure('different_mapChg40')
    avg = gridavg2(pos(1,:),pos(2,:),DifPK40,xtick,ytick,false,N,density);
    avg=avg';     
    avg2=interp2(yi1,xi1,avg,yii,xii);
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yi,xi,avg); axis image;axis ij;
    colormap(jet)
    colorbar