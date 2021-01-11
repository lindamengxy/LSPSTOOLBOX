function [cells,handles,flag] = readdata_Gui(handles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cells=[];
%close all
%% data properties
flag=0;

if strcmp(handles.experimenttype,'Whole cell')
    
    if ~isfield(handles,'filepath')||isempty(handles.filepath)
        if ~isfield(handles,'data_folder')||isempty(handles.data_folder)||sum(ishandle(handles.data_folder))
            set(handles.error_run,'string','Data folder is empty. Please input the name of the folder where you saved all experimental data and rerun')
            return;
            
        end
        if ~isfield(handles,'mapname')||isempty(handles.mapname)
            set(handles.error_run,'string','Mapname is empty. Please fill in the mapname and rerun')
            return;
        end
        cell_folder=handles.data_folder;
        mapname=handles.mapname;
        %experiment=handles.experiment;
        path_file=fullfile(cell_folder,mapname);
        
        
    else
        path_file=handles.filepath;
        
        
    end
    
    [pth fname ext] = fileparts(path_file);
    handles.mapname=fname;
    mapname=fname;
    handles.mapfolder=pth;
    cell_folder=pth;
    set(handles.mapname1,'string',handles.mapname)
    if ~isfield(handles, 'laserstimValue')||isempty(handles.laserstimValue)
        handles.laserstim=1;
    end
    
    if ~isfield(handles,'eventwindowValue')||isempty(handles.eventwindowValue)
        handles.eventwindow=50;
    end
    if ~isfield(handles,'thresholdValue')||isempty(handles.thresholdValue)
        handles.threshold=10;
    end
    
    
%     if ~isfield(handles,'Tp')|isempty(handles.Tp)
%         handles.Tp=1;
%     end
    
    
    if ~isfield(handles,'age')||isempty(handles.age)
        set(handles.error_run,'string','Age is empty. Please fill in the mouse age and rerun')
        return;
    else
        
        set(handles.error_run,'string','')
    end
    if ~isfield(handles,'flipimgValue')||isempty(handles.flipimgValue)
        handles.flipimg=0;
    end
    if ~isfield(handles,'flipimgValue2')||isempty(handles.flipimgValue2)
        handles.flipimg2=0;
    end
    if ~isfield(handles,'darkexpValue')||isempty(handles.darkexpValue)
        handles.darkexp=0;
        
    end
      if ~isfield(handles,'BDValue')||isempty(handles.BDValue)
        handles.BD=5;
    end
    if ~isfield(handles,'save_path')||isempty(handles.save_path)||sum(ishandle(handles.save_path))
        
        set(handles.error_run,'string','Save path is empty. Please fill in the path way that you want to restore your data and rerun it')
        return;
    else
        set(handles.error_run,'string','')
    end
    
    if ~isfield(handles,'experimenttypeValue')||isempty(handles.experimenttypeValue)
        set(handles.error_run,'string','Experiment type is empty. Please select one type for your experiment and rerun it')
        return;
    else
        set(handles.error_run,'string','')
    end
    if ~isfield(handles,'direct_tValue')||isempty(handles.direct_tValue)
        handles.direct_t=8;
    else
        set(handles.error_run,'string','')
    end
      if ~isfield(handles,'Tp')||isempty(handles.Tp)||sum(ishandle(handles.drug))
        set(handles.error_run,'string','drug type is empty. Please select one type for your experiment and rerun it')
        return;
    else
        set(handles.error_run,'string','')
    end
    
  
foldername=sprintf('%s%i','eventwindow',handles.eventwindow);
save_folder=fullfile(handles.save_path,foldername);
handles.save_folder=save_folder;
%set(handles.save_path,'string',save_folder)
    
    Tp=handles.Tp;
    LaserStim=handles.laserstim;
    eventwindow=handles.eventwindow;
    flipimg=handles.flipimg;
    flipimg2=handles.flipimg2;
    darkexp=handles.darkexp;
    age=handles.age;
    % save data path
    path_save=handles.save_folder;
    
    direct_t1=handles.direct_t;
    direct_t2=direct_t1;
    %% read data
    flag=1;
    
    header = load(fullfile(pth, [fname '.mat']));
    
    fid = fopen(fullfile(pth, [fname '.dat']),'r');
    nTraces = fread(fid,1,'uint32');
    nSamples = fread(fid,1,'uint32');
    [data cnt] = fread(fid,'float64');
    fclose(fid);
    nPts = cnt./nTraces./nSamples;
    data1 = reshape(data,[nSamples nTraces nPts]);
    % traces: (1) amplifier output (2) photodiode output
    data = squeeze(data1(:,1,:));
    header.nTraces = nTraces;
    header.nSamples = nSamples;
    header.nPts = nPts;
    
    
    header.stimStartSamp = find(diff(header.digitalOutput),1)+1;
    %     if nTraces==3
    %         LFP=squeeze(data1(:,3,:));
    %         header.LFP=LFP;
    %         LFPfilter=LSPSfilter_LFP(header,LFP);
    %
    %     else
    %         header.LFP=NaN;
    %
    %     end
    
    if strcmp(mapname,'Map_111130_1954') header.holdingPotential=50;
    end
    if strcmp(mapname,'Map_111208_1504') header.holdingPotential=50;
    end
    if strcmp(mapname,'Map_120327_1501') header.holdingPotential=-60;
    end
    if strcmp(mapname,'Map_120403_1901') header.holdingPotential=10;
    end
    if strcmp(mapname,'Map_120510_1401') header.holdingPotential=10;
    end
    if strcmp(mapname,'Map_120512_1445') header.holdingPotential=-60;
    end
    if strcmp(mapname,'Map_120118_2306') header.holdingPotential=50;
    end
    if strcmp(mapname,'Map_120118_1952') header.holdingPotential=50;
    end
    if strcmp(mapname,'Map_120306_1600') header.holdingPotential=50;
    end
    if strcmp(mapname,'Map_120314_2040') header.holdingPotential=50;
    end
    
    
    
    thr60=handles.threshold;thr10=handles.threshold;
    %% filter data
    fdata = LSPSfilter(header,data,500);
    cells.mapname=mapname;
    cells.cell_folder=cell_folder;
    %% smooth data
    % fdata=smooth(fdata,11);
    for i=1:1:size(fdata,2)
        fdata(:,i)=smooth2(fdata(:,i),7);
    end
    % if laser is off
    if ~LaserStim
        [fbaseline fbaseline_rng fbaseline_fit] = LSPSgetbaseline2_spont(header,fdata,thr60,thr10);% get spontaneous baseline (the mean of the medians of small-window data set)
        meanb=mean(fbaseline(:));
        
        % judge whether the holding potential is recorded correctly
        if meanb<0&header.holdingPotential>0
           if ~isfield(handles,'RectifyHoldingPotential1')|isempty(handles.RectifyHoldingPotential1)
                set(handles.error_run,'string',sprintf('Current holding potential is %i. Please make sure it is right. If not, input the correct value in the upper blanket rectify holding potential and rerun it',header.holdingPotential))
                return;
            else   header.holdingPotential=handles.RectifyHoldingPotential1; set(handles.error_run,'string',''); handles.RectifyHoldingPotential1=[];
            end
        end
        
        if meanb>0&header.holdingPotential<0
             if ~isfield(handles,'RectifyHoldingPotential1')|isempty(handles.RectifyHoldingPotential1)
                set(handles.error_run,'string',sprintf('Current holding potential is %i. Please make sure it is right. If not, input the correct value in the upper blanket rectify holding potential and rerun it',header.holdingPotential))
                return;
            else   header.holdingPotential=handles.RectifyHoldingPotential1; set(handles.error_run,'string',''); handles.RectifyHoldingPotential1=[];
            end
        end
        
        [events fiterror header]= LSPSgetevents3_spont(header,fdata,fbaseline,fbaseline_rng,data);% get spontaneous events
        Boundry=boundary(header,flipimg,handles.BD,flipimg2);
        EventPeak=[];
        for i=1:header.nPts
            %R=R+sum(events.flag{i}(:)>0);
            EventPeak=[EventPeak events.peakAmp{i}(find(~isnan(events.flag{i}(:))&events.flag{i}(:)>0))];
        end
        xx=0:2:100;
        SpontPeak=[];R=0;
        for i=1:header.nPts
            R=R+sum(events.flag{i}(:)>0);
            SpontPeak=[SpontPeak events.peakAmp{i}(find(~isnan(events.peakAmp{i}(:))&events.peakAmp{i}(:)>0))];
        end
        
        %
        cells.Boundary=Boundry;
        cells.SpontPeak=SpontPeak;
        cells.eventsrate=R;
        cells.events=events;
        cells.header=header;
        cells.baseline=fbaseline;
        cells.experimenttype=handles.experimenttype;
        
        %cells.experiment=experiment;
        
        %cells.subdir=cell;
        cells.LaserStim=LaserStim;
        
        cells.flipimg=flipimg;
         cells.flipimg2=flipimg2;
        cells.age=age;
        cells.data=fdata;
        cells.darkexp=darkexp;
        cells.Tp=Tp;
        
         cells.eventwindow=handles.eventwindow;
        cells.direct_t=handles.direct_t
        foname=handles.mapfolder;
        foldername=fullfile(foname);
        if (exist(foldername) == 0)
            mkdir (foldername);
        end
%         f1=fullfile(foldername,experiment);
%         if (exist(f1) == 0)
%             mkdir (f1);
%         end
%         f2=fullfile(f1,cell);
%         if (exist(f2) == 0)
%             mkdir (f2);
%         end
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
        f4=fullfile(foldername,typefolder);
        if (exist(f4) == 0)
            mkdir (f4);
        end
        
        filesname1=fullfile(f4,sprintf('hist_spontpeak%i.png',header.holdingPotential));
        filesname2=fullfile(f4,sprintf('hist_spontpeak%i.fig',header.holdingPotential));
        filesname3=fullfile(f4,sprintf('hist_spontpeak%i.eps',header.holdingPotential));
        saveImages=1;
        if saveImages == 1
            xx=0:2:200;
            
            find_figure('SpontaneousPeakDistribution')
            yy=hist(SpontPeak,xx);
            yy=yy./length(SpontPeak);
            bar(xx,yy,'k')
            title(sprintf('spontaneous peak, rate=%i, holdingpotential=%i',R/header.nPts,header.holdingPotential))
            hold on
            yy=[];
            for j=1:1:length(xx)
                yy=[yy;sum(SpontPeak<xx(j))./length(SpontPeak)];
            end
            plot(xx,yy,'r--')
            j=2:length(yy);
            k=find(yy(j-1)<0.9&yy(j)>=0.9);
            
            line([xx(k),xx(k)],[0,1],'Color','g','LineStyle','--')
            text(20,0.5,sprintf('x=%i',xx(k)))
            hold off
            saveas(gcf,filesname1,'png');
            saveas(gcf,filesname2,'fig');
            print(gcf, '-depsc2', filesname3);
            %         print(gcf, '-dpng', sprintf('peak_%i.png',header.holdingPotential));
        end
        
    else
        if Tp==2
            [fbaseline fbaseline_rng fbaseline_fit] = LSPSgetbaseline2_PTX(header,fdata,thr60,thr10);% get spontaneous baseline (the mean of the medians of small-window data set)
            meanb=mean(fbaseline(:));
        else
            [fbaseline fbaseline_rng fbaseline_fit] = LSPSgetbaseline2(header,fdata,thr60,thr10);% get spontaneous baseline (the mean of the medians of small-window data set)
            meanb=mean(fbaseline(:));
        end
        % judge whether the holding potential is recorded correctly
       if meanb<0&header.holdingPotential>0
            if ~isfield(handles,'RectifyHoldingPotential1')||isempty(handles.RectifyHoldingPotential1)
                set(handles.error_run,'string',sprintf('Current holding potential is %i. Please make sure it is right. If not, input the correct value in the upper blanket rectify holding potential and rerun it',header.holdingPotential))
                return;
            else   header.holdingPotential=handles.RectifyHoldingPotential1; set(handles.error_run,'string',''); handles.RectifyHoldingPotential1=[];
            end
        end
        
        if meanb>0&header.holdingPotential<0
             if ~isfield(handles,'RectifyHoldingPotential1')||isempty(handles.RectifyHoldingPotential1)
                set(handles.error_run,'string',sprintf('Current holding potential is %i. Please make sure it is right. If not, input the correct value in the upper blanket rectify holding potential and rerun it',header.holdingPotential))
                return;
            else   header.holdingPotential=handles.RectifyHoldingPotential1; set(handles.error_run,'string',''); handles.RectifyHoldingPotential1=[];
            end
        end
        
        [events fiterror header]= LSPSgetevents3(header,fdata,fbaseline,fbaseline_rng,direct_t1,direct_t2,eventwindow);% get spontaneous events
        Boundry=boundary(header,flipimg,handles.BD,flipimg2);
        EventPeak=[];
        for i=1:header.nPts
            %R=R+sum(events.flag{i}(:)>0);
             hh=find(~isnan(events.flag{i}(:))&events.flag{i}(:)>1);
            if length(hh)
            EventPeak=[EventPeak events.peakAmp{i}(hh(1))];
            else  EventPeak=[EventPeak 0];
            end
        end
    
        cells.Boundary=Boundry;
        cells.EventPeak=EventPeak;
        cells.events=events;
        cells.header=header;
        cells.baseline=fbaseline;
        %cells.experiment=experiment;
        %cells.subdir=cell;
        cells.LaserStim=LaserStim;
        cells.mapname=mapname;
        cells.flipimg=flipimg;
         cells.flipimg2=flipimg2;
        cells.age=age;
        cells.data=fdata;
        cells.darkexp=darkexp;
        cells.Tp=Tp;
        cells.cell_folder=cell_folder;
         cells.eventwindow=handles.eventwindow;
        cells.direct_t=handles.direct_t
         cells.experimenttype=handles.experimenttype;
        if (exist(path_save) == 0)
            mkdir (path_save);
        end
        datafile=sprintf('LSPS_%s.mat',mapname);
        
        datafile=fullfile(path_save,datafile)
        if abs(header.nPts-length(header.StimCoordinates(1,:)))<0.1
            save(datafile,'cells','-mat')
            set(handles.error_run,'string','')
            return;
        else
            set(handles.error_run,'string','Map is not finished')
            return;
        end
        %      foldername=fullfile(foname);
        %     if (exist(foldername) == 0)
        %         mkdir (foldername);
        %     end
        %     f1=fullfile(foldername,experiment);
        %     if (exist(f1) == 0)
        %         mkdir (f1);
        %     end
        %     f2=fullfile(f1,cell);
        %     if (exist(f2) == 0)
        %         mkdir (f2);
        %     end
        %     if Tp==2
        %         typefolder='PTX';
        %     elseif Tp==3
        %         typefolder='APV_PTX';
        %     elseif Tp==1
        %         typefolder='highMg';
        %     else typefolder='internalMK801';
        %     end
        %     f4=fullfile(f2,typefolder);
        %     if (exist(f4) == 0)
        %         mkdir (f4);
        %     end
        %
        %     filesname1=fullfile(f4,sprintf('hist_eventpeak%i.png',header.holdingPotential));
        %     filesname2=fullfile(f4,sprintf('hist_eventpeak%i.fig',header.holdingPotential));
        %     filesname3=fullfile(f4,sprintf('hist_eventpeak%i.eps',header.holdingPotential));
        %     saveImages=1;
        %     if saveImages == 1
        %         saveas(gcf,filesname1,'png');
        %         saveas(gcf,filesname2,'fig');
        %         print(gcf, '-depsc2', filesname3);
        %         %         print(gcf, '-dpng', sprintf('peak_%i.png',header.holdingPotential));
        %     end
    end
else
     if ~isfield(handles,'filepath')||isempty(handles.filepath)
        if ~isfield(handles,'data_folder')||isempty(handles.data_folder)||sum(ishandle(handles.data_folder))
            set(handles.error_run,'string','Data folder is empty. Please input the name of the folder where you saved all experimental data and rerun')
            return;
            
        end
        if ~isfield(handles,'mapname')||isempty(handles.mapname)
            set(handles.error_run,'string','Mapname is empty. Please fill in the mapname and rerun')
            return;
        end
        cell_folder=handles.data_folder;
        mapname=handles.mapname;
        %experiment=handles.experiment;
        path_file=fullfile(cell_folder,mapname);
        
        
    else
        path_file=handles.filepath;
        
        
    end
    
    [pth fname ext] = fileparts(path_file);
    handles.mapname=fname;
    mapname=fname;
    handles.mapfolder=pth;
    cell_folder=pth;
    set(handles.mapname1,'string',handles.mapname)
    if ~isfield(handles, 'laserstimValue')||isempty(handles.laserstimValue)
        handles.laserstim=1;
    end
    
    if ~isfield(handles,'eventwindowValue')||isempty(handles.eventwindowValue)
        handles.eventwindow=50;
    end
    if ~isfield(handles,'thresholdValue')||isempty(handles.thresholdValue)
        handles.threshold=-20;
    end
    
    
%     if ~isfield(handles,'Tp')|isempty(handles.Tp)
%         handles.Tp=1;
%     end
    
    
    if ~isfield(handles,'age')||isempty(handles.age)
        set(handles.error_run,'string','Age is empty. Please fill in the mouse age and rerun')
        return;
    else
        set(handles.error_run,'string','')
    end
    if ~isfield(handles,'flipimg')||isempty(handles.flipimg)||ishandle(handles.flipimg)
        handles.flipimg=0;
    end
    if ~isfield(handles,'flipimg2')||isempty(handles.flipimg)||ishandle(handles.flipimg2)
        handles.flipimg2=0;
    end
    if ~isfield(handles,'darkexpValue')||isempty(handles.darkexpValue)
        handles.darkexp=0;
        
    end
    if ~isfield(handles,'save_path')||isempty(handles.save_path)||sum(ishandle(handles.save_path))
        
        set(handles.error_run,'string','Save path is empty. Please fill in the path way that you want to restore your data and rerun it')
        return;
    else
        set(handles.error_run,'string','')
    end
    
    if ~isfield(handles,'experimenttypeValue')||isempty(handles.experimenttypeValue)
        set(handles.error_run,'string','Experiment type is empty. Please select one type for your experiment and rerun it')
        return;
    else
        set(handles.error_run,'string','')
    end
    if ~isfield(handles,'direct_tValue')||isempty(handles.direct_tValue)
        handles.direct_t=8;
    else
        set(handles.error_run,'string','')
    end
      if ~isfield(handles,'Tp')||isempty(handles.Tp)||sum(ishandle(handles.drug))
        set(handles.error_run,'string','drug type is empty. Please select one type for your experiment and rerun it')
        return;
    else
        set(handles.error_run,'string','')
    end
    
    
     
foldername=sprintf('%s%i','eventwindow',handles.eventwindow);
save_folder=fullfile(handles.save_path,foldername);
handles.save_folder=save_folder;
%set(handles.save_path,'string',save_folder)
    %experiment=handles.experiment;
    
    
    %cell=handles.cell_folder;
    
    LaserStim=handles.laserstim;
    
    flipimg=handles.flipimg;
    flipimg2=handles.flipimg2;
    darkexp=handles.darkexp;
    mapname=handles.mapname;
    age=handles.age;
%     mappath=fullfile(experiment,cell,mapname);
    % read data path
%     path=handles.data_folder;
    % save data path
    path_save=handles.save_folder;
    foname=cell_folder;
    flag=1;
    %% read data
%     path_file=fullfile(path,mappath);
%     [pth fname ext] = fileparts(path_file);
    header = load(fullfile(pth, [fname '.mat']));
    fid = fopen(fullfile(pth, [fname '.dat']),'r');
    nTraces = fread(fid,1,'uint32');
    nSamples = fread(fid,1,'uint32');
    [data cnt] = fread(fid,'float64');
    fclose(fid);
    nPts = cnt./nTraces./nSamples;
    data1 = reshape(data,[nSamples nTraces nPts]);
    % traces: (1) amplifier output (2) photodiode output
    data = squeeze(data1(:,1,:));
    header.nTraces = nTraces;
    header.nSamples = nSamples;
    header.nPts = nPts;
    
    
    header.stimStartSamp = find(diff(header.digitalOutput),1)+1;
%     if nTraces==3
%         LFP=squeeze(data1(:,3,:));
%         header.LFP=LFP;
%         LFPfilter=LSPSfilter_LFP(header,LFP);
%         
%     else
%         header.LFP=NaN;
%         
%     end
    
    
    
    
    
    recordingdur=fix(header.nSamples/header.sampleRate)*1000;
    fdata = LSPSfilter(header,data,1000);
    % LFPfilter=LSPSfilter_LFP(header,LFP);
    
    %stimulus start time
    st=find(header.stimulatorOutput(4,:));
    for i=1:1:size(fdata,2)
        fdata(:,i)=smooth2(fdata(:,i),7);
    end
    Latency=zeros(size(fdata,2),1);
    %get baseline
    EventPeak=zeros(size(fdata,2),1);
    thr1=handles.threshold;thr0=handles.threshold;
    [fbaseline fbaseline_rng fbaseline_fit] = LSPSgetbaseline2(header,fdata,thr1,thr0);
    meanb=mean(fbaseline(:));
    for k=1:size(fdata,2)
        
        data_=(fdata(:,k)-fbaseline(:,k));
        negtive_event=sign(data_)<0;
        data_=abs(data_);
        
        pot_spikes=find(data_(st(1):end)>fbaseline_rng(k)&negtive_event(st(1):end));
        if ~isempty(pot_spikes)
            Latency(k)=(pot_spikes(1)).*recordingdur./header.sampleRate;
            EventPeak(k)=max(data_(pot_spikes));
        else Latency(k)=NaN;EventPeak(k)=NaN;
        end
    end
    Tp=handles.Tp;
    cells.Tp=Tp;
    cells.Latency=Latency;
    cells.EventPeak=EventPeak;
    cells.cell_folder=cell_folder;
    cells.header=header;
    cells.baseline=fbaseline;
   % cells.experiment=experiment;
   % cells.subdir=cell;
    cells.LaserStim=LaserStim;
    cells.mapname=mapname;
    cells.flipimg=flipimg;
    cells.flipimg2=flipimg2;
    cells.age=age;
    cells.data=fdata;
    cells.darkexp=darkexp;
    cells.eventwindow=handles.eventwindow;
    cells.direct_t=handles.direct_t;
    cells.experimenttype=handles.experimenttype;
 if (exist(path_save) == 0)
        mkdir (path_save);
 end
    
datafile=sprintf('LSPSCellattached_%s.mat',mapname);

datafile=fullfile(path_save,datafile)

if abs(header.nPts-length(header.StimCoordinates(1,:)))<0.1
    save(datafile,'cells','-mat')
    set(handles.error_run,'string','')
    return;
else
    set(handles.error_run,'string','Map is not finished')
    return;
end
   
    
end




% draw maps



end

