
% expect that this function is called with the filtered data
function [events bestfit_error hder] = LSPSgetevents3(header,data,baseline,baseline_rng,direct_t1,direct_t2,eventwindow)
% Detect the events, ignoring the event less than 2 ms; If the distance from the two events is
% less than 5 ms, they will be connected and supposed to be one events.
% xxx - expose parameters?
%Flag==1.....events<10
%Flag==2.....events<50
%Flag==3.....50<events<500


% get soma

if length(header.Soma1Coordinates)
    somax=header.Soma1Coordinates(1,1);
    somay=header.Soma1Coordinates(1,2);


else
   
   currfig=find_figure('ClickCell')
    [imd imdIdx] = sort(header.ImageData(:));
    imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
    imd = uint8(reshape(imd, size(header.ImageData)));
    
    
    ImageScale=header.ImageScale;
    
    
    image([-1 1].*ImageScale(1)/2, [-1 1].*ImageScale(2)/2, repmat(imd,[1 1 3]));
    hold on
    title('click the cell you patch')
    axis xy;
    set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, ...
        'DataAspectRatio',[ 1 1 1]);
    [x,y]=ginput(1);
    plot(x,y,'rO')
    hold off
    header.Soma1Coordinates=zeros(1,2);
    header.Soma1Coordinates(1,1)=x;
    header.Soma1Coordinates(1,2)=y;
    somax=x;
    somay=y;
    close(currfig)
end



recordingdur=fix(header.nSamples/header.sampleRate)*1000;

bestfit_error=[];

maxGapSamps = fix(0.002*header.sampleRate); % sample rate is 10000; 0.002*10000=2 ms
minEventSamps = fix(0.002*header.sampleRate);% 2ms.  
aldur=10;

% stimulus start
stim_start=find(header.stimulatorOutput(4,:)); 
% shutter
wsize=length(1:1:stim_start(1));
wrng = fix([0.01 0.99]*wsize);
% select for events as areas that are outside of the baseline range
x = data-baseline;
trace1=x(1:10:end,:);
%Flag=ploteventtraces(trace1,header);
negative_event = (sign(x)<0); 

if abs(header.holdingPotential-10)<20
    x=x.*(~negative_event);
end
x = abs(x);

% use 2 times std of baseline to be the cutoff boundary.

sel = bsxfun(@gt, x, baseline_rng);
Eve1=zeros(size(data));
Eve2=zeros(size(data));
% label continuous events that occur after stimulation
events.nEvents = zeros(1,header.nPts);
% events.negflag=cell(1,header.nPts);
% events.posflag=cell(1,header.nPts);
events.startSamp = cell(1,header.nPts);
events.stopSamp = cell(1,header.nPts);
events.peakSamp = cell(1,header.nPts);
events.peakAmp = cell(1,header.nPts);
events.negative = cell(1,header.nPts);
events.area = cell(1,header.nPts);
events.duration = cell(1,header.nPts);
events.direct=cell(1,header.nPts);
events.data=cell(1,header.nPts);
events.eve1=cell(1,1);
events.eve2=cell(1,1);
events.tau=cell(1,header.nPts);
events.Numevents=0;
events.dpri=cell(1,header.nPts);
events.flag=cell(1,header.nPts);
events.flag_llatency=cell(1,header.nPts);
events.rng=zeros(1,header.nPts);
baseline_rng2=zeros(1,header.nPts);

for i=1:header.nPts
 
    %z score
    events.dpri{i}(1)=(mean(x(stim_start(1):1:stim_start(1)+eventwindow*fix(header.sampleRate/recordingdur),i))-mean(x(5:1:stim_start(1)-1,i)))/std(x(5:1:stim_start(1)-1,i))*sqrt(eventwindow*fix(header.sampleRate/recordingdur));
      % if events happens before l
 % plot(x(1:10:5000,i))
  %trace bound 99% of baseline;
x_sortbaseline=sort(x(1:1:stim_start(1),i));

rng=max([x_sortbaseline(wrng(1)),x_sortbaseline(wrng(2))]);

%rng=mean(baseline_rng);
events.rng(i)=rng;
if header.holdingPotential>20
    direct_dt=direct_t2;
    baseline_rng1=10;
else baseline_rng1=5;
    direct_dt=direct_t1;
end
baseline_rng1=max(rng,baseline_rng1);
baseline_rng1=min(2*rng,baseline_rng1);
baseline_rng2(i)=baseline_rng1; 
      
  cur = sel(:,i);
  cur(1:header.stimStartSamp)=false;

  % fill in any small holes between same-type events
  cur1 = LSPSfillgaps(~(cur&negative_event(:,i)),maxGapSamps);
  cur2 = LSPSfillgaps(~(cur&~negative_event(:,i)),maxGapSamps);
  cur = (~cur1|~cur2);
  
  % remove very short duration events
  cur = bwareaopen(cur,minEventSamps);
  
  % identify events
  CC = bwconncomp(cur);
  events.nEvents(i) = length(CC.PixelIdxList);
%   if events.dpri{i}(1)>=50&&~events.nEvents(i)
%       sel(:,i)=bsxfun(@gt,x(:,i),baseline_rng1);
%       cur = sel(:,i);
%       cur(1:header.stimStartSamp)=false;
%       
%       % fill in any small holes between same-type events
%       cur1 = LSPSfillgaps(~(cur&negative_event(:,i)),maxGapSamps);
%       cur2 = LSPSfillgaps(~(cur&~negative_event(:,i)),maxGapSamps);
%       cur = (~cur1|~cur2);
%       
%       % remove very short duration events
%       cur = bwareaopen(cur,minEventSamps);
%       
%       % identify events
%       CC = bwconncomp(cur);
%       events.nEvents(i) = length(CC.PixelIdxList);
%       
%   end

%  std_b=std(x(1:1:stim_start(1),i));
% mean_b=0;
 
  %events.dpri(i)=(-median(x(stim_start(1)-599:1:stim_start(1),i))+median(x(stim_start(1)+1:1:stim_start(1)+600,i)));
  if ~events.nEvents(i)
%       events.negflag{i}=0;
%       events.posflag{i}=0;
       events.flag{i}=0;
       events.peakAmp{i}=0;
       events.tau{i}=NaN;
       events.startSamp{i} =0;  
       events.stopSamp{i} =0;
       events.peakSamp{i} = 0; 
       events.negative{i} =0;
       events.area{i} =0;
       events.duration{i} =0; 
       events.tau{i}=0;
       events.flag_llatency{i}=0;
        
       %events.dpri{i}(1)=mean(x(stim_start(1):1:stim_start(1)+60*fix(header.sampleRate/recordingdur),i))-mean(x(stim_start(1)-length(stim_start(1):1:stim_start(1)+60*fix(header.sampleRate/recordingdur)):1:stim_start(1)-1,i))/sqrt(std(x(stim_start(1):1:stim_start(1)+60*fix(header.sampleRate/recordingdur),i))+std(x(stim_start(1)-length(stim_start(1):1:stim_start(1)+60*fix(header.sampleRate/recordingdur)):1:stim_start(1)-1,i)));
  else
   
      %plot trace
%       find_figure('single trace');
%       clf;
%       plot(baseline_rng(i).*ones(size(x(1:10:10000,i))))
%       hold on
%       plot(events.rng(i).*ones(size(x(1:10:10000,i))),'r--')
%   
%      
%       plot(x(1:10:10000,i),'g')       
      
        
      
     stdard=0; 
%   events.negflag{i}=ones(1,events.nEvents(i));
%   events.posflag{i}=ones(1,events.nEvents(i));
     events.flag{i}=zeros(1,events.nEvents(i));
     events.flag_llatency{i}=zeros(1,events.nEvents(i));
     events.peakAmp{i} = zeros(1,events.nEvents(i));
     events.startSamp{i} = zeros(1,events.nEvents(i));  
     events.stopSamp{i} = zeros(1,events.nEvents(i));
     events.peakSamp{i} = zeros(1,events.nEvents(i)); 
     events.negative{i} = zeros(1,events.nEvents(i));
     events.area{i} = zeros(1,events.nEvents(i));
     events.duration{i} = zeros(1,events.nEvents(i)); 
     events.tau{i}=zeros(1,events.nEvents(i));     
     events.dpri{i}=(0:1:events.nEvents(i)-1).*stdard;
     
     events.rise{i} = zeros(1,events.nEvents(i));
      meanbaseline=abs(mean(data(stim_start(1)-20:stim_start(1),i))-mean(data(1:stim_start(1),i)));
      %z score
    events.dpri{i}(1)=(mean(x(stim_start(1):1:stim_start(1)+eventwindow*fix(header.sampleRate/recordingdur),i))-mean(x(5:1:stim_start(1)-1,i)))/std(x(5:1:stim_start(1)-1,i))*sqrt(eventwindow*fix(header.sampleRate/recordingdur));
      % if events happens before laser stimuli time, it will treated not
      % to be stimuli driven
     if meanbaseline>10

         events.flag{i}(:)=0;
         events.flag_llatency{i}(:)=0;
         events.dpri{i}(:)=0;
         for j=1:events.nEvents(i)
          events.startSamp{i}(j) = NaN;
           events.stopSamp{i}(j) = NaN;     
          events.peakAmp{i}(j) = NaN;
          events.peakSamp{i}(j) =NaN;
          events.negative{i}(j) = NaN;        
          events.duration{i}(j) = NaN;          
          events.area{i}(j) = NaN;
          events.rise{i}(j)=NaN;
         end
     else
         times=1;    
         Ind_pre_end=stim_start(1);
         %z-score
       
         
         for j=1:events.nEvents(i) 
             
             
             mincc=min(CC.PixelIdxList{j});
             maxcc=max(CC.PixelIdxList{j});
             if mincc<stim_start(1)%if events begin before laser onset, set everything to NaN
                 flagindrange=0;
                 events.startSamp{i}(j)=NaN;
                 events.stopSamp{i}(j) = NaN;     
                 events.peakAmp{i}(j)= NaN;
                 events.peakSamp{i}(j) =NaN;
                 events.negative{i}(j) = NaN;        
                 events.duration{i}(j) = NaN;          
                 events.area{i}(j) = NaN;
                 events.tau{i}(j)=NaN;
                 events.rise{i}(j)=NaN;
             else
                 
               Ind1=find2stdindex(x,mincc,rng,maxGapSamps,minEventSamps,-1,i,header);% trace back to 99% range of baseline
               Ind2=find2stdindex(x,maxcc,rng,maxGapSamps,minEventSamps,1,i,header);
%                Ind1=mincc;
%                Ind2=maxcc;
               indrange1=Ind1:1:Ind2;
               %integration window size is 200ms
               if Ind2*recordingdur/header.sampleRate>200+stim_start(1)*recordingdur/header.sampleRate
                   Ind2=200*fix(header.sampleRate/recordingdur)+stim_start(1);
               end
               
               if Ind1>=Ind_pre_end&Ind1<200*fix(header.sampleRate/recordingdur)+stim_start(1)
                events.startSamp{i}(j) =Ind1;           
                events.stopSamp{i}(j) = Ind2; 
            
                indrange=Ind1:1:Ind2;
                [events.peakAmp{i}(j) k] = max(x(indrange,i));
                events.peakSamp{i}(j) = indrange(k);
                events.negative{i}(j) = negative_event(events.peakSamp{i}(j),i);
                events.area{i}(j) = sum(x(indrange,i))/header.sampleRate;
                events.duration{i}(j) = (events.stopSamp{i}(j)-events.startSamp{i}(j));
                events.rise{i}(j)=(indrange(k)-Ind1)*(recordingdur/header.sampleRate);
                flagindrange=1;%sign for meaningful responses
               else
                 indrange=NaN;
                 flagindrange=0;
                 events.startSamp{i}(j) = NaN;
                 events.stopSamp{i}(j) = NaN;     
                 events.peakAmp{i}(j) = NaN;
                 events.peakSamp{i}(j) = NaN;
                 events.negative{i}(j) =NaN;        
                 events.duration{i}(j) = NaN;          
                 events.area{i}(j) = NaN;
                 events.tau{i}(j)=NaN;
                  events.rise{i}(j)=NaN;
              end
               Ind_pre_end=Ind2;
            end
          
          % fit exponential
           if flagindrange==0  %  exclude all the events with latency larger than 200 ms or earlier than the stimuli
              events.flag{i}(j)=NaN;
              events.tau{i}(j)=NaN;
              events.dpri{i}(j)=NaN;
              events.flag_llatency{i}(j)=NaN;
               events.peakAmp{i}(j) = NaN;
                events.area{i}(j) =NaN;
               
           else
               
               % fit events wiht exponential 
              ind90=find(x(indrange1,i)>0.75*events.peakAmp{i}(j));
              if length(indrange1(ind90(end)):indrange1(end))>70   %if the decay time from 75% peak to baseline is less than 5ms, it will treat as noise not response       
              
                  D=x(indrange1(ind90(end)):indrange1(end),i);
                  D_norm=D./max(D);          
                  xx = (0:1:length(D_norm)-1).*recordingdur./header.sampleRate; yy = D_norm';
                  [xfit yerror]=lsqcurvefit(@(xfit,xdata) xfit(1).*exp(-xdata./xfit(2))+xfit(3),[1,2,0],xx,yy,[],[],optimset('display','off'));
                  vaf_r2 = 1 - yerror / sum((D_norm - mean(D_norm)).^2);
                  fit = (xfit(1).*exp(-xx./xfit(2))+xfit(3)).*max(D);   
              
              else vaf_r2=NaN;fit=[];
              end
              bestfit_error=[bestfit_error vaf_r2];
              fitstand=0.7; % standard for measure how good the fit is
        
              if header.holdingPotential>20
                  peakthreshold=20;
              else peakthreshold=20;
              end
                
             if ((vaf_r2<fitstand&&events.peakAmp{i}(j)<peakthreshold&&events.area{i}(j)<0.5)||(isnan(vaf_r2)&&events.peakAmp{i}(j)<peakthreshold&&events.area{i}(j)<0.5))% if tau=1ms, amp=10, y=amp*x/tau*exp(1-x/tau), area=0.027
                  events.flag{i}(j)=0;
                  events.flag_llatency{i}(j)=0;
                  events.tau{i}(j)=NaN;
                  events.startSamp{i}(j)=0;
                 % events.area{i}(j)=0;
              elseif events.startSamp{i}(j).*recordingdur./header.sampleRate<direct_dt+stim_start(1).*recordingdur./header.sampleRate&&events.startSamp{i}(j).*recordingdur./header.sampleRate>=stim_start(1).*recordingdur./header.sampleRate
                   events.flag{i}(j)=1;
                    events.flag_llatency{i}(j)=0;
                   if vaf_r2>fitstand
                    events.tau{i}(j)=xfit(2);
                   end
              elseif events.startSamp{i}(j).*recordingdur./header.sampleRate<eventwindow+stim_start(1).*recordingdur./header.sampleRate&&events.startSamp{i}(j).*recordingdur./header.sampleRate>=stim_start(1).*recordingdur./header.sampleRate+direct_dt
                  events.flag{i}(j)=2;
                   events.flag_llatency{i}(j)=0;
                  if vaf_r2>fitstand
                  events.tau{i}(j)=xfit(2);
                  end
              elseif events.startSamp{i}(j).*recordingdur./header.sampleRate<500+stim_start(1).*recordingdur./header.sampleRate&&events.startSamp{i}(j).*recordingdur./header.sampleRate>=stim_start(1).*recordingdur./header.sampleRate+eventwindow
                  events.flag{i}(j)=0;
                   events.flag_llatency{i}(j)=1;
                  
                   events.area{i}(j)=0;
                   events.tau{i}(j)=NaN;
%                   if vaf_r2>fitstand
%                   events.tau{i}(j)=xfit(2);
%                   end
              elseif events.startSamp{i}(j).*recordingdur./header.sampleRate>=500
                  events.flag{i}(j)=0;
                   events.flag_llatency{i}(j)=0;
                   events.tau{i}(j)=NaN;
%                   if vaf_r2>fitstand
%                   events.tau{i}(j)=xfit(2);
%                   end
              else events.flag{i}(j)=0;
                   events.flag_llatency{i}(j)=0;
                   %events.area{i}(j)=0;
                   events.tau{i}(j)=NaN;
%                  if vaf_r2>fitstand
%                   events.tau{i}(j)=xfit(2);
%                  end
              end
              
           end
        if ~isnan(events.negative{i}(j))
           if header.holdingPotential<0&&~events.negative{i}(j)
              events.flag{i}(j)=0; 
               events.flag_llatency{i}(j)=0;
               events.area{i}(j)=NaN;
                events.tau{i}(j)=NaN;
                events.startSamp{i}(j) = NaN;
                 events.stopSamp{i}(j) = NaN;     
                 events.peakAmp{i}(j) = NaN;
                 events.peakSamp{i}(j) = NaN;
                 events.negative{i}(j) =NaN;        
                 events.duration{i}(j) = NaN;   
                events.rise{i}(j)=NaN;
           end
           if header.holdingPotential>0&&events.negative{i}(j)
              events.flag_llatency{i}(j)=0;
              events.area{i}(j)=NaN;
               events.tau{i}(j)=NaN;
               events.startSamp{i}(j) = NaN;
                 events.stopSamp{i}(j) = NaN;     
                 events.peakAmp{i}(j) = NaN;
                 events.peakSamp{i}(j) = NaN;
                 events.negative{i}(j) =NaN;        
                 events.duration{i}(j) = NaN; 
                  events.rise{i}(j)=NaN;
                  events.flag{i}(j)=0; 
                  
           end  
           
        end

%           if events.flag{i}(j)>=1&&length(fit)
%             Eve1(indrange1,i)=x(indrange1,i);
%             Eve2(indrange1,i)=[x(indrange1(1:ind90(end)-1),i)' fit];
%          end 
        
         end % for j
            
    end       
      
%    if events.nEvents(i)&& length(find(events.flag{i}(:)==2))
%       events.Numevents=events.Numevents+1;
%    end
  end
end
%  events.eve1=Eve1;
%  events.eve2=Eve2;
hder=header;
% 
%   for i=1:1:header.nPts
%       find_figure('single trace')
%       plot(baseline_rng(i).*ones(size(x(1:10:5000,i))))
%       hold on
%       plot(baseline_rng2(i).*ones(size(x(1:10:5000,i))),'r--')
%       hold on
%       %plot(events.rng(i).*ones(size(x(1:10:5000,i))),'r--')
%       rectangle('Position',[stim_start(1).*recordingdur./header.sampleRate,-30,eventwindow,60],'LineStyle','--')
%        text(stim_start(1).*recordingdur./header.sampleRate, 20,num2str(events.dpri{i}(1)))
%        text(stim_start(1).*recordingdur./header.sampleRate, 30,num2str(events.area{i}(1)))
%       if length(find(events.flag{i}(:)>0))
%           
%         plot(x(1:10:5000,i),'g')       
%         hold on
%         plot(data(1:10:5000,i),'k')
%         plot(baseline(1:10:5000,i),'k')
%      
%       
%        
%         plot(events.eve2(1:10:5000,i),'r')
%         
% %       elseif length(find(events.flag_llatency{i}(:)>0))
% %           plot(x(1:10:5000,i),'y')  
% %            text(100,20,num2str(events.startSamp{i}(find(events.flag_llatency{i}(:)>0))))  
%           
%       else  plot(x(1:10:5000,i),'b')
%           
%          hold on
%          plot(data(1:10:5000,i),'k')
%          plot(baseline(1:10:5000,i),'k')
%       end
%       pause(2)
%       hold off
% %      
%   end
% %   

end % function
  

function sel = LSPSfillgaps(sel,maxGapSamps)

  CC1 = bwconncomp(sel);
  for j=1:length(CC1.PixelIdxList)
     if length(CC1.PixelIdxList{j}) <= maxGapSamps
        sel(CC1.PixelIdxList{j}) = false;
     end
  end
end

function Ind=find2stdindex(x,Indx,rng,maxGapSamps,minEventSamps,dir,i,header)
   % flag=0;
   Ind=0;
    if x(Indx,i)<=rng
        Ind=Indx;
    else
        %fill small gaps
        sel1 = bsxfun(@gt, x, rng);
        cur11 = sel1(:,i);
        cur11(1:header.stimStartSamp)=false;

       % fill in any small holes between two events
        cur11 = LSPSfillgaps(~cur11,maxGapSamps);
         cur11 = ~cur11;  
       % remove very short duration events
         %cur = bwareaopen(cur,minEventSamps);
  
  % identify events
       Cc = bwconncomp(cur11);  
        for jj=1:1:length(Cc.PixelIdxList)
           if Indx>=min(Cc.PixelIdxList{jj})&&Indx<=max(Cc.PixelIdxList{jj})
            if dir<0
            Ind=min(Cc.PixelIdxList{jj});
            else Ind=max(Cc.PixelIdxList{jj});
            end
            break;
           end
        end
    
        
            
%         if dir<0
%           for k=Indx-1:-1:stt
%             if x(k,i)<=times*stdb+mb&&x(k+1,i)>times*stdb+mb
%                 Ind=k;
%                 flag=1;
%                 break;
%             end
%             if flag==0
%                  Ind=stt;
%              end
%           end
%         else for k=Indx:1:length(x(:,i))-1
%                 if x(k,i)>=times*stdb+mb&&x(k+1,i)<times*stdb+mb
%                    Ind=k;
%                    flag=1;
%                    break;
%                 end
%             end
%              if flag==0
%                  Ind=length(x(:,i));
%              end
%         end
    end
end
                
%end



