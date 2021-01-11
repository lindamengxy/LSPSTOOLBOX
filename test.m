function varargout = test(varargin)
% TEST MATLAB code for test.fig
%      TEST, by itself, creates a new TEST or raises the existing
%      singleton*.
%
%      H = TEST returns the handle to a new TEST or the handle to
%      the existing singleton*.
%
%      TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST.M with the given input arguments.
%
%      TEST('Property','Value',...) creates a new TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test

% Last Modified by GUIDE v2.5 27-Apr-2015 15:41:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_OpeningFcn, ...
                   'gui_OutputFcn',  @test_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before test is made visible.
function test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test (see VARARGIN)
%set up a toggle button


% Choose default command line output for test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles)

% UIWAIT makes test wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% hb=uicontrol('style','togglebutton');
% set(hb,'Position',[520 6 120 50]);
% set(hb,'String','Run');
% set(hb,'callback',{@readdata_Gui,handles});
% guidata(hObject,handles)
% --- Outputs from this function are returned to the command line.
function varargout = test_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function data_folder_Callback(hObject, eventdata, handles)
% hObject    handle to data_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data_folder=get(hObject,'string')
guidata(hObject, handles);



% Hints: get(hObject,'String') returns contents of data_folder as text
%        str2double(get(hObject,'String')) returns contents of data_folder as a double


% --- Executes during object creation, after setting all properties.
function data_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function experiment_Callback(hObject, eventdata, handles)
% hObject    handle to experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.experiment=get(hObject,'string');
% Hints: get(hObject,'String') returns contents of experiment as text
%        str2double(get(hObject,'String')) returns contents of experiment as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in drug.
function drug_Callback(hObject, eventdata, handles)
% hObject    handle to drug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 contents = cellstr(get(hObject,'String'));
 handles.drug=contents{get(hObject,'Value')};
switch contents{get(hObject,'Value')}
    case {'High Mg'}

        handles.Tp=1;
    case {'High Mg + PTX'}
        handles.Tp=2;
    case{'High Mg + PTX +APV'}
        handles.Tp=3;
    case{'High Mg+TTX'}
        handles.Tp=4;
    otherwise
            handles.Tp=5;
end

guidata(hObject, handles);
        
% Hints: contents = cellstr(get(hObject,'String')) returns drug contents as cell array
%        contents{get(hObject,'Value')} returns selected item from drug


% --- Executes during object creation, after setting all properties.
function drug_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function laserstim_Callback(hObject, eventdata, handles)
% hObject    handle to panel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of panel1 as text
%        str2double(get(hObject,'String')) returns contents of panel1 as a double
handles.laserstim=str2num(get(hObject,'string'));
handles.laserstimValue=handles.laserstim;
    
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function panel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function mapname1_Callback(hObject, eventdata, handles)
% hObject    handle to mapname1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mapname=get(hObject,'string');

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of mapname1 as text
%        str2double(get(hObject,'String')) returns contents of mapname1 as a double


% --- Executes during object creation, after setting all properties.
function mapname1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mapname1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function direct_t_Callback(hObject, eventdata, handles)
% hObject    handle to direct_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of direct_t as text
%        str2double(get(hObject,'String')) returns contents of direct_t as a double
handles.direct_t=str2num(get(hObject,'string'));
handles.direct_tValue=handles.direct_t;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
% function direct_t_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to direct_t (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end






function eventwindow_Callback(hObject, eventdata, handles)
% hObject    handle to eventwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.eventwindow=str2num(get(hObject,'string'));
handles.eventwindowValue=handles.eventwindow;
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of eventwindow as text
%        str2double(get(hObject,'String')) returns contents of eventwindow as a double


% --- Executes during object creation, after setting all properties.
function eventwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over data_folder.



function save_folder_Callback(hObject, eventdata, handles)
% hObject    handle to save_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.save_path=get(hObject,'string')
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of save_folder as text
%        str2double(get(hObject,'String')) returns contents of save_folder as a double


% --- Executes during object creation, after setting all properties.
function save_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.


% Hint: place code in OpeningFcn to populate axes1


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 contents = cellstr(get(hObject,'String'));
 
 handles.datatype=contents{get(hObject,'Value')};
    
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot1.
function plot1_Callback(hObject,eventdata,handles)
% hObject    handle to plot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% panel1=handles.panel1;
eventwindow=handles.eventwindow;
darkexp=handles.darkexp;
mapname=handles.mapname;

% load data
% path_save=handles.save_path;
foname=handles.mapfolder;
% experiment=handles.experiment
% subdir=handles.cell_folder;
flipimg=handles.flipimg;
flipimg2=handles.flipimg2;
header=handles.cells.header;
fdata=handles.cells.data;
BD=handles.cells.Boundary;
age=handles.age;
plot_fig=1;
Tp=handles.Tp;
LaserStim=handles.laserstim;

if strcmp(handles.experimenttype,'Whole cell')
    direct_t1=handles.direct_t;
    direct_t2=handles.direct_t;
    events=handles.cells.events;
end
if LaserStim
      
      [Stim_coor]=Matrixrotate_singlemap(header,flipimg,flipimg2);
      X=Stim_coor(1,:);
      Y=Stim_coor(2,:);
      
 
      
       hold on;
      if strcmp(handles.datatype,'peak+DIC image')
      handles=plotsinglemap2_Gui(X,Y,events.peakAmp,mapname,foname,Tp,header,plot_fig,'peak+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles);
      elseif strcmp(handles.datatype,'charge')
      [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.area,mapname,foname,Tp,header,plot_fig,'charge',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles,BD);
      elseif strcmp(handles.datatype,'latency')
      [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.startSamp,mapname,foname,Tp,header,plot_fig,'Latency',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles,BD);
      elseif strcmp(handles.datatype,'peak')
      [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.peakAmp,mapname,foname,Tp,header,plot_fig,'peak',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles,BD);
      elseif strcmp(handles.datatype,'charge+DIC image')
      handles=plotsinglemap2_Gui(X,Y,events.area,mapname,foname,Tp,header,plot_fig,'charge+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles);
      elseif strcmp(handles.datatype,'Latency+DIC image')
      handles=plotsinglemap2_Gui(X,Y,events.startSamp,mapname,foname,Tp,header,plot_fig,'Latency+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles);
     
      elseif strcmp(handles.datatype,'traces+DIC image')
          handles=maptraces_Gui(fdata,header,mapname,foname,Tp,eventwindow,handles.axes_plot1,hObject,handles)
          
      elseif strcmp(handles.experimenttype,'Cell attached') && (strcmp(handles.datatype,'SpikeLatency+DIC image')|strcmp(handles.datatype,'SpikeLatency inside eventswindow'))
          SpikeLatencyImage_Gui(mapname,foname,Tp,header,handles.datatype,age,handles.cells,handles.axes_plot1,eventwindow);
      elseif strcmp(handles.experimenttype,'Cell attached') && strcmp(handles.datatype,'AcummulateDstributionSpikelatency')
          Spikeacummulation_Gui(mapname,foname,Tp,header,handles.datatype,age,handles.cells,handles.axes_plot1,eventwindow)
      elseif strcmp(handles.experimenttype,'Cell attached') && strcmp(handles.datatype,'AcummulateDistributionSpikelatencyEventwindow')
          Spikeacummulation_Gui(mapname,foname,Tp,header,handles.datatype,age,handles.cells,handles.axes_plot1,eventwindow)
      end
      
      guidata(hObject, handles);
     
      

end




% --- Executes on selection change in datatype2.
function datatype2_Callback(hObject, eventdata, handles)
% hObject    handle to datatype2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 contents = cellstr(get(hObject,'String'));
 
 handles.datatype=contents{get(hObject,'Value')}
    
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns datatype2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from datatype2


% --- Executes during object creation, after setting all properties.
function datatype2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datatype2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double
handles.threshold=str2num(get(hObject,'string'));
handles.thresholdValue=handles.threshold;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot_button2.
function plot_button2_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ispushed=get(hObject,'Value');
if ispushed
eventwindow=handles.eventwindow;
darkexp=handles.darkexp;
mapname=handles.mapname;

% load data
path_save=handles.save_path;
foname=handles.mapfolder;



% 
% experiment=handles.experiment
% subdir=handles.cell_folder;
flipimg=handles.flipimg;
header=handles.cells.header;
BD=handles.cells.Boundary;
LaserStim=handles.laserstim;
fdata=handles.cells.data;
age=handles.age;
plot_fig=1;
Tp=handles.Tp;
if strcmp(handles.experimenttype,'Whole cell')
    direct_t1=handles.direct_t;
    direct_t2=handles.direct_t;
    events=handles.cells.events;
end

if LaserStim
      
      [Stim_coor]=Matrixrotate_singlemap(header,flipimg);
      X=Stim_coor(1,:);
      Y=Stim_coor(2,:);
      
 
      
       
      if strcmp(handles.datatype,'peak+DIC image')
      handles=plotsinglemap2_Gui(X,Y,events.peakAmp,mapname,foname,Tp,header,plot_fig,'peak+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles);
      elseif strcmp(handles.datatype,'charge')
      [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.area,mapname,foname,Tp,header,plot_fig,'charge',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles,BD);
      elseif strcmp(handles.datatype,'latency')
      [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.startSamp,mapname,foname,Tp,header,plot_fig,'Latency',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles,BD);
      elseif strcmp(handles.datatype,'peak')
      [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.peakAmp,mapname,foname,Tp,header,plot_fig,'peak',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles,BD);
      elseif strcmp(handles.datatype,'charge+DIC image')
      handles=plotsinglemap2_Gui(X,Y,events.area,mapname,foname,Tp,header,plot_fig,'charge+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles);
      elseif strcmp(handles.datatype,'Latency+DIC image')
      handles=plotsinglemap2_Gui(X,Y,events.startSamp,mapname,foname,Tp,header,plot_fig,'Latency+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles);
     
      elseif strcmp(handles.datatype,'traces+DIC image')
          maptraces_Gui(fdata,header,mapname,foname,Tp,eventwindow,handles.axes_plot2,events)
      elseif strcmp(handles.experimenttype,'Cell attached') && (strcmp(handles.datatype,'SpikeLatency+DIC image')|strcmp(handles.datatype,'SpikeLatency inside eventswindow'))
          SpikeLatencyImage_Gui(mapname,foname,Tp,header,handles.datatype,age,handles.cells,handles.axes_plot2,eventwindow);
      elseif strcmp(handles.experimenttype,'Cell attached') && strcmp(handles.datatype,'AcummulateDstributionSpikelatency')
          Spikeacummulation_Gui(mapname,foname,Tp,header,handles.datatype,age,handles.cells,handles.axes_plot2,eventwindow)
      elseif strcmp(handles.experimenttype,'Cell attached') && strcmp(handles.datatype,'AcummulateDistributionSpikelatencyEventwindow')
          Spikeacummulation_Gui(mapname,foname,Tp,header,handles.datatype,age,handles.cells,handles.axes_plot2,eventwindow)
 
      end
      
end    
      %guidata(hObject,handles);
 set(hObject,'Value',0);     
      
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in run_analysis.
% function run_analysis_Callback(hObject, eventdata, handles)
% % hObject    handle to run_analysis (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% readata_Gui(handles)



function cell_folder_Callback(hObject, eventdata, handles)
% hObject    handle to cell_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cell_folder=get(hObject,'string')
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of cell_folder as text
%        str2double(get(hObject,'String')) returns contents of cell_folder as a double


% --- Executes during object creation, after setting all properties.
function cell_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Age_Callback(hObject, eventdata, handles)
% hObject    handle to Age (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Age=str2num(get(hObject,'string'));
handles.age=handles.Age;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Age as text
%        str2double(get(hObject,'String')) returns contents of Age as a double


% --- Executes during object creation, after setting all properties.
function Age_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Age (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function flipimg_Callback(hObject, eventdata, handles)
% hObject    handle to flipimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flipimg=str2num(get(hObject,'string'));
handles.flipimgValue=handles.flipimg;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of flipimg as text
%        str2double(get(hObject,'String')) returns contents of flipimg as a double


% --- Executes during object creation, after setting all properties.
function flipimg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flipimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function darkexp_Callback(hObject, eventdata, handles)
% hObject    handle to darkexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.darkexp=str2num(get(hObject,'string'))
handles.darkexpValue=handles.darkexp;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of darkexp as text
%        str2double(get(hObject,'String')) returns contents of darkexp as a double


% --- Executes during object creation, after setting all properties.
function darkexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darkexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in runanalysis.
% function runanalysis_Callback(hObject, eventdata, handles)
% % hObject    handle to runanalysis (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runanalysis


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ispushed=get(hObject,'Value');
if ispushed
    
   set(hObject, 'String', 'Calculating...')
   drawnow
%    guidata(hObject, handles);
    flag=0;
    
    [cells,handles,flag]=readdata_Gui(handles);
    if flag
    display('finished')
    set(hObject, 'String', 'Ready.');
    handles.cells=cells;
    handles.mapfolder=cells.cell_folder;
    end
    pause(2);
      drawnow update
      set(hObject, 'String', 'Run');
      set(hObject,'Value',0);
end
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of Run

% --- Executes on button press in togglebutton6.
function togglebutton6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton6


% % --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on button press in save_figure2.
function save_figure2_Callback(hObject, eventdata, handles)
% hObject    handle to save_figure2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in save1.
function save1_Callback(hObject, eventdata, handles)
% hObject    handle to save1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

      
  


% --- Executes during object creation, after setting all properties.
function traceplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to traceplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
axes(hObject)
imshow('panda.jpg')
guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate traceplot


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eventwindow=handles.eventwindow;
path_save=handles.mapfolder;
foname=path_save;


if get(hObject,'Value')
    
   if strcmp(handles.experimenttype,'Whole cell')
       direct_t1=handles.direct_t;
   plottrace2kawen_Gui(handles.cells,foname,direct_t1,eventwindow,handles)
   elseif strcmp(handles.experimenttype,'Cell attached')
       plottrace2kawenSpike_Gui(handles.cells,foname,eventwindow,handles)
   end
end
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes during object creation, after setting all properties.
function axes11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
axes(hObject)
imshow('dongwu17.jpg')
guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate axes11


% --- Executes during object creation, after setting all properties.
function axes12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes12


% --- Executes during object creation, after setting all properties.
function axes_plot1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_plot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_plot1


% --- Executes during object creation, after setting all properties.
function axes_plot2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_plot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_plot2


% --- Executes during object creation, after setting all properties.
function axes13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
Fcn to populate axes13


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_Gui.
function load_Gui_Callback(hObject, eventdata, handles)
% hObject    handle to load_Gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helloWorld = gcf;
hgload(handles.savedguidata);
close(helloWorld);
% Hint: get(hObject,'Value') returns toggle state of load_Gui


% --- Executes on button press in save_GUI.
function save_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to save_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=sprintf('mydata_%s',handles.mapname)
filename=fullfile(handles.save_path,filename)
 hgsave(filename)
% Hint: get(hObject,'Value') returns toggle state of save_GUI



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.savedguidata=get(hObject,'string');
 guidata(hObject, handles);
 
% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.


%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function panda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panda_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns
% called
% set(gcf,'CurrentAxes',hObject);
%     
% imshow('panda_new.jpg')
% guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate panda_new
set(gcf,'CurrentAxes',hObject);
imshow('check-in-minion.jpg')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function error_run_CreateFcn(hObject, eventdata, handles)
% hObject    handle to error_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in experimenttype.
function experimenttype_Callback(hObject, eventdata, handles)
% hObject    handle to experimenttype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 contents = cellstr(get(hObject,'String'));
 
 handles.experimenttype=contents{get(hObject,'Value')};
 handles.experimenttypeValue=handles.experimenttype;
 guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns experimenttype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from experimenttype


% --- Executes during object creation, after setting all properties.
function experimenttype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to experimenttype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function panda_new_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panda_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate panda_new

% --- Executes during object creation, after setting all properties.
function showfilepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showfilepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 guidata(hObject, handles);

% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename1,filepath]=uigetfile({'*.*', 'All file'},'Select the map to be analyzed','/Users/lindameng/Documents/RoarData/Myelin/wholecell')
 path_file=fullfile(filepath,filename1);
 handles.filepath=path_file;
 set(hObject,'string',handles.filepath)
 guidata(hObject, handles);

   


% --- Executes on button press in save_path.
function save_path_Callback(hObject, eventdata, handles)
% hObject    handle to save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.save_path=uigetdir;

if ~isfield(handles,'eventwindowValue')||isempty(handles.eventwindowValue)
        handles.eventwindow=50;
end
foldername=sprintf('%s%i','eventwindow',handles.eventwindow)
save_folder=fullfile(handles.save_path,foldername)

handles.save_folder=save_folder;
set(hObject,'string',handles.save_path)

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function direct_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to direct_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minage_Callback(hObject, eventdata, handles)
% hObject    handle to minage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.minage=str2num(get(hObject,'string'));
handles.minageValue=handles.minage;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of minage as text
%        str2double(get(hObject,'String')) returns contents of minage as a double


% --- Executes during object creation, after setting all properties.
function minage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxage_Callback(hObject, eventdata, handles)
% hObject    handle to maxage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxage=str2num(get(hObject,'string'));
handles.maxageValue=handles.maxage;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of maxage as text
%        str2double(get(hObject,'String')) returns contents of maxage as a double


% --- Executes during object creation, after setting all properties.
function maxage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in averagedrugtype.
function averagedrugtype_Callback(hObject, eventdata, handles)
% hObject    handle to averagedrugtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 contents = cellstr(get(hObject,'String'));
 
switch contents{get(hObject,'Value')}
    case {'High Mg'}

        handles.avgTp=1;
    case {'High Mg + PTX'}
        handles.avgTp=2;
    case{'High Mg + PTX +APV'}
        handles.avgTp=3;
    case{'High Mg+TTX'}
        handles.avgTp=4;
    otherwise
            handles.avgTp=5;
end

guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns averagedrugtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from averagedrugtype


% --- Executes during object creation, after setting all properties.
function averagedrugtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to averagedrugtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgdarkexp_Callback(hObject, eventdata, handles)
% hObject    handle to avgdarkexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.avgdarkexp=str2num(get(hObject,'string'));
handles.avgdarkexpValue=handles.avgdarkexp;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of avgdarkexp as text
%        str2double(get(hObject,'String')) returns contents of avgdarkexp as a double


% --- Executes during object creation, after setting all properties.
function avgdarkexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgdarkexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgHoldingPotential_Callback(hObject, eventdata, handles)
% hObject    handle to avgHoldingPotential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.avgHoldingPotential=str2num(get(hObject,'string'));
handles.avgHoldingPotentialValue=handles.avgHoldingPotential;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of avgHoldingPotential as text
%        str2double(get(hObject,'String')) returns contents of avgHoldingPotential as a double


% --- Executes during object creation, after setting all properties.
function avgHoldingPotential_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgHoldingPotential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ispushed=get(hObject,'Value');
if ispushed
     if ~isfield(handles,'eventwindowValue')||isempty(handles.eventwindowValue)
        handles.eventwindow=50;
        handles.eventwindowValue=50;
        set(handles.error_run,'string','Make sure the event window is 50ms. If not, please enter direct window on the right side panel')
        
     elseif ~isfield(handles,'direct_tValue')||isempty(handles.direct_tValue)
        handles.direct_t=8;
        handles.direct_tValue=handles.direct_t;
        set(handles.error_run,'string','Make sure the direct window is 8ms. If not, please enter direct window on the right side panel')
    else
        set(handles.error_run,'string','')
        
        
        
        set(hObject, 'String', 'Calculating...')
        drawnow
        %    guidata(hObject, handles);
        flag=0;
        
        handles= mapave_gui2(handles.avgfolder,handles)
        if flag
            display('finished')
            set(hObject, 'String', 'Ready.');
        end
        pause(2);
        drawnow update
        set(hObject, 'String', 'RUN AVERAGE MAP');
        set(hObject,'Value',0);
     end
    
end
guidata(hObject, handles);
% --- Executes on selection change in avgdata.
function avgdata_Callback(hObject, eventdata, handles)
% hObject    handle to avgdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
 
 handles.avgdata=contents{get(hObject,'Value')};
 handles.avgdataValue=handles.avgdata;
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns avgdata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from avgdata


% --- Executes during object creation, after setting all properties.
function avgdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in avgfolder.
function avgfolder_Callback(hObject, eventdata, handles)
% hObject    handle to avgfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
avgfolder=uigetdir;
handles.avgfolder=avgfolder;
set(hObject,'string',handles.avgfolder)
guidata(hObject, handles);


% --- Executes on button press in plotavg.
function plotavg_Callback(hObject, eventdata, handles)
% hObject    handle to plotavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ispushed=get(hObject,'Value');
if ispushed
    
   set(hObject, 'String', 'Calculating...')
   drawnow
%    guidata(hObject, handles);
    
handles= MeanValueCalculationSiRatio_TMC1( handles );

% if strcmp(handles.avgdata,'AveragePeak')
%     
%     Avgevent_Gui(handles.avgxaxis,handles.avgyaxis,handles.avgDensity,handles.avgPk,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
% elseif strcmp(handles.avgdata,'AverageCharge')
%     Avgevent_Gui(handles.avgxaxis,handles.avgyaxis,handles.avgDensity,handles.avgChg,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
% elseif strcmp(handles.avgdata,'AverageLatency')
%     Avgevent_Gui(handles.avgxaxis,handles.avgyaxis,handles.avgDensity,handles.Lat,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
% elseif strcmp(handles.avgdata,'Density')
%     Avgevent_Gui(handles.avgxaxis,handles.avgyaxis,handles.avgDensity,handles.avgDensity,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
% elseif strcmp(handles.avgdata,'boxplotPercentagePeak')
%     boxplot_GUI(handles.Ppeaklayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
% elseif strcmp(handles.avgdata,'boxplotPercentageCharge')
%     boxplot_GUI(handles.Pchglayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
% elseif strcmp(handles.avgdata,'boxplotSkewPeak')
%     boxplot_GUI(handles.SkPeaklayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
% elseif strcmp(handles.avgdata,'boxplotSkewCharge')
%     boxplot_GUI(handles.Skchglayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
% elseif strcmp(handles.avgdata,'boxplotDistributionCharge')
%     boxplot_GUI(handles.DistWidthCharge,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
% elseif strcmp(handles.avgdata,'boxplotDistributionPeak')
%     boxplot_GUI(handles.DistWidthPeak,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
% elseif strcmp(handles.avgdata,'AgeVsL4Bnd')    
%     ScatterPlot_AgVsBnd(handles.AvgAge,handles.Bnd,handles.avgfolder,handles.avgdata,handles.avgmap,handles)
% end

pause(2);
      drawnow update
      set(hObject, 'String', 'Map Average');
      set(hObject,'Value',0);
end

guidata(hObject, handles);



function cutoffdensity_Callback(hObject, eventdata, handles)
% hObject    handle to cutoffdensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoffdensity as text
%        str2double(get(hObject,'String')) returns contents of cutoffdensity as a double
handles.cutoffdensity=str2num(get(hObject,'string'));
handles.cutoffdensityValue=handles.cutoffdensity;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cutoffdensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoffdensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function avgmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate avgmap


% --- Executes during object creation, after setting all properties.
function axes17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 imshow('panda.jpg')
guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate axes17


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function RectifyHoldingPotential_Callback(hObject, eventdata, handles)
% hObject    handle to RectifyHoldingPotential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RectifyHoldingPotential1=str2num(get(hObject,'string'));
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of RectifyHoldingPotential as text
%        str2double(get(hObject,'String')) returns contents of RectifyHoldingPotential as a double


% --- Executes during object creation, after setting all properties.
function RectifyHoldingPotential_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RectifyHoldingPotential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function additiondrugs_Callback(hObject, eventdata, handles)
% hObject    handle to additiondrugs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.additiondrugs=get(hObject,'String');
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of additiondrugs as text
%        str2double(get(hObject,'String')) returns contents of additiondrugs as a double


% --- Executes during object creation, after setting all properties.
function additiondrugs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to additiondrugs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function laserstim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laserstim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
% function minions_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to minions (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% set(gcf,'CurrentAxes',hObject);
% imshow('check-in-minion.jpg')
% guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate minions


% --- Executes on button press in allmaps.
function allmaps_Callback(hObject, eventdata, handles)
% hObject    handle to allmaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  get(hObject,'Value')
    eventwindow=handles.eventwindow;
    darkexp=handles.darkexp;
    mapname=handles.mapname;
    
    
%     path_save=handles.save_path;
    foname=handles.mapfolder;
    
    
    flipimg=handles.flipimg;
    header=handles.cells.header;
   
    LaserStim=handles.laserstim;
    fdata=handles.cells.data;
    age=handles.age;
    plot_fig=1;
    Tp=handles.Tp;
    flipimg2=handles.flipimg2;
    
    if LaserStim
        
        [Stim_coor,header]=Matrixrotate_singlemap(header,flipimg,flipimg2);
        X=Stim_coor(1,:);
        Y=Stim_coor(2,:);
        
        
        if strcmp(handles.experimenttype,'Cell attached')
            
            SpikeLatencyImage_Gui(mapname,foname,Tp,header,'SpikeLatency+DIC image',age,handles.cells,handles.axes_plot2,eventwindow);
             SpikeLatencyImage_Gui(mapname,foname,Tp,header,'SpikeLatency inside eventswindow',age,handles.cells,handles.axes_plot2,eventwindow);
            Spikeacummulation_Gui(mapname,foname,Tp,header,'AcummulateDistributionSpikelatencyEventwindow',age,handles.cells,handles.axes_plot1,eventwindow)
            
            Spikeacummulation_Gui(mapname,foname,Tp,header,'AcummulateDstributionSpikelatency',age,handles.cells,handles.axes_plot1,eventwindow)
            
            plottrace2kawenSpike_Gui(handles.cells,foname,eventwindow,handles);
        elseif strcmp(handles.experimenttype,'Whole cell')
             BD=handles.cells.Boundary;
            direct_t1=handles.direct_t;
            direct_t2=handles.direct_t;
            events=handles.cells.events;
           
            handles=plotsinglemap2_Gui(X,Y,events.peakAmp,mapname,foname,Tp,header,plot_fig,'peak+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles);
            
            [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.area,mapname,foname,Tp,header,plot_fig,'charge',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles,BD);
            
            [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.startSamp,mapname,foname,Tp,header,plot_fig,'Latency',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles,BD);
            
            [avg,dist_x_chg,dist_y_chg,handles]=singlemapkawen_GUI(X,Y,events.peakAmp,mapname,foname,Tp,header,plot_fig,'peak',events,direct_t1,direct_t2,eventwindow,handles.axes_plot1,hObject,handles,BD);
            
            handles=plotsinglemap2_Gui(X,Y,events.area,mapname,foname,Tp,header,plot_fig,'charge+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles);
            
            handles=plotsinglemap2_Gui(X,Y,events.startSamp,mapname,foname,Tp,header,plot_fig,'Latency+Direct',events,direct_t1,direct_t2,eventwindow,handles.axes_plot2,hObject,handles);
            
            
            maptraces_Gui(fdata,header,mapname,foname,Tp,eventwindow,handles.axes_plot2,events)
             plottrace2kawen_Gui(handles.cells,foname,direct_t1,eventwindow,handles)
        end
    end
end
      
% Hint: get(hObject,'Value') returns toggle state of allmaps


% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ispushed=get(hObject,'Value');
if ispushed
    set(gcf,'CurrentAxes',handles.axes_plot2);
    hold on
    [x,y] = ginput;
    plot(x,y,'Oy')
    header=handles.cells.header;
    fdata=handles.cells.data;
    flag=handles.cells.events.flag;
    k=length(x);Indx=[];
    for i=1:k
        dist=(x(i)-header.StimCoordinates(2,:)).^2+(y(i)-header.StimCoordinates(1,:)).^2;
        Id=find(dist<=55);
        Idf=[];
        for ii=1:length(Id)
        if (flag{Id(ii)}(1)>0);
        Idf=[Idf,Id(ii)];
        end
        end
        if length(Idf)
        [mxy,I]=min(dist(Idf));
        Indx=[Indx;Idf(I(1))];
        else
        [mxy,I]=min(dist);
        Indx=[Indx;(I(1))];
        end
        
    end
    currfig=find_figure('SelectTraces')
    for i=1:k
        subplot(k,1,i)
        plot(fdata(1:10:end,Indx(i)))
    end
    hold off
    
end
        


% --- Executes during object creation, after setting all properties.
function pandaPic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pandaPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 imshow('panda.jpg')
guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate pandaPic



function mindist_Callback(hObject, eventdata, handles)
% hObject    handle to mindist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mindist=str2num(get(hObject,'string'));
handles.mindistValue=handles.mindist;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of mindist as text
%        str2double(get(hObject,'String')) returns contents of mindist as a double


% --- Executes during object creation, after setting all properties.
function mindist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mindist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function maxdist_Callback(hObject, eventdata, handles)
% hObject    handle to maxdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxdist=str2num(get(hObject,'string'));
handles.maxdistValue=handles.maxdist;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of maxdist as text
%        str2double(get(hObject,'String')) returns contents of maxdist as a double


% --- Executes during object creation, after setting all properties.
function maxdist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in layeravg.
function layeravg_Callback(hObject, eventdata, handles)
% hObject    handle to layeravg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
 handles.layeravg=contents{get(hObject,'Value')};
switch contents{get(hObject,'Value')}
    case {'Layer1'}

        handles.indb=2;
    case {'Layer2/3'}
        handles.indb=2;
    case{'Layer4'}
        handles.indb=2;
    case{'Layer5/6'}
        handles.indb=3;
    case{'Subplate'}
        handles.indb=4;
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns layeravg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from layeravg


% --- Executes during object creation, after setting all properties.
function layeravg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layeravg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minmaps_Callback(hObject, eventdata, handles)
% hObject    handle to minmaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.minmaps=str2num(get(hObject,'string'));
handles.minmapsValue=handles.minmaps;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of minmaps as text
%        str2double(get(hObject,'String')) returns contents of minmaps as a double


% --- Executes during object creation, after setting all properties.
function minmaps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minmaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function panda_pic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panda_pic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 imshow('panda.jpg')
guidata(hObject, handles);
% Hint: place code in OpeningFcn to populate panda_pic


% --- Executes on button press in plotaverage.
function plotaverage_Callback(hObject, eventdata, handles)
% hObject    handle to plotaverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ispushed=get(hObject,'Value');
if ispushed
    
   set(hObject, 'String', 'Calculating...')
   drawnow
%    guidata(hObject, handles);
    
% handles= MeanValueCalculationSi( handles );
datafolder=handles.avgfolder;
if isfield(handles,'folderMeanpath')
    folderMeanpath=handles.folderMeanpath;
else
    folderMeanname=sprintf('meanValueeachcell_eventwindow%i_direct%i',handles.eventwindow,handles.direct_t);
    folderMeanpath=fullfile(datafolder,folderMeanname);
end
filesum=sprintf('SilentsynapseSummaryAge%ito%i.mat',handles.minage,handles.maxage)
load(fullfile(folderMeanpath,filesum));
hd=Cellsum;
handles.dX=hd.dX;
handles.avgCellNum=hd.avgCellNum;
if strcmp(handles.avgdata,'AveragePeak')
    
    Avgevent_Gui(hd.avgxaxis,hd.avgyaxis,hd.avgDensity,hd.avgPk,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
elseif strcmp(handles.avgdata,'AverageCharge')
    Avgevent_Gui(hd.avgxaxis,hd.avgyaxis,hd.avgDensity,hd.avgChg,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
elseif strcmp(handles.avgdata,'AverageLatency')
    Avgevent_Gui(hd.avgxaxis,hd.avgyaxis,hd.avgDensity,hd.Lat,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
elseif strcmp(handles.avgdata,'Density')
    Avgevent_Gui(hd.avgxaxis,hd.avgyaxis,hd.avgDensity,hd.avgDensity,handles.cutoffdensity,handles.avgdata,handles.avgfolder,handles,handles.avgmap,handles.eventwindow)
elseif strcmp(handles.avgdata,'boxplotPercentagePeak')
    boxplot_GUI(hd.Ppeaklayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'boxplotPercentageCharge')
    boxplot_GUI(hd.Pchglayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'boxplotSkewPeak')
    boxplot_GUI(hd.SkPeaklayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'boxplotDirArea')
    boxplot_GUI(hd.DirArea,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'boxplotSkewCharge')
    boxplot_GUI(hd.Skchglayer,handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'ScatterplotDistributionCharge')
    
    scatterdist(hd.DistWidthCharge(:,1),hd.DistWidthCharge(:,3),handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'ScatterplotDistributionPeak')
   scatterdist(hd.DistWidthPeak(:,1),hd.DistWidthPeak(:,3),handles.avgfolder,handles.avgdata,handles,handles.avgmap)
elseif strcmp(handles.avgdata,'AgeVsL4Bnd')    
    ScatterPlot_AgVsBnd(hd.AvgAge,hd.Bnd,handles.avgfolder,handles.avgdata,handles.avgmap,handles)
end

pause(2);
      drawnow update
      set(hObject, 'String', 'Map Average');
      set(hObject,'Value',0);
end
guidata(hObject, handles);


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
 load(handles.filepath);
 handles.cells=cells;
 handles.eventwindow=cells.eventwindow;
handles.darkexp=cells.darkexp;
handles.mapname=cells.mapname;

% load data
% handles.save_path=cells.cell_folder;
handles.mapfolder=cells.cell_folder;
% experiment=handles.experiment
% subdir=handles.cell_folder;
handles.flipimg=cells.flipimg;
handles.age=cells.age;

handles.Tp=cells.Tp;
handles.laserstim=cells.LaserStim;
handles.direct_t=cells.direct_t;
handles.experimenttype=cells.experimenttype
% if cmpstr(handles.experimenttype,'cellattach')
%     handles.threshold=-20;
% else
%     handles.threshold
% end
if ~isfield(handles,'thresholdValue')|| isempty(handles.thresholdValue)
    if strcmp(handles.experimenttype,'Cell attached')
        handles.threshold=-20;
    elseif strcmp(handles.experimenttype,'Whole cell')
    handles.threshold=10;
    end
else
    handles.threshold=handles.thresholdValue;

end
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox5



function BD_Callback(hObject, eventdata, handles)
% hObject    handle to BD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BD=str2num(get(hObject,'string'));
handles.BDValue=handles.BD;
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of BD as text
%        str2double(get(hObject,'String')) returns contents of BD as a double


% --- Executes during object creation, after setting all properties.
function BD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flipimg2=str2num(get(hObject,'string'));
handles.flipimgValue2=handles.flipimg2;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
