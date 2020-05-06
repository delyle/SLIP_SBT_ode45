function varargout = hodograph_ui2(varargin)
% hodograph_ui2 MATLAB code for hodograph_ui2.fig
%      hodograph_ui2, by itself, creates a new hodograph_ui2 or raises the existing
%      singleton*.
%
%      H = hodograph_ui2 returns the handle to a new hodograph_ui2 or the handle to
%      the existing singleton*.
%
%      hodograph_ui2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in hodograph_ui2.M with the given input arguments.
%
%      hodograph_ui2('Property','Value',...) creates a new hodograph_ui2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hodograph_ui2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hodograph_ui2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hodograph_ui2

% Last Modified by GUIDE v2.5 08-May-2015 17:01:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hodograph_ui2_OpeningFcn, ...
                   'gui_OutputFcn',  @hodograph_ui2_OutputFcn, ...
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


% --- Executes just before hodograph_ui2 is made visible.
function hodograph_ui2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hodograph_ui2 (see VARARGIN)
f1 = str2double(get(handles.f1,'String'));
f2 = str2double(get(handles.f2,'String'));
offset = get(handles.offset,'Value');
    t = 0:0.001:1;
    vy = sin(2*pi*f1*t);
    vx = sin(2*pi*(f2*t+offset));
axes(handles.axes1)
    plot(t,[vx;vy])
    legend('Vx','Vy','Location','northoutside','orientation','horizontal')
axes(handles.axes2)
    cla
    hold on
    plot(vx,vy)
    h(:,1) = plot(vx(1),vy(1), 'ko');
    u = vx(10)-vx(1);
    v = vy(10)-vy(1);
    h(:,2) = quiver(vx(1),vy(1),u,v,'color','k','maxheadsize',1);
    box on
    axis([-1 1 -1 1])
    legend(h,{'Start','Direction'},'Location','northoutside','orientation','horizontal')
    
% Choose default command line output for hodograph_ui2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes hodograph_ui2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hodograph_ui2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f1 = str2double(get(handles.f1,'String'));
f2 = str2double(get(handles.f2,'String'));
offset = get(hObject,'Value');
    t = 0:0.001:1;
    vy = sin(2*pi*f1*t);
    vx = sin(2*pi*(f2*t+offset));
axes(handles.axes1)
    plot(t,[vx;vy])
    legend('Vx','Vy','Location','northoutside','orientation','horizontal')
axes(handles.axes2)
    cla
    hold on
    plot(vx,vy)
    h(:,1) = plot(vx(1),vy(1), 'ko');
    u = vx(10)-vx(1);
    v = vy(10)-vy(1);
    h(:,2) = quiver(vx(1),vy(1),u,v,'color','k','maxheadsize',1);
    box on
    axis([-1 1 -1 1])
    legend(h,{'Start','Direction'},'Location','northoutside','orientation','horizontal')
    




% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function f1_Callback(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Hints: get(hObject,'String') returns contents of f1 as text
%        str2double(get(hObject,'String')) returns contents of f1 as a double


% --- Executes during object creation, after setting all properties.
function f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f2_Callback(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f2 as text
%        str2double(get(hObject,'String')) returns contents of f2 as a double


% --- Executes during object creation, after setting all properties.
function f2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
