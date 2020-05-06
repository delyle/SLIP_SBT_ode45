function varargout = QuadWalkUI(varargin)
%QUADWALKUI M-file for QuadWalkUI.fig
%      QUADWALKUI, by itself, creates a new QUADWALKUI or raises the existing
%      singleton*.
%
%      H = QUADWALKUI returns the handle to a new QUADWALKUI or the handle to
%      the existing singleton*.
%
%      QUADWALKUI('Property','Value',...) creates a new QUADWALKUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to QuadWalkUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      QUADWALKUI('CALLBACK') and QUADWALKUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in QUADWALKUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QuadWalkUI

% Last Modified by GUIDE v2.5 18-Jun-2015 17:45:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuadWalkUI_OpeningFcn, ...
                   'gui_OutputFcn',  @QuadWalkUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before QuadWalkUI is made visible.
function QuadWalkUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

plot_stuff(hObject, eventdata, handles)
handles.HAR21fz_out.String = num2str(get(handles.HAR21fz,'Value'));
handles.HAR21hz_out.String = num2str(get(handles.HAR21hz,'Value'));
handles.HAR31fz_out.String = num2str(get(handles.HAR31fz,'Value'));
handles.HAR31hz_out.String = num2str(get(handles.HAR31hz,'Value'));
handles.HAR42fy_out.String = num2str(get(handles.HAR42fy,'Value'));
handles.HAR42hy_out.String = num2str(get(handles.HAR42hy,'Value'));
handles.DFf_out.String = num2str(get(handles.DFf,'Value'));
handles.DFh_out.String = num2str(get(handles.DFh,'Value'));
handles.DFh_out.String = num2str(get(handles.DFh,'Value'));
handles.Td_off.String = num2str(get(handles.Ttdh,'Value'));

% Choose default command line output for QuadWalkUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QuadWalkUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = QuadWalkUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function HAR21fz_Callback(hObject, eventdata, handles)
% hObject    handle to HAR21fz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.HAR21fz_out.String = num2str(get(handles.HAR21fz,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function HAR21fz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HAR21fz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function HAR21hz_Callback(hObject, eventdata, handles)
% hObject    handle to HAR21hz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.HAR21hz_out.String = num2str(get(handles.HAR21hz,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function HAR21hz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HAR21hz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function HAR31fz_Callback(hObject, eventdata, handles)
% hObject    handle to HAR31fz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.HAR31fz_out.String = num2str(get(handles.HAR31fz,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function HAR31fz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HAR31fz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function HAR42fy_Callback(hObject, eventdata, handles)
% hObject    handle to HAR42fy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.HAR42fy_out.String = num2str(get(handles.HAR42fy,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function HAR42fy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HAR42fy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function DFf_Callback(hObject, eventdata, handles)
% hObject    handle to DFf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.DFf_out.String = num2str(get(handles.DFf,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function DFf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DFf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function HAR31hz_Callback(hObject, eventdata, handles)
% hObject    handle to HAR31hz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.HAR31hz_out.String = num2str(get(handles.HAR31hz,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function HAR31hz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HAR31hz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function HAR42hy_Callback(hObject, eventdata, handles)
% hObject    handle to HAR42hy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.HAR42hy_out.String = num2str(get(handles.HAR42hy,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function HAR42hy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HAR42hy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function DFh_Callback(hObject, eventdata, handles)
% hObject    handle to DFh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.DFh_out.String = num2str(get(handles.DFh,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function DFh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DFh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function HAR21fz_out_Callback(hObject, eventdata, handles)
% hObject    handle to HAR21fz_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.HAR21fz.Value = str2double(handles.HAR21fz_out.String);
plot_stuff(hObject, eventdata, handles)


% Hints: get(hObject,'String') returns contents of HAR21fz_out as text
%        str2double(get(hObject,'String')) returns contents of HAR21fz_out as a double



function FHRz_Callback(hObject, eventdata, handles)
% hObject    handle to FHRz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of FHRz as text
%        str2double(get(hObject,'String')) returns contents of FHRz as a double


% --- Executes during object creation, after setting all properties.
function FHRz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FHRz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FHRy_Callback(hObject, eventdata, handles)
% hObject    handle to FHRy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of FHRy as text
%        str2double(get(hObject,'String')) returns contents of FHRy as a double


% --- Executes during object creation, after setting all properties.
function FHRy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FHRy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Ttdh_Callback(hObject, eventdata, handles)
% hObject    handle to Ttdh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)
handles.Td_off.String = num2str(get(handles.Ttdh,'Value'));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Ttdh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ttdh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Td_off_Callback(hObject, eventdata, handles)
% hObject    handle to Td_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ttdh.Value = str2double(handles.Td_off.String);
plot_stuff(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Td_off as text
%        str2double(get(hObject,'String')) returns contents of Td_off as a double


% --- Executes during object creation, after setting all properties.
function Td_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Td_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- This function is only here so that I don't get an error
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Td_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when called, makes computations and plots the data --- %
function plot_stuff(hObject, eventdata, handles)
% Vertical Forces
FHRz = str2double(get(handles.FHRz,'String')); % The ratio of front to hind peak vertical forces
sh21Rzf = get(handles.HAR21fz,'Value'); % Forelimbs: The relative amplitude of the second to first harmonic
sh31Rzf = get(handles.HAR31fz,'Value'); % Forelimbs: The relative amplitude of the third to first harmonic
sh21Rzh = get(handles.HAR21hz,'Value'); % Hindlimbs: The relative amplitude of the second to first harmonic
sh31Rzh = get(handles.HAR31hz,'Value'); % Hindlimbs: The relative amplitude of the third to first harmonic

% Anterior forces
FHRy = str2double(get(handles.FHRy,'String')); % The ratio of front to hind peak anterior forces
sh42Ryf = get(handles.HAR42fy,'Value'); % Forelimbs: The relative amplitude of the fourth to second harmonic
sh42Ryh = get(handles.HAR42hy,'Value'); % Hindlimbs: The relative amplitude of the fourth to second harmonic

% Temporal variables
Tsf = get(handles.DFf,'Value'); % The duty factor of the forelimbs
Tsh = get(handles.DFh,'Value'); % The duty factor of the hindlimbs
Ttdh = get(handles.Ttdh,'Value'); % The touchdown time of the hindlimbs, relative to previous forelimb touchdown
                                  % normalized by stride period

% Average forward velocity
vy_bar = str2double(get(handles.Vy_bar,'String'));

% --- Establish forces during stance --- %:

t = (0:0.001:1)';
Fz_f = sin(pi*t)+sh21Rzf*sin(2*pi*t)+sh31Rzf*sin(3*pi*t);
Fz_h = (sin(pi*t)+sh21Rzh*sin(2*pi*t)+sh31Rzh*sin(3*pi*t))/FHRz;

Fy_f = 0.15*(-sin(2*pi*t)+sh42Ryf*sin(4*pi*t));
Fy_h = 0.15*(-sin(2*pi*t)+sh42Ryh*sin(4*pi*t))/FHRy;


% --- Translate forces into normalized stride time --- %

% First need to convert forces to the proper format:
Fz_f = [t, Fz_f]; Fz_h = [t, Fz_h];
Fy_f = [t, Fy_f]; Fy_h = [t, Fy_h];

% convert to stride time. Note that RF touchdown is t = 0;

t0R = 0; tfR = Tsf;
Fz_Rf = stance2gait(Fz_f,t0R,tfR,t); Fy_Rf = stance2gait(Fy_f,t0R,tfR,t);
t0L = 0.5; tfL = 0.5+Tsf;
Fz_Lf = stance2gait(Fz_f,t0L,tfL,t); Fy_Lf = stance2gait(Fy_f,t0L,tfL,t);
t0L = Ttdh; tfL = Ttdh+Tsh;
Fz_Lh = stance2gait(Fz_h,t0L,tfL,t); Fy_Lh = stance2gait(Fy_h,t0L,tfL,t);
t0R = 0.5 + Ttdh; tfR = 0.5+Ttdh+Tsh;
Fz_Rh = stance2gait(Fz_h,t0R,tfR,t); Fy_Rh = stance2gait(Fy_h,t0R,tfR,t);


% --- Determine total forces --- %

Fz_tot = Fz_Rf+Fz_Rh+Fz_Lf+Fz_Lh;  Fy_tot = Fy_Rf+Fy_Rh+Fy_Lf+Fy_Lh;

% scale the vertical forces so that the animal is supported
scale = mean(Fz_tot);
Fz_tot = Fz_tot/scale;
Fz_Lf = Fz_Lf/scale; Fz_Rf = Fz_Rf/scale;
Fz_Lh = Fz_Lh/scale; Fz_Rh = Fz_Rh/scale;


% --- Determine the velocities --- %


% Now calculate the CoM velocity

vy = cumtrapz(t,Fy_tot)';
vy0 = -trapz(t,vy);
vy = vy + vy0;

vz = cumtrapz(t,Fz_tot)';
vz0 = 0 + 1/2 - trapz(t,vz);
vz = vz0 + vz -t;

% Calculate CoM vertical position
    % normalize z so that zo = 0
z = cumtrapz(t,vz);

% Finally, calculate CoM energy
EP = -z;
EK = 1/2*sqrt(vy.^2+vz.^2);

% Calculate abs(phi)

V = [vy'+vy_bar;vz'];
%F = [Fy_tot;Fz_tot]; 
%dCoT = abs(dot(F,V));
dCoT1 = subplus(dot([Fy_Lf;Fz_Lf],V));
dCoT2 = subplus(dot([Fy_Rf;Fz_Rf],V));
dCoT3 = subplus(dot([Fy_Lh;Fz_Lh],V));
dCoT4 = subplus(dot([Fy_Rh;Fz_Rh],V));
dCoT = dCoT1 + dCoT2 + dCoT3 + dCoT4;
CoT = trapz(t,dCoT)/vy_bar; 
%CoT = (trapz(t,dCoT1)+trapz(t,dCoT2)+trapz(t,dCoT3)+trapz(t,dCoT4))/vy_bar;

handles.CoT.String = ['Cost of Transport: ', num2str(CoT,5)];

% ----- Make a figure ----- %

% Plot vert forces
axes(handles.axes1)
    cla
    hold on
    plot(t,Fz_Lf)%plot(t,Fz_Lf)
    plot(t,Fz_Rf)
    plot(t,Fz_Lh, '--')
    plot(t,Fz_Rh, '--')
    ylabel('Vertical F^*')
    box on
    legend('LF','RF','LH','RH','location', 'southwest')
    
% Plot total forces
 
axes(handles.axes2)
    cla
    hold on
    plot(t,[Fz_tot-1;Fy_tot])
    plot([t(1) t(end)],[0 0],'k--')
    ylabel('Total F^*')
    legend('F_z - 1','F_y','location','southwest')
    box on

% Plot anterior forces from each leg
axes(handles.axes3)
    cla
    hold on
    plot(t,Fy_Lf)
    plot(t,Fy_Rf)
    plot(t,Fy_Lh, '--')
    plot(t,Fy_Rh, '--')
    ylabel('Anterior F^*')
    box on
    
% Plot COM Energy
axes(handles.axes4)
    cla
    hold on
    plot(t,[EP,EK+EP])
    ylabel('E^*')
    xlabel('t^*')
    box on
    legend('E_p','E_{tot}','location','southwest')
    
% Plot hodograph

axes(handles.axes5)
    cla
    hold on
    plot(vy,vz);
    h(:,1) = plot(vy(1),vz(1), 'go', 'markerfacecolor', 'g'); % start marker
    u = vy(20)-vy(1);
    v = vz(20)-vz(1);
    h(:,2) = quiver(vy(1),vz(1),u,v,'color','k','maxheadsize',1);
    box on
    xlabel('v_y^* - v^*_{yAv}')
    ylabel('v^*_z')
    legend(h,{'start','direction'},'location','southeast')

% Plot COM velocities

axes(handles.axes6)
    cla
    hold on
    plot(t,vz/max(vz))
    plot(t,vy/max(vy))
    plot([t(1) t(end)],[0 0],'k--')
    ylabel('(v^* - v_{Av}^*) / v^*_{max}')
    xlabel('t^*')
    legend('V_z','V_y','location','southeast')
    box on

% Plot dCoT

axes(handles.axes7)
    cla
    hold on
    plot(t(t >= 0 & t <= 0.5),dCoT(t >= 0 & t <= 0.5))
    xlabel('t^*')
    ylabel('|\bf{F.V}|')
    box on


function HAR31fz_out_Callback(hObject, eventdata, handles)
% hObject    handle to HAR31fz_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HAR31fz.Value = str2double(handles.HAR31fz_out.String);
plot_stuff(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of HAR31fz_out as text
%        str2double(get(hObject,'String')) returns contents of HAR31fz_out as a double



function HAR21hz_out_Callback(hObject, eventdata, handles)
% hObject    handle to HAR21hz_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HAR21hz.Value = str2double(handles.HAR21hz_out.String);
plot_stuff(hObject, eventdata, handles)


% Hints: get(hObject,'String') returns contents of HAR21hz_out as text
%        str2double(get(hObject,'String')) returns contents of HAR21hz_out as a double



function HAR42fy_out_Callback(hObject, eventdata, handles)
% hObject    handle to HAR42fy_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HAR42fy.Value = str2double(handles.HAR42fy_out.String);
plot_stuff(hObject, eventdata, handles)


% Hints: get(hObject,'String') returns contents of HAR42fy_out as text
%        str2double(get(hObject,'String')) returns contents of HAR42fy_out as a double



function HAR31hz_out_Callback(hObject, eventdata, handles)
% hObject    handle to HAR31hz_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HAR31hz.Value = str2double(handles.HAR31hz_out.String);
plot_stuff(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of HAR31hz_out as text
%        str2double(get(hObject,'String')) returns contents of HAR31hz_out as a double



function HAR42hy_out_Callback(hObject, eventdata, handles)
% hObject    handle to HAR42hy_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HAR42hy.Value = str2double(handles.HAR42hy_out.String);
plot_stuff(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of HAR42hy_out as text
%        str2double(get(hObject,'String')) returns contents of HAR42hy_out as a double



function DFf_out_Callback(hObject, eventdata, handles)
% hObject    handle to DFf_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DFf.Value = str2double(handles.DFf_out.String);
plot_stuff(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of DFf_out as text
%        str2double(get(hObject,'String')) returns contents of DFf_out as a double



function DFh_out_Callback(hObject, eventdata, handles)
% hObject    handle to DFh_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DFh.Value = str2double(handles.DFh_out.String);
plot_stuff(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of DFh_out as text
%        str2double(get(hObject,'String')) returns contents of DFh_out as a double


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Vertical Forces
Forces.Vertical.FHRz = str2double(get(handles.FHRz,'String')); % The ratio of front to hind peak vertical forces
Forces.Vertical.sh21Rzf = get(handles.HAR21fz,'Value'); % Forelimbs: The relative amplitude of the second to first harmonic
Forces.Vertical.sh31Rzf = get(handles.HAR31fz,'Value'); % Forelimbs: The relative amplitude of the third to first harmonic
Forces.Vertical.sh21Rzh = get(handles.HAR21hz,'Value'); % Hindlimbs: The relative amplitude of the second to first harmonic
Forces.Vertical.sh31Rzh = get(handles.HAR31hz,'Value'); % Hindlimbs: The relative amplitude of the third to first harmonic

% Anterior forces
Forces.Horizontal.FHRy = str2double(get(handles.FHRy,'String')); % The ratio of front to hind peak anterior forces
Forces.Horizontal.sh42Ryf = get(handles.HAR42fy,'Value'); % Forelimbs: The relative amplitude of the fourth to second harmonic
Forces.Horizontal.sh42Ryh = get(handles.HAR42hy,'Value'); % Hindlimbs: The relative amplitude of the fourth to second harmonic

Times.Tsf = get(handles.DFf,'Value'); % The duty factor of the forelimbs
Times.Tsh = get(handles.DFh,'Value'); % The duty factor of the hindlimbs
Times.Ttdh = get(handles.Ttdh,'Value'); % The touchdown time of the hindlimbs, relative to previous forelimb touchdown
                                  % normalized by stride period
save_name = handles.edit4.String;
make_figure_synthetic_vars(Forces,Times,save_name)



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Reset_button.
function Reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1)
load('save_handles.mat')
save_axes = save_handles;
handles = save_handles;
plot_stuff(hObject, eventdata, handles)

% --- Executes on button press in SCS_button.
function SCS_button_Callback(hObject, eventdata, handles)
% hObject    handle to SCS_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_handles = handles;
save('save_handles.mat','save_handles')
close(handles.figure1)
load('save_handles.mat')


% --- Executes during object creation.
function CoT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function CoT_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CoT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Vy_bar_Callback(hObject, eventdata, handles)
% hObject    handle to Vy_bar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_stuff(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Vy_bar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vy_bar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_test_data.
function Save_test_data_Callback(hObject, eventdata, handles)
% hObject    handle to Save_test_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Vertical Forces
FHRz = str2double(get(handles.FHRz,'String')); % The ratio of front to hind peak vertical forces
sh21Rzf = get(handles.HAR21fz,'Value'); % Forelimbs: The relative amplitude of the second to first harmonic
sh31Rzf = get(handles.HAR31fz,'Value'); % Forelimbs: The relative amplitude of the third to first harmonic
sh21Rzh = get(handles.HAR21hz,'Value'); % Hindlimbs: The relative amplitude of the second to first harmonic
sh31Rzh = get(handles.HAR31hz,'Value'); % Hindlimbs: The relative amplitude of the third to first harmonic

% Anterior forces
FHRy = str2double(get(handles.FHRy,'String')); % The ratio of front to hind peak anterior forces
sh42Ryf = get(handles.HAR42fy,'Value'); % Forelimbs: The relative amplitude of the fourth to second harmonic
sh42Ryh = get(handles.HAR42hy,'Value'); % Hindlimbs: The relative amplitude of the fourth to second harmonic

% Average forward velocity
vy_bar = str2double(get(handles.Vy_bar,'String'));

% --- Establish forces during stance --- %:

t = (0:0.01:1)';
Fz_f = [t,sin(pi*t)+sh21Rzf*sin(2*pi*t)+sh31Rzf*sin(3*pi*t)];
Fz_h = [t,(sin(pi*t)+sh21Rzh*sin(2*pi*t)+sh31Rzh*sin(3*pi*t))/FHRz];

Fy_f = [t,0.15*(-sin(2*pi*t)+sh42Ryf*sin(4*pi*t))];
Fy_h = [t,0.15*(-sin(2*pi*t)+sh42Ryh*sin(4*pi*t))/FHRy];

prefix = datestr(datetime('now'),'yyyyMMdd');
filename = [prefix,'_test_data.mat'];
save(filename,'Fz_f','Fz_h','Fy_f','Fy_h','vy_bar')
disp(['Saved data to ', filename])
