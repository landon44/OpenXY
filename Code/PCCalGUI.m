function varargout = PCCalGUI(varargin)
% PCCALGUI MATLAB code for PCCalGUI.fig
%      PCCALGUI, by itself, creates a new PCCALGUI or raises the existing
%      singleton*.
%
%      H = PCCALGUI returns the handle to a new PCCALGUI or the handle to
%      the existing singleton*.
%
%      PCCALGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCCALGUI.M with the given input arguments.
%
%      PCCALGUI('Property','Value',...) creates a new PCCALGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCCalGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCCalGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCCalGUI

% Last Modified by GUIDE v2.5 08-Jul-2015 12:01:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PCCalGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PCCalGUI_OutputFcn, ...
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

% --- Executes just before PCCalGUI is made visible.
function PCCalGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCCalGUI (see VARARGIN)
% Choose default command line output for PCCalGUI


%get Settings sent in if this was opened before and is being reopened.
if length(varargin) == 1
   Settings = varargin{1}; 

else
    %otherwise load the Settings(should be default)
    load Settings
end

%[SquareFileVals, ScanParams] = ReadScanFile(Settings.ScanFilePath);

IQ = Settings.IQ;
Nx = Settings.Nx;
Ny = Settings.Ny;
ScanType = Settings.ScanType;

switch ScanType

    case 'Square'
        IQPlot = reshape(IQ, Nx,Ny)';
        if Ny == 1 %Lines Scans
            IQPlot = repmat(IQPlot,floor(Settings.ScanLength/4),1);
        end

        axes(handles.axes1)
        imagesc(IQPlot)
        axis image %scales to natural width and height
        
    case 'Hexagonal'
        NumColsEven = floor(Nx/2);
        NumColsOdd = ceil(Nx/2);
        x = Settings.XData;
        y = Settings.YData;
        x = Hex2Array(x,NumColsOdd,NumColsEven);
        y = Hex2Array(y,NumColsOdd,NumColsEven);
        iq = Hex2Array(IQ,NumColsOdd,NumColsEven);
        
        StdDev = std(iq(:));
        Mean = mean(iq(:));
        Limits(1) = Mean - 3*StdDev;
        Limits(2) = Mean + 3*StdDev;
        
        axes(handles.axes1)
        surf(x,y,iq, 'EdgeColor','none');
        view(2);
        caxis(Limits);
        IQPlot.x = x;
        IQPlot.y = y;
        IQPlot.iq = iq;
        IQPlot.Limits = Limits;
        
end
axes(handles.axes2)
plot3(Settings.ScanParams.xstar,Settings.ScanParams.ystar,Settings.ScanParams.zstar,'bo')
colormap jet

if Settings.ImageTag
    set(handles.fromtiffbutton,'Enable','on');
else   
    set(handles.fromtiffbutton,'Enable','off');
end
set(handles.fromtiff,'Enable','off');

set(handles.SavePCCal,'Value',1);
handles.SaveAllPC = 1;

% Update handles structure
handles.Settings = Settings;
handles.ScanParams = Settings.ScanParams;
handles.VanPont = 0;
handles.calibrated = 0;
handles.tiffread = 0;
handles.IQPlot = IQPlot;
guidata(hObject, handles);
%uiwait(handles.figure1)

% --- Outputs from this function are returned to the command line.
function varargout = PCCalGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = 1;

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in savenclose.
function savenclose_Callback(hObject, eventdata, handles)
% hObject    handle to savenclose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
planefitpanel_SelectionChangeFcn(handles.planefitpanel, eventdata, handles);
handles = guidata(hObject);
YesNo = 'Yes';
if ~handles.calibrated && ~handles.tiffread
    YesNo = questdlg({'No pattern center calibration has been performed.'; 'Are you sure you want to continue?'},'No Calibration');
end
if ~strcmp(YesNo,'No')
    Settings = handles.Settings;
    SaveAllPC = handles.SaveAllPC;
    
    if ~handles.calibrated
        Settings.DoPCStrainMin = 0;
    elseif SaveAllPC
        PCCal.MeanXstar = handles.MeanXstar;
        PCCal.MeanYstar = handles.MeanYstar;
        PCCal.MeanZstar = handles.MeanZstar;
        PCCal.FitXstar = handles.FitXstar;
        PCCal.FitYstar = handles.FitYstar;
        PCCal.FitZstar = handles.FitZstar;
        PCCal.NaiveXstar = handles.NaiveXstar;
        PCCal.NaiveYstar = handles.NaiveYstar;
        PCCal.NaiveZstar = handles.NaiveZstar;
        if handles.tiffread
            PCCal.TiffXstar = handles.TiffXstar;
            PCCal.TiffYstar = handles.TiffYstar;
            PCCal.TiffZstar = handles.TiffZstar;
        end
        Settings.PCCal = PCCal;
    end
    if strcmp(YesNo,'Cancel')
        Settings.Exit = 1;
    end
    save('Settings.mat','Settings');
    delete(handles.figure1);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
try
    dlg = 'Continue';
    if handles.calibrated
        dlg = questdlg('Please select an option:', 'Continue','Cancel','Continue');
    end
    switch dlg
        case 'Continue'
            savenclose_Callback(handles.savenclose,eventdata,handles);
        case 'Cancel'
            Settings = handles.Settings;
            Settings.Exit = 1;
            save('Settings.mat','Settings');
            delete(handles.figure1);
    end
catch
    delete(hObject)
end

% --- Executes on button press in calibratebutton.
function calibratebutton_Callback(hObject, eventdata, handles)
% hObject    handle to calibratebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.VanPont

    
    Settings = handles.Settings;
    ScanParams = handles.ScanParams;

    npoints = length(Settings.CalibrationPointIndecies);

    CalibrationPointsPC = zeros(npoints,3);
    
    NumCores = Settings.DoParallel;
    try
        ppool = gcp('nocreate');
        if isempty(ppool) && NumCores ~= 1
            parpool(NumCores);
        end
    catch
        ppool = matlabpool('size');
        if ~ppool
            matlabpool('local',NumCores); 
        end
    end
    M = NumCores;
    
    pctRunOnAll javaaddpath('java')
    ppm = ParforProgMon( 'Point Calibration ', npoints,1,400,50 );
%     profile on
    parfor (i=1:npoints,M)
        PCref = PCMinSinglePattern(Settings, ScanParams, Settings.CalibrationPointIndecies(i));
        disp(['Point: ' num2str(i)])
        CalibrationPointsPC(i,:) = PCref';
        ppm.increment();
    end
%     profile off
%     profile viewer
    ppm.delete();
    Settings.CalibrationPointsPC = CalibrationPointsPC;

    psize = Settings.PhosphorSize;
    
    handles.MeanXstar = mean(Settings.CalibrationPointsPC(:,1)+(Settings.XData(CalibrationPointIndecies))/psize);
    handles.MeanYstar = mean(Settings.CalibrationPointsPC(:,2)-(Settings.YData(CalibrationPointIndecies))/psize*sin(Settings.SampleTilt));
    handles.MeanZstar = mean(Settings.CalibrationPointsPC(:,3)-(Settings.YData(CalibrationPointIndecies))/psize*cos(Settings.SampleTilt));
%     disp(['xstar: ' num2str(handles.MeanXstar(1))]);
%     disp(['ystar: ' num2str(handles.MeanYstar(1))]);
%     disp(['zstar: ' num2str(handles.MeanZstar(1))]);

    handles.NaiveXstar = handles.MeanXstar-(Settings.XData)/psize;
    handles.NaiveYstar = handles.MeanYstar+(Settings.YData)/psize*sin(Settings.SampleTilt);
    handles.NaiveZstar = handles.MeanZstar+(Settings.YData)/psize*cos(Settings.SampleTilt);

    [n,V,p] = affine_fit(Settings.CalibrationPointsPC);

    line = zeros(npoints,6);
    for i=1:npoints
        line(i,:)=[Settings.CalibrationPointsPC(i,:) n'];
    end

    plane = [p V(:,1)' V(:,2)'];

    X=intersectLinePlane(line,plane);

    C = zeros(npoints,2);
    for i=1:npoints
        C(i,1)=Settings.XData(Settings.CalibrationPointIndecies(i));
        C(i,2)=Settings.YData(Settings.CalibrationPointIndecies(i));
    end

    A=zeros(3*npoints,9);
    PC=zeros(3*npoints,1);

    for i=1:npoints
        x=C(i,1);
        y=C(i,2);
        for k=1:3
            A(k+(i-1)*3,(1+3*(k-1)))=x;
            A(k+(i-1)*3,(2+3*(k-1)))=y;
            A(k+(i-1)*3,(3+3*(k-1)))=1;
            PC((i-1)*3+k)=X(i,k);
        end
    end

    coeffs=pinv(A)*PC;

    handles.FitXstar = coeffs(1)*Settings.XData+coeffs(2)*Settings.YData+coeffs(3);
    handles.FitYstar = coeffs(4)*Settings.XData+coeffs(5)*Settings.YData+coeffs(6);
    handles.FitZstar = coeffs(7)*Settings.XData+coeffs(8)*Settings.YData+coeffs(9);

    Settings.PhosphorSize = psize;
    handles.Settings = Settings;
    guidata(hObject, handles);
    
    cla(handles.axes2)
    handles.calibrated = 1;
    
    planefitpanel_SelectionChangeFcn(handles.planefitpanel, eventdata, handles);
    
    handles = guidata(hObject);
    if get(handles.autorunbox,'Value')
        savenclose_Callback(handles.savenclose,eventdata,handles);
    end
    
else
    mboxhandle = msgbox('You need to select points to calibrate first');
    pause(1)
    close(mboxhandle)
    
end


% --- Executes on button press in selectpointsbutton.
function selectpointsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectpointsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Right-click last point or use RETURN key to exit.

Settings = handles.Settings;
Nx = Settings.Nx;
Ny = Settings.Ny;
IQ = Settings.IQ;
XData = Settings.XData;
YData = Settings.YData;
NumColsEven = floor(Nx/2);
NumColsOdd = ceil(Nx/2);
IQPlot = handles.IQPlot;
XStep = XData(2)-XData(1);
YStep = YData(YData > 0);
if isempty(YStep) %Lines Scans
    YStep = 0;
else
    YStep = YStep(1);
end
Title = 'Press RETURN key or right-click last point to exit';
if Settings.Ny > 1
    MinPoints = 3;
elseif Settings.Ny == 1
    MinPoints = 1;
end

%Create Correct Indice Matrix
switch Settings.ScanType
    case 'Square'
        indi = 1:1:Settings.Nx*Settings.Ny;
        indi = reshape(indi, Settings.Nx,Settings.Ny)';
        
        %Create the Figure
        selectfig = figure;
        imagesc(IQPlot)
        axis image
        title(Title);
        
        npoints = 1;
        morepoints = 1;
        while morepoints   
            %Gets X,Y data from user
            [x,y, button] = ginput(1);
            if x > Nx
                x = Nx;
            elseif x < 1
                x = 1;
            end
            if y > Ny
                y = Ny;
            elseif y < 1
                y = 1;
            end
            
            if ~isempty(x)
                CalibrationPointIndecies(npoints) = indi(round(y), round(x));
                Xind(npoints) = round(x);
                Yind(npoints) = round(y);
                
                if button ~= 1 && npoints >= MinPoints
                    morepoints = 0;
                end
                npoints = npoints + 1;
            elseif npoints > MinPoints %When RETURN key is pressed
                morepoints = 0;
            end
            
            hold on
            plot(round(x),round(y),'kd','MarkerFaceColor','k');
        end
        close (selectfig);
        hold off
        axes(handles.axes1)
        imagesc(IQPlot)
        axis image
        hold on
        plot(Xind,Yind,'kd','MarkerFaceColor','k');
        %CalibrationPointIndecies
        
    case 'Hexagonal'
        indi = 1:length(XData);
        indi = Hex2Array(indi, NumColsOdd, NumColsEven);
        color = IQPlot.iq;
        
        selectfig = figure;
                
        npoints = 1;
        morepoints = 1;
        while morepoints
            %Create the Figure
            surf(IQPlot.x,IQPlot.y,IQPlot.iq,color,'EdgeColor','none');
            view(2);
            caxis(IQPlot.Limits);
            title(Title);
            colormap jet
            
            %Gets X,Y data from user
            [x,y, button] = ginput(1);
            Xind = round(x/XStep);
            Yind = round(y/YStep);
            if Xind > NumColsEven
                Xind = NumColsEven;
            elseif Xind < 1
                Xind = 1;
            end
            if Yind > Ny
                Yind = Ny;
            elseif Yind < 1
                Yind = 1;
            end
            
            if ~isempty(x)
            
                CalibrationPointIndecies(npoints) = indi(Yind, Xind);

                hold on
                scatter3(XData(CalibrationPointIndecies(npoints)),YData(CalibrationPointIndecies(npoints)),max(IQ));
                color(Yind,Xind) = max(IQ)*1000;

                if button ~= 1 && npoints > 2
                    morepoints = 0;
                end
                npoints = npoints + 1;
            elseif npoints > 3 %When RETRUN key is pressed
                morepoints = 0;
            end    
        end
        close (selectfig);
        hold off
        axes(handles.axes1)
        surf(IQPlot.x,IQPlot.y,IQPlot.iq,'EdgeColor','none');
        view(2);
        caxis(IQPlot.Limits);
        hold on
        height(1:length(CalibrationPointIndecies)) = max(IQ);
        scatter3(XData(CalibrationPointIndecies),YData(CalibrationPointIndecies),height);
                
end
Settings.CalibrationPointIndecies = CalibrationPointIndecies;
handles.Settings = Settings;
handles.VanPont = 1;
guidata(hObject, handles);

% --- Executes on button press in autorunbox.
function autorunbox_Callback(hObject, eventdata, handles)
% hObject    handle to autorunbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of autorunbox

% --------------------------------------------------------------------
function planefitpanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to planefitpanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Settings = handles.Settings;
Nx = Settings.Nx;
NumColsEven = floor(Nx/2);
NumColsOdd = ceil(Nx/2);
ScanType = Settings.ScanType;

if handles.calibrated

    cla(handles.axes2)
    
    % Naive Plane Fit
    if get(handles.naivebutton,'Value')
        Settings.XStar = handles.NaiveXstar;
        Settings.YStar = handles.NaiveYstar;
        Settings.ZStar = handles.NaiveZstar;

        max(handles.NaiveXstar);
        min(handles.NaiveXstar);

        axes(handles.axes2)
        plot3(handles.ScanParams.xstar,handles.ScanParams.ystar,handles.ScanParams.zstar,'bo')
        hold on
        plot3(Settings.CalibrationPointsPC(:,1),Settings.CalibrationPointsPC(:,2),Settings.CalibrationPointsPC(:,3),'ro')
        
        switch ScanType
            case 'Square'
                NaiveXstar = reshape(handles.NaiveXstar,Settings.Nx,Settings.Ny)';
                NaiveYstar = reshape(handles.NaiveYstar,Settings.Nx,Settings.Ny)';
                NaiveZstar = reshape(handles.NaiveZstar,Settings.Nx,Settings.Ny)';
            case 'Hexagonal'
                NaiveXstar = Hex2Array(handles.NaiveXstar(1:end-1), NumColsOdd, NumColsEven);
                NaiveYstar = Hex2Array(handles.NaiveYstar(1:end-1), NumColsOdd, NumColsEven);
                NaiveZstar = Hex2Array(handles.NaiveZstar(1:end-1), NumColsOdd, NumColsEven);
        end
        
        if Settings.Ny > 1
            surf(NaiveXstar,NaiveYstar,NaiveZstar,zeros(size(NaiveZstar)))
        else
            plot3(NaiveXstar,NaiveYstar,NaiveZstar)
        end
        %surf(reshape(handles.NaiveXstar,Settings.Nx,Settings.Ny)',reshape(handles.NaiveYstar,Settings.Nx,Settings.Ny)',reshape(handles.NaiveZstar,Settings.Nx,Settings.Ny)',zeros(Settings.Ny,Settings.Nx))
        shading flat
    end
    
    % PC Data Fit
    if get(handles.pcplanefit,'Value')
        Settings.XStar = handles.FitXstar;
        Settings.YStar = handles.FitYstar;
        Settings.ZStar = handles.FitZstar;

        axes(handles.axes2)
        plot3(handles.ScanParams.xstar,handles.ScanParams.ystar,handles.ScanParams.zstar,'bo')
        hold on
        plot3(Settings.CalibrationPointsPC(:,1),Settings.CalibrationPointsPC(:,2),Settings.CalibrationPointsPC(:,3),'ro')
        switch ScanType
            case 'Square'
                FitXstar = reshape(handles.FitXstar,Settings.Nx,Settings.Ny)';
                FitYstar = reshape(handles.FitYstar,Settings.Nx,Settings.Ny)';
                FitZstar = reshape(handles.FitZstar,Settings.Nx,Settings.Ny)';
            case 'Hexagonal'
                FitXstar = Hex2Array(handles.FitXstar, NumColsOdd, NumColsEven);
                FitYstar = Hex2Array(handles.FitYstar, NumColsOdd, NumColsEven);
                FitZstar = Hex2Array(handles.FitZstar, NumColsOdd, NumColsEven);
                FitXstar = FitXstar(1:length(FitXstar)-1,:);
                FitYstar = FitYstar(1:length(FitYstar)-1,:);
                FitZstar = FitZstar(1:length(FitZstar)-1,:);
        end
        if Settings.Ny > 1
            surf(FitXstar,FitYstar,FitZstar,.5*ones(size(FitXstar)));
        else
            plot3(FitXstar,FitYstar,FitZstar);
        end
        shading flat
    end

    if get(handles.nofit,'Value')
        Settings.XStar = handles.MeanXstar*ones(size(handles.NaiveXstar));
        Settings.YStar = handles.MeanYstar*ones(size(handles.NaiveXstar));
        Settings.ZStar = handles.MeanZstar*ones(size(handles.NaiveXstar));

        axes(handles.axes2)
        plot3(handles.ScanParams.xstar,handles.ScanParams.ystar,handles.ScanParams.zstar,'bo')
        hold on
        plot3(Settings.CalibrationPointsPC(:,1),Settings.CalibrationPointsPC(:,2),Settings.CalibrationPointsPC(:,3),'ro')
        plot3(handles.MeanXstar,handles.MeanYstar,handles.MeanZstar,'go')
    end
    

    handles.Settings = Settings;
    guidata(hObject, handles);

end

if handles.tiffread && get(handles.fromtiff,'Value')
    Settings.Xstar = handles.TiffXstar;
    Settings.Ystar = handles.TiffYstar;
    Settings.Zstar = handles.TiffZstar;
    
    axes(handles.axes2)
    plot3(handles.ScanParams.xstar,handles.ScanParams.ystar,handles.ScanParams.zstar,'bo')
    hold on
    try
        plot3(Settings.CalibrationPointsPC(:,1),Settings.CalibrationPointsPC(:,2),Settings.CalibrationPointsPC(:,3),'ro')
    catch me
    end
    
    TiffXstar = reshape(handles.TiffXstar,Settings.Nx,Settings.Ny)';
    TiffYstar = reshape(handles.TiffYstar,Settings.Nx,Settings.Ny)';
    TiffZstar = reshape(handles.TiffZstar,Settings.Nx,Settings.Ny)';
    
    if Settings.Ny > 1
        surf(TiffXstar,TiffYstar,TiffZstar,zeros(size(TiffZstar)))
    else
        plot3(TiffXstar,TiffYstar,TiffZstar)
    end
    shading flat
    
    handles.Settings = Settings;
    guidata(hObject, handles);  
    
end

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on button press in LoadPCCalc.
function LoadPCCalc_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPCCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
w = cd;
Settings = handles.Settings;
cd(fileparts(Settings.ScanFilePath));
[filename, filepath] = uigetfile('*mat','Select an Analysis_Params file');
cd(w);
matfile = load(fullfile(filepath,filename));
if isfield(matfile,'Settings')
    loadedSettings = matfile.Settings;
    valid = 1;
    if (Settings.AccelVoltage ~= loadedSettings.AccelVoltage)...
        ||(Settings.SampleTilt ~= loadedSettings.SampleTilt)...
        ||(Settings.SampleAzimuthal ~= loadedSettings.SampleAzimuthal)...
        ||(Settings.CameraElevation ~= loadedSettings.CameraElevation)...
        ||(Settings.CameraAzimuthal ~= loadedSettings.CameraAzimuthal)...
        ||(~strcmp(Settings.ScanFilePath,loadedSettings.ScanFilePath))...
        ||(~strcmp(Settings.ScanType,loadedSettings.ScanType))...
        ||(sum(Settings.IQ ~= loadedSettings.IQ)>0)
        valid = 0;
    end
    if valid == 0
        warndlg('Selected Scan not compatible with current scan');
    end
    % Validation
    if isfield(loadedSettings,'PCCal') && valid
        PCCal = loadedSettings.PCCal;
        handles.MeanXstar = PCCal.MeanXstar;
        handles.MeanYstar = PCCal.MeanYstar;
        handles.MeanZstar = PCCal.MeanZstar;
        handles.FitXstar = PCCal.FitXstar;
        handles.FitYstar = PCCal.FitYstar;
        handles.FitZstar = PCCal.FitZstar;
        handles.NaiveXstar = PCCal.NaiveXstar;
        handles.NaiveYstar = PCCal.NaiveYstar;
        handles.NaiveZstar = PCCal.NaiveZstar;
        handles.calibrated = 1;
        handles.Settings.CalibrationPointIndecies = loadedSettings.CalibrationPointIndecies;
        handles.Settings.CalibrationPointsPC = loadedSettings.CalibrationPointsPC;
        handles.Settings.PhosphorSize = loadedSettings.PhosphorSize;
        
        IQPlot = handles.IQPlot;
        CalibrationPointIndecies = loadedSettings.CalibrationPointIndecies;
        Nx = Settings.Nx;
        Ny = Settings.Ny;
        IQ = Settings.IQ;
        XData = Settings.XData;
        YData = Settings.YData;
        NumColsEven = floor(Nx/2);
        NumColsOdd = ceil(Nx/2);
        switch Settings.ScanType
            case 'Square'
                indi = 1:1:Nx*Ny;
                indi = reshape(indi,Nx,Ny)';
                for i = 1:length(CalibrationPointIndecies)
                    Xind(i) = mod(CalibrationPointIndecies(i),Nx);
                    Yind(i) = ceil(CalibrationPointIndecies(i)/Nx);
                end
                axes(handles.axes1)
                imagesc(IQPlot)
                axis image
                hold on
                plot(Xind,Yind,'kd','MarkerFaceColor','k');
            case 'Hexagonal'
                indi = 1:length(XData);
                indi = Hex2Array(indi, NumColsOdd, NumColsEven);
                axes(handles.axes1)
                surf(IQPlot.x,IQPlot.y,IQPlot.iq,'EdgeColor','none');
                view(2);
                caxis(IQPlot.Limits);
                hold on
                height(1:length(CalibrationPointIndecies)) = max(IQ);
                scatter3(XData(CalibrationPointIndecies),YData(CalibrationPointIndecies),height);
        end
        guidata(hObject,handles);
    end
end
clear matfile;


% --- Executes on button press in SavePCCalc.
function SavePCCalc_Callback(hObject, eventdata, handles)
% hObject    handle to SavePCCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[path, name] = fileparts(handles.Settings.OutputPath);
Settings = handles.Settings;
Settings.PCCal.MeanXstar = handles.MeanXstar;
Settings.PCCal.MeanYstar = handles.MeanYstar;
Settings.PCCal.MeanZstar = handles.MeanZstar;
Settings.PCCal.FitXstar = handles.FitXstar;
Settings.PCCal.FitYstar = handles.FitYstar;
Settings.PCCal.FitZstar = handles.FitZstar;
Settings.PCCal.NaiveXstar = handles.NaiveXstar;
Settings.PCCal.NaiveYstar = handles.NaiveYstar;
Settings.PCCal.NaiveZstar = handles.NaiveZstar;

save(fullfile(path,name),'Settings');


% --- Executes on button press in SavePCCal.
function SavePCCal_Callback(hObject, eventdata, handles)
% hObject    handle to SavePCCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SavePCCal
handles.SaveAllPC = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in fromtiffbutton.
function fromtiffbutton_Callback(hObject, eventdata, handles)
% hObject    handle to fromtiffbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
Settings = handles.Settings;
info = imfinfo(Settings.ImageNamesList{1});

if Settings.ImageTag %See MainGUI.m SetImageFields

    handles.TiffXstar = zeros(size(Settings.XData));
    handles.TiffYstar = handles.TiffXstar;
    handles.TiffZstar = handles.TiffXstar;
    
    VHRatio = Settings.VHRatio;
    lscan = Settings.ScanLength;

    h = waitbar(0,'Reading Tiff Files');
    for i=1:lscan
        try
            info = imfinfo(Settings.ImageNamesList{i});

            xistart = strfind(info.UnknownTags.Value,'<pattern-center-x-pu>');
            xifinish = strfind(info.UnknownTags.Value,'</pattern-center-x-pu>');

            thisx = str2double(info.UnknownTags.Value(xistart+length('<pattern-center-x-pu>'):xifinish-1));
            handles.TiffXstar(i) = (thisx - (1-VHRatio)/2)/VHRatio;

            yistart = strfind(info.UnknownTags.Value,'<pattern-center-y-pu>');
            yifinish = strfind(info.UnknownTags.Value,'</pattern-center-y-pu>');

            handles.TiffYstar(i) = str2double(info.UnknownTags.Value(yistart+length('<pattern-center-y-pu>'):yifinish-1));

            zistart = strfind(info.UnknownTags.Value,'<detector-distance-pu>');
            zifinish = strfind(info.UnknownTags.Value,'</detector-distance-pu>');

            handles.TiffZstar(i) = str2double(info.UnknownTags.Value(zistart+length('<detector-distance-pu>'):zifinish-1))/VHRatio;
        catch
            handles.TiffXstar(i) = handles.ScanParams.xstar;
            handles.TiffYstar(i) = handles.ScanParams.ystar;
            handles.TiffZstar(i) = handles.ScanParams.zstar;
        end
        waitbar(i/lscan,h)
    end
    close(h)

    handles.tiffread = 1;
    set(handles.fromtiff,'Enable','on');

    if ~handles.calibrated
        Settings.Xstar = handles.TiffXstar;
        Settings.Ystar = handles.TiffYstar;
        Settings.Zstar = handles.TiffZstar;
        
        handles.Settings = Settings;
    end
    
    guidata(hObject, handles);
    if get(handles.autorunbox,'Value')
        set(handles.fromtiff,'Value',1)
    end
    planefitpanel_SelectionChangeFcn(handles.planefitpanel, eventdata, handles);
    
    if get(handles.autorunbox,'Value')
        savenclose_Callback(handles.savenclose,eventdata,handles);
    end

else
    msgbox('No PC data found in image files','Read PC from TIFF')
end

