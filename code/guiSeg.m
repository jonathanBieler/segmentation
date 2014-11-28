function varargout = guiSeg(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiSeg_OpeningFcn, ...
                   'gui_OutputFcn',  @guiSeg_OutputFcn, ...
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


% --- Executes just before guiSeg is made visible.
function guiSeg_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for guiSeg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


global moviePath imInd Nframes segPara segParaFrames; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moviePath = [pwd() '/'] %'/Users/bieler/Desktop/matlab/13sept13_longdaysin1to5uM/movie35/';

imInd = 1;

cd(moviePath)

Nframes = length(dir('zStackedYFP/*.png'));

set(gcf, 'WindowKeyPressFcn', @KeyPress)

if(exist('segPara.mat','file'))
    load('segPara.mat');
else
    segPara = struct('thresh',{},'minSize',{},'maxSize',{},'openSize',{},'ratioThresh',{});

end

if(exist('segParaFrames.mat','file'))
    load('segParaFrames.mat');
else
    segParaFrames = [];
end


updateImage(handles)
disp(segParaFrames)

if( ismember( imInd,segParaFrames ) )

    paraInd = find(segParaFrames==imInd);

    setParameters(handles, paraInd );   
    
end

set(gcf,'WindowStyle','normal')


function KeyPress(Source, EventData)
    

global moviePath imInd Nframes segPara segParaFrames; 

handles = guidata(gcf);
   
  
switch EventData.Character
    
    case 'q'
        
            imInd = max(imInd-20,1);
                  
            if( ismember( imInd,segParaFrames ))
                setParameters(handles, find(segParaFrames==imInd) );
            end
            
            updateImage(handles)
            
            disp(segParaFrames);
                    

    case 'w'
        
            imInd = min(imInd+20,Nframes);
        
            
            if( ismember( imInd,segParaFrames ))
                setParameters(handles, find(segParaFrames==imInd) );
            end
            
            updateImage(handles)
            
    case 's'
            doSegment(handles);
        
end


set(handles.text6,'String',num2str(imInd));


function updateImage(handles)

    global imInd;

    if(exist(['zStackedThresh/' num2str(imInd) '.png'],'file'))

        a = imread( ['zStackedYFP/' num2str(imInd) '.png']);

        b = imread( ['zStackedThresh/' num2str(imInd) '.png'] );

        %superpose both        
        b = bwmorph(b,'remove');

        a(b==1) = min(a(:));

        axes(handles.axes1)
        imagesc(imnorm(a).^0.5)
        colormap jet


    end

function [] =  setParameters(handles,k)

    global moviePath imInd Nframes segPara segParaFrames; 
        
    set(handles.edit1,'String',num2str(segPara(k).thresh));
    set(handles.edit2,'String',num2str(segPara(k).minSize));
    set(handles.edit3,'String',num2str(segPara(k).maxSize));
    set(handles.edit4,'String',num2str(segPara(k).openSize));
    set(handles.edit5,'String',num2str(segPara(k).ratioThresh));
    
    
function para =  getParameters(handles)
    
    para = [];
        
    para.thresh = str2double( get( handles.edit1,'String') );
    para.minSize = str2double( get( handles.edit2,'String'));
    para.maxSize = str2double( get( handles.edit3,'String'));
    para.openSize = str2double( get( handles.edit4,'String'));
    para.ratioThresh = str2double( get(handles.edit5,'String'));


% --- Outputs from this function are returned to the command line.
function varargout = guiSeg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [] = doSegment(handles)

    global moviePath imInd Nframes segPara segParaFrames; 

    a = imread([moviePath 'zStackedYFP/' num2str(imInd) '.png']);

    inputFolder = 'zStackedYFP/';
        
    %get parameters
    if( ismember( imInd,segParaFrames ) )
        
        paraInd = find(segParaFrames==imInd);
        
        setParameters(handles, paraInd );        
                
    else
        
        para =  getParameters(handles);
        segParaFrames = [segParaFrames imInd]; 
        
        paraInd = find(segParaFrames==imInd);
        
        segPara(paraInd) = para;
        
    end
    
    para = segPara(paraInd);
    
    disp(segParaFrames)
    
        
    
    save('segPara.mat','segPara');
    save('segParaFrames.mat','segParaFrames');
    

    doDraw = 0;

    tic
    [filters openFilter] = generateFilters(para,doDraw);
    segmentImage(imInd,para,inputFolder,filters,openFilter,doDraw);
    toc

    b = imread( ['zStackedThresh/' num2str(imInd) '.png'] );

    %superpose both
    %b = (edge(b,'log',0));
    b = bwmorph(b,'remove');


    a(b==1) = min(a(:));

    axes(handles.axes1)
    imagesc(imnorm(a).^0.5)
    colormap jet


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

    doSegment(handles)
    

    
function onParaChange(handles)
        
        
   global moviePath imInd Nframes segPara segParaFrames; 

        
    %get parameters
    if( ismember( imInd,segParaFrames ) )
        
        paraInd = find(segParaFrames==imInd);
        
        para =  getParameters(handles);
        segPara(paraInd) = para;
                
    else
        
        para =  getParameters(handles);
        segParaFrames = [segParaFrames imInd]; 
        
        paraInd = find(segParaFrames==imInd);
        
        segPara(paraInd) = para;
        
    end
        
    disp(segParaFrames)
        


function edit1_Callback(hObject, eventdata, handles)

onParaChange(handles);

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



function edit2_Callback(hObject, eventdata, handles)

onParaChange(handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)

onParaChange(handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)

onParaChange(handles);

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



function edit5_Callback(hObject, eventdata, handles)

onParaChange(handles);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
