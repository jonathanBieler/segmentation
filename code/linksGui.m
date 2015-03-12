function varargout = linksGui(varargin)
% LINKSGUI MATLAB code for linksGui.fig
%      LINKSGUI, by itself, creates a new LINKSGUI or raises the existing
%      singleton*.
%
%      H = LINKSGUI returns the handle to a new LINKSGUI or the handle to
%      the existing singleton*.
%
%      LINKSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LINKSGUI.M with the given input arguments.
%
%      LINKSGUI('Property','Value',...) creates a new LINKSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before linksGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to linksGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help linksGui

% Last Modified by GUIDE v2.5 14-Jan-2015 15:38:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @linksGui_OpeningFcn, ...
                   'gui_OutputFcn',  @linksGui_OutputFcn, ...
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


% --- Executes just before linksGui is made visible.
function linksGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to linksGui (see VARARGIN)

    % Choose default command line output for linksGui
    handles.output = hObject;
    
    disp(handles)

    % Update handles structure
    guidata(hObject, handles);
    
    
    set(hObject,'Interruptible','off');
    
    global axesData ;
    axesData = handles.axes1;

    global axesSeg ;
    axesSeg = handles.axes2;


    global aKey dKey key1 key2;
    global numMouseMove canChange;
    numMouseMove = 0;
    canChange =1;
    
    global k t;
    k=1;
    t=0;
    
    global moviePath Nframes;     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% MAIN PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    moviePath = [pwd() '/'] %'/Users/bieler/Desktop/matlab/13sept13_longdaysin1to5uM/movie34/';
    
    
    if(exist([moviePath 'zStackedThreshCorrected/1.png'],'file'))
        seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
    else
        copyFileButton_Callback(hObject, eventdata, handles);
        seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
        
    end
    
    
    global a; 
    
    a = getImage(k,handles);
        
    axes(axesSeg)
    segImage = imagesc(seg);
    
    set(axesSeg, 'xlimmode','manual',...
    'ylimmode','manual',...
    'zlimmode','manual',...
    'climmode','manual',...
    'alimmode','manual');
    
    %save all data
    handles = guidata(gcf);
    
    handles.axesSeg=axesSeg;
    handles.dKey=0;
    handles.previousMouseX=0; 
    handles.previousMouseY=0;
    handles.seg = seg;    
    handles.segImage = segImage;
 
    guidata(gcf,handles);
    
    global links nuclei;
    
    if( ~exist('links.mat','file') )
        handles.seg = 0;
        doLinksOnly = 1;
        
        Nframes = length(dir('zStackedThreshSplit/*.png'));
        
        disp('Doing the initial tracking')
        updateTracking_Callback(0, eventdata, handles,doLinksOnly);    
        %edit4_Callback(hObject, eventdata, handles)
    else
        
        load links.mat;
        load nuclei.mat;
    end


    aKey=0;
    dKey=0;
    
    key1 = 1;
    key2 = 1;
    

    Nframes = length(links);
    
    %moviePath =  get(handles.edit1,'string');    
        
    set(handles.edit1,'string',moviePath);        
    set(handles.edit2,'string',Nframes);    
    
    %Nframes = str2num( get(handles.edit2,'string') );

    
    updateImageData(0,0);
    updateImageSeg(seg);

    global previousMouseX previousMouseY;
    previousMouseX=0;
    previousMouseY=0;
    
            
    set(gcf,'doublebuffer','off');
    set(gcf, 'WindowButtonMotionFcn', @mouseMove);
    set(gcf, 'WindowKeyPressFcn', @KeyPress)
    set(gcf, 'WindowScrollWheelFcn', @MouseScroll)
    
    set(gcf,'WindowButtonDownFcn',@MouseClick);
        
    set(gcf, 'WindowKeyReleaseFcn', @KeyRelease)

    % UIWAIT makes linksGui wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    

function a = getImage(k,handles)
    
    global moviePath;

    expe = experimentPara();
    
    if( expe.hasTrans && get(handles.transCheckBox,'value')  )
        a = abs( imread([moviePath 'img/' getImageName(expe.transName,k)]) );
        a = repmat(a,[1 1 3]);
    else
        a = imread([moviePath 'zStackedYFP/' num2str(k) '.png']);
    end
    



% --- Outputs from this function are returned to the command line.
function varargout = linksGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

    
function updateImageSeg(seg,previousSeg)
    
    if nargin ==1
        previousSeg = seg;
    end

    handles = guidata(gcf);
    %%seg = handles.seg;
    axesSeg = handles.axesSeg;
    
    
    axes(axesSeg);
    hold off
    imagesc( double(seg)+0.2*double(previousSeg) );
    %set(handles.segImage,'Cdata',double(seg)+0.2*double(previousSeg))
    %set(handles.segImage,'hold','on');
    hold on
    
    global links k nuclei;
    
    if( get(handles.dispLinks,'value') )
    
       
        if(~isempty(links)) 
            tmp = links{k};    

            for i=1:size(tmp,1)

                x1 = nuclei(k).carth( tmp(i,1),: );                
                x2 = nuclei(tmp(i,3)).carth( tmp(i,2),: );

                plot([x1(1) x2(1)],[x1(2) x2(2)],'w')
                plot([x1(1) x2(1)],[x1(2) x2(2)],'w.')


            end
        end
    end
    
    %imagesc(seg)
    %set(gca,'Ydir','normal');
    %axis([1 512 1 512])
    

function updateImageData(x,y)
    
    global a;
    global axesData;

    axes(axesData);
    cla;
    hold on;
    d = imnorm(a).^0.5;
    imagesc(d)
    d=d(:);
    
    caxis([quantile(d,0.01),quantile(d,0.99)])

    plot(x,y,'w+')
    
    N1 = size(a,1);
    N2 = size(a,2);

    set(gca,'Ydir','reverse');
    axis([1 N2 1 N1])
    
    
function KeyRelease(Source,EventData)

    global key1 key2
    
    switch EventData.Character
       
        case '1'
            key1=1;
            
        case '2'
            key2=1;
        
    end
    
        
    
function KeyPress(Source, EventData)
    

global k a previousSeg dKey sKey aKey moviePath Nframes canChange key1 key2;   

handles = guidata(gcf);
seg = handles.seg;
    
EventData.Character

N1 = size(seg,1);
N2 = size(seg,2);


  
switch EventData.Character
    
    case ' '
        
        set(handles.transCheckBox,'value', ~get(handles.transCheckBox,'value') );
        
        
        C = get(handles.axesSeg, 'CurrentPoint');
        i = ceil(C(1,1));
        j = ceil(C(1,2));

        i = max(i,1); j = max(j,1);             
        i = min(i,512); j = min(j,512);
        
        a = getImage(k,handles);
        updateImageData(i,j);
    
    case 'u'
        
        updateTracking_Callback(0, EventData, handles);
    
    case '+'
        
        th = str2double( get(handles.edit3,'String') );
        th = th+0.05; 
        set(handles.edit3,'String',num2str(th));
    
    case '-'
        
        th = str2double( get(handles.edit3,'String') );
        th = th-0.05; 
        set(handles.edit3,'String',num2str(th));
        
        
    %split 
    case 'f'
        
        
       C = get(handles.axesSeg, 'CurrentPoint');
       i = ceil(C(1,1));
       j = ceil(C(1,2));
        
       i = max(i,1); j = max(j,1);             
       i = min(i,512); j = min(j,512);
       
       [seli selj] = getNeiInd(j,i,round(N1/8) ,N1,N2);
       
       sub_seg = imnorm( seg(seli,selj) );
              
       seg(seli,selj) = split(sub_seg,handles);
       
       if(k>1)
           previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
        else
           previousSeg = seg; 
        end    
             
       updateImageSeg(seg,previousSeg);
        
            
    case 'd'
       
       if(dKey == 1)
            dKey = 0;
           
            se = strel('disk',3);
            seg = imfill(seg,'holes');
            seg = imopen(seg, se);     

            if(k>1)
               previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
            else
               previousSeg = seg; 
            end    

            updateImageSeg(seg,previousSeg);
       else
           dKey = 1;
           aKey = 0;
       end    
       
    case 'a'
       
       if(aKey == 1)
           aKey = 0;
           
            se = strel('disk',3);
            seg = imfill(seg,'holes');
            seg = imopen(seg, se);        
            
            if(k>1)
               previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
            else
               previousSeg = seg; 
            end    

            updateImageSeg(seg,previousSeg);
                      
       else
           aKey = 1;
           dKey = 0;
       end   
       
   case 's'
       
       %imwrite(seg,[moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
       %disp(['saved frame:' num2str(k)])
       
       C = get(handles.axesSeg, 'CurrentPoint');
       i = ceil(C(1,1)); j = ceil(C(1,2));
               
       i = max(i,1); j = max(j,1);             
       i = min(i,N2); j = min(j,N1);
       
       [seli selj] = getNeiInd(j,i,round(N1/8),N2,N1);
              
       sub_seg = imnorm( a(seli,selj) );                   
       th = quantile(sub_seg(:),str2double( get(handles.edit3,'String')) );       
       sub_seg = im2bw(sub_seg,th);
       
              
       %circular mask
       x = linspace(-1,1,length(seli));
       [x y] = meshgrid(x,x); r = x.^2 + y.^2;
        
       out = (r>=1).*seg(seli,selj) + (r<1).*sub_seg;
        
       se = strel('disk',4);
       out = imfill(out,'holes'); out = imopen(out, se);
              
       seg(seli,selj) = out;
       
       
       if(k>1)
           previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
        else
           previousSeg = seg; 
        end    
             
       updateImageSeg(seg,previousSeg);       
       updateImageData(i,j);
       
    %segment three frames
    case 'q'

       
       for ii = 1:4
           
                   
           C = get(handles.axesSeg, 'CurrentPoint');
           i = ceil(C(1,1)); j = ceil(C(1,2));

           i = max(i,1); j = max(j,1);             
           i = min(i,N2); j = min(j,N1);
       
           [seli selj] = getNeiInd(j,i,round(N1/8),N1,N2);

           sub_seg = imnorm( a(seli,selj) );                   
           th = quantile(sub_seg(:),str2double( get(handles.edit3,'String')) );       
           sub_seg = im2bw(sub_seg,th);


           %circular mask
           x = linspace(-1,1,length(seli));
           [x y] = meshgrid(x,x); r = x.^2 + y.^2;

           out = (r>=1).*seg(seli,selj) + (r<1).*sub_seg;

           se = strel('disk',4);
           out = imfill(out,'holes'); out = imopen(out, se);

           seg(seli,selj) = out;


           if(k>1)
               previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
            else
               previousSeg = seg; 
            end    

           updateImageSeg(seg,previousSeg);       
           updateImageData(i,j);
           
           pause(0.25)
           
           %move to next frame
            imwrite(seg,[moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
            disp(['saved frame:' num2str(k)])

            if(k==Nframes)
                break;
            end
            
            k=min(k+1,Nframes);

            a = getImage(k,handles);
            seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);

            if(k>1)
               previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
            else
               previousSeg = seg; 
            end

            updateImageData(0,0);
            updateImageSeg(seg,previousSeg);
            pause(0.2)
          
       end
       
       
       %go back
       if(0)
        imwrite(seg,[moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
        disp(['saved frame:' num2str(k)])

       
        k=k-4;

        a = getImage(k,handles);
        seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);

        if(k>1)
           previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
        else
           previousSeg = seg; 
        end

        updateImageData(0,0);
        updateImageSeg(seg,previousSeg);
        pause(0.15)
       end
       
    case 'g' %draw circle
        
        
       C = get(handles.axesSeg, 'CurrentPoint');
       i = ceil(C(1,1));
       j = ceil(C(1,2));
        
       i = max(i,1); j = max(j,1);             
       i = min(i,N2); j = min(j,N1);
       
       [seli selj] = getNeiInd(j,i,round(N1/8),N2,N1);
              
       sub = ( seg(seli,selj) );
       
            
              
       %circular mask
        x = linspace(-1,1,length(seli));
        [x y] = meshgrid(x,x);
        r = sqrt(x.^2 + y.^2);
        
        sub(r<0.3)=1;
        
       
       seg(seli,selj) = sub;
       
       updateImageSeg(seg,previousSeg);
       
       %erase area
    case 'e'
              
       C = get(handles.axesSeg, 'CurrentPoint');
       i = ceil(C(1,1));
       j = ceil(C(1,2));
        
       i = max(i,1); j = max(j,1);             
       i = min(i,N2); j = min(j,N1);
       
       [seli selj] = getNeiInd(j,i,round(N1/8),N2,N1);
              
       seg(seli,selj) = 0;
       
       if(k>1)
           previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
        else
           previousSeg = seg; 
        end    
             
       updateImageSeg(seg,previousSeg);
       
    %erase object
    case 'r' 
        
       C = get(handles.axesSeg, 'CurrentPoint');
       i = ceil(C(1,1));
       j = ceil(C(1,2));
        
       i = max(i,1); j = max(j,1);             
       i = min(i,N2); j = min(j,N1);
       
       dist = max(min([ i, N2-i,81,j,N1-j]),9);
       
       [seli selj] = getNeiInd(j,i,dist,N2,N1);
       
       seli = unique(seli);
       selj = unique(selj); %remove mirror on borders
       
        m = imnorm( seg(seli,selj) );
                            
        labeledImage = bwlabel(m, 8);    
        val = labeledImage(round(size(m,1)/2),round(size(m,2)/2));
        m(  labeledImage== val)=0;
                
        seg(seli,selj) = m;
        
        if(k>1)
           previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
        else
           previousSeg = seg; 
        end    
             
       updateImageSeg(seg,previousSeg);
       
       
    %do a bit of smoothing    
    case 'w'
               
        se = strel('disk',3);
        seg = imfill(seg,'holes');
        seg = imopen(seg, se);        
        %seg = imclose(seg, se);
        
        %H = fspecial('gaussian',50,1);
        %seg = imfilter(double(seg),H,'symmetric') > 0.5;

        guidata(gcf,handles);
        updateImageSeg(seg);

    
    case '1'
        
        if(key1==1)
        
            key1=0;
            imwrite(seg,[moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);

            %uiresume(gcf)
            disp(['saved frame:' num2str(k)])

            k = max(k-1,1);

            a = getImage(k,handles);
            seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);

            if(k>1)
               previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
            else
               previousSeg = seg; 
            end        

            updateImageData(0,0);
            updateImageSeg(seg,previousSeg);   
        end

        
    case '2'
        
        if(key2==1)
            
            
            
            key2=0;
        
            imwrite(seg,[moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
            disp(['saved frame:' num2str(k)])

            k=min(k+1,Nframes);

            a = getImage(k,handles);
            seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);

            if(k>1)
               previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
            else
               previousSeg = seg; 
            end

            updateImageData(0,0);
            updateImageSeg(seg,previousSeg);
            

        end

end

    %update crosshair
    if(dKey == 1)
        set(gcf,'Pointer','cross');
    elseif(aKey==1)
        set(gcf,'Pointer','crosshair');
    else
        set(gcf,'Pointer','arrow');
    end


    handles = guidata(gcf);
        
    C = get(handles.axesSeg, 'CurrentPoint');
    handles.previousMouseX=C(1,1); 
    handles.previousMouseY=C(1,2);
        

    handles.dKey=dKey;
    handles.seg = seg;
    
    guidata(gcf,handles);
    
    
    set(handles.text5,'String',num2str(k));
    
function MouseScroll(hObject,eventdata)
    
    handles = guidata(gcf);
    
    dTh = eventdata.VerticalScrollCount;
    
    th = str2double( get(handles.edit3,'String') );
    th = th - 0.025*dTh; 
    th = min(th,1);
    th = max(th,0);
    set(handles.edit3,'String',num2str(th));
    
function MouseClick(hObject,eventdata)


    EventData = struct();
    EventData.Character = 's';

    KeyPress(0, EventData)


function mouseMove(hObject,eventdata)

   
    %global axesSeg dKey previousMouseX previousMouseY;
    global aKey axesData k a moviePath;
    global numMouseMove Nframes;
    numMouseMove = numMouseMove+1;
    if( numMouseMove > 2 ||  (aKey == 1) )
        
    numMouseMove = 0;
        
    %load data
    handles = guidata(gcf);
    
    seg = handles.seg;
    axesSeg = handles.axesSeg;
    dKey = handles.dKey;
    previousMouseX = handles.previousMouseX;
    previousMouseY = handles.previousMouseY;
    
    
    N1 = size(seg,1);
    N2 = size(seg,2);
    
    %%change cursor
    if(dKey == 1)
        set(gcf,'Pointer','cross');
    elseif(aKey==1)
        set(gcf,'Pointer','crosshair');
    else
        set(gcf,'Pointer','arrow');
    end
        
    C = get (axesSeg, 'CurrentPoint');
    %title(axesSeg, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);

    mouseX = C(1,1);
    mouseY = C(1,2);
    
    % move data image in time
    C2 = get (axesData, 'CurrentPoint');
    mouseX_data = ceil( C2(1,1));
    mouseY_data = ceil( C2(1,2) );
    %
    if(mouseX_data>0 && mouseX_data<512 && mouseY_data >0 && mouseY_data < 512 )
       
        kk = round(20*mouseX_data/512);
        
        dataIx = min(k + kk,Nframes);
        
        a = getImage(dataIx,handles);
                
        updateImageData(0,0);
        
        
    end
    
    
    
    if(dKey == 1 || aKey == 1)
       
        i = ceil(mouseX);
        j = ceil(mouseY);
        
        i = max(i,1); j = max(j,1);             
        i = min(i,N2); j = min(j,N1);
        
        ip = ceil(previousMouseX);
        jp = ceil(previousMouseY);
        
        previousMouseX = mouseX;
        previousMouseY = mouseY;
        
        ip = max(ip,1); jp = max(jp,1);             
        ip = min(ip,N2); jp = min(jp,N1);
        
        dist = sqrt( (ip-i)^2 + (jp-j)^2);
                
        line_i = linspace(ip,i,25*dist);
        line_j = linspace(jp,j,25*dist);
        
        line_i = line_i + 0.5*randn(size(line_i));
        line_j = line_j + 0.5*randn(size(line_j));
        
        line_i = ceil(line_i);
        line_j = ceil(line_j);
        
        line_i(line_i<1) = 1; line_i(line_i>N2) = N2;
        line_j(line_j<1) = 1; line_j(line_j>N1) = N1;
                
        ind = sub2ind([N1 N2],line_j,line_i);
                
        if(aKey)
            seg(ind)=0;
        end
        if(dKey)
            seg(ind)=1;
        end
                
        handles.seg = seg;
        guidata(gcf,handles);        
                
        updateImageSeg(seg);
        
    else
        
        i = ceil(mouseX);
        j = ceil(mouseY);
        
        if(rand <0.6)
            updateImageData(i,j);
        end
        
    end
    
    handles.previousMouseX=previousMouseX; 
    handles.previousMouseY=previousMouseY;
    
    guidata(gcf,handles);
    
    end
    
%    updateImageData();



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

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


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

% --- Executes on button press in copyFileButton.
function copyFileButton_Callback(hObject, eventdata, handles)

    global moviePath Nframes; 
        
    %Nframes = str2num( get(handles.edit2,'string') );

    mkdirIfNotExist([moviePath 'zStackedThreshCorrected/']);
        
    system(['cp -v ' moviePath 'zStackedThreshSplit/*.png ' moviePath 'zStackedThreshCorrected/'  ]);
    



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)

    global moviePath Nframes k; 

    system(['cp -v ' moviePath 'zStackedThreshSplit/' num2str(k) '.png ' moviePath 'zStackedThreshCorrected/'  ]);
    
    handles = guidata(gcf);
    
    seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
    handles.seg = seg;    
    guidata(gcf,handles);
    updateImageSeg(seg);
    

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

    global k moviePath;
    
    k = str2double( get(handles.edit4,'string') );
    

    seg = imread([moviePath 'zStackedThreshCorrected/' num2str(k) '.png']);
    handles.seg = seg;
    
    updateImageSeg(seg);
    guidata(gcf,handles);



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

   
function mOri = split(m,handles)

                        
%        m = m(selj,seli);   
        mOri=m;

        %is there more than 2 objects?
        labeledImage = bwlabel(m, 8);    
        val = labeledImage(round(size(m,1)/2),round(size(m,2)/2));
        m( labeledImage>0 & labeledImage~= val)=0;

        if( length(unique(labeledImage)) >0 )
                 
            [y x] = find(m==1);
            
            center=zeros(1,2);
            
            center(1) = mean(x);
            center(2) = mean(y);
                
            %m = (edge(m,'log',0));
            m = bwperim(m);
                    
            [y x] = find(m==1);
                        
            [thO rO]= cart2pol(x-center(1),y-center(2)); 
            [v ind] = sort(thO); thO=thO(ind); rO=rO(ind);
            
            chull = convhull(x,y);
            x = x(chull);   y = y(chull);
            
            x = interp1(linspace(0,1,length(x)),x,linspace(0,1,length(thO)));
            y = interp1(linspace(0,1,length(y)),y,linspace(0,1,length(thO)));
                        
            %go into polar coordinates                                    
            [th r]= cart2pol(x-center(1),y-center(2));             
            [v ind] = sort(th); th=th(ind); r=r(ind);
            
            [th,ia] = unique(th); r = r(ia);
            
            rInterp = interp1(th,r,thO);
            
            nanVal = isnan(rInterp);
            rInterp(nanVal) = 0;
                                    
            err = ( abs((rInterp - rO)) ) ;
            
            err(nanVal)=0;
                        
            %look for specific pattern
            thPat = linspace(-2*pi,2*pi,200);
            pat = exp(- (thPat+0*pi/2).^2 /0.1 );% +exp(- (thPat -pi/2).^2 /0.3 );
            pat = pat.^2;
            pat = pat ./ sum(pat);
                      
            peaks = getPeakFromShape(err, pat,40,str2double( get(handles.edit5,'string') ),0)
                     
            if( length(peaks) >= 2)
                
                rSep = linspace(-8,8,100);
                sepx = [rO(peaks(1)).*cos(thO(peaks(1))); rO(peaks(2)).*cos(thO(peaks(2)))] + center(1);
                sepy = [rO(peaks(1)).*sin(thO(peaks(1))); rO(peaks(2)).*sin(thO(peaks(2)))] + center(2);
                
                line = drawline([sepy(1) sepx(1)],[sepy(2) sepx(2)],size(mOri));
                sep = zeros(size(mOri));
                sep(line) = 1;
                
                 se = strel('disk',1);
                 sep = imdilate(sep,se);
                 
                 mOri = mOri .* (sep <1);  
                
                se = strel('disk',3);
                mOri = imopen(mOri,se);
                

            end
            
        end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


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


% --- Executes on button press in dispLinks.
function dispLinks_Callback(hObject, eventdata, handles)
% hObject    handle to dispLinks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dispLinks


% --- Executes on button press in updateTracking.
function updateTracking_Callback(hObject, eventdata, handles, doLinksOnly)
% hObject    handle to updateTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if(nargin == 3)
        
        doLinksOnly = get(handles.doLinksOnly,'value');
        
    end
        

    global Nframes k moviePath;  
    seg = handles.seg;
    
       
    threshFolder = 'zStackedThreshCorrected/';
    
    expe = experimentPara();
    doMeasures(Nframes,threshFolder,expe);   
    
    Me = loadMeasures(Nframes);
    
    [tracks, signal, traj, ind, divisions] = doTracking(Nframes, Me,doLinksOnly);
    
    if(k>1)
       previousSeg = imread([moviePath 'zStackedThreshCorrected/' num2str(k-1) '.png']);
    else
       previousSeg = seg; 
    end
    
    global links nuclei;
    
    load links.mat;
    load nuclei.mat;
    
    updateImageSeg(seg,previousSeg)


% --- Executes on button press in doLinksOnly.
function doLinksOnly_Callback(hObject, eventdata, handles)
% hObject    handle to doLinksOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doLinksOnly


% --- Executes on button press in transCheckBox.
function transCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to transCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transCheckBox
