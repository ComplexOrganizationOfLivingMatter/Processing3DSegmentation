function varargout = window(varargin)
% WINDOW MATLAB code for window.fig
%      WINDOW, by itself, creates a new WINDOW or raises the existing
%      singleton*.
%
%      H = WINDOW returns the handle to a new WINDOW or the handle to
%      the existing singleton*.
%
%      WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINDOW.M with the given input arguments.
%
%      WINDOW('Property','Value',...) creates a new WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help window

% Last Modified by GUIDE v2.5 23-Mar-2020 15:47:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @window_OpeningFcn, ...
                   'gui_OutputFcn',  @window_OutputFcn, ...
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

% --- Executes just before window is made visible.
function window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to window (see VARARGIN)

% Choose default command line output for window
handles.output = hObject;

set(0, 'currentfigure', hObject); 

set(handles.missingApical,'string', strjoin(arrayfun(@num2str, getappdata(0, 'notFoundCellsApical'), 'UniformOutput', false), ', '));
set(handles.missingBasal,'string', strjoin(arrayfun(@num2str, getappdata(0, 'notFoundCellsBasal'), 'UniformOutput', false), ', '))

setappdata(0, 'labelledImageTemp', getappdata(0, 'labelledImage'));
setappdata(0, 'lumenImageTemp', getappdata(0, 'lumenImage'));
resizeImg = getappdata(0,'resizeImg');
labelledImage = getappdata(0, 'labelledImage');
originalSize = size(labelledImage);
sizeResized = originalSize * resizeImg;
sizeResized(3) = originalSize(3);

setappdata(0, 'labelledImageTemp_Resized', imresize3(labelledImage, sizeResized, 'nearest'));
setappdata(0, 'lumenImageTemp_Resized', imresize3(double(getappdata(0, 'lumenImage')), sizeResized, 'nearest')>0);

setappdata(0, 'selectedZ', 1);
setappdata(0, 'cellId', 1);
setappdata(0,'showAllCells',0);
setappdata(0,'windowListener',1);
setappdata(0, 'showBackground', 0)

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using window.
if strcmp(get(hObject,'Visible'),'off')
    showSelectedCell()
end

% UIWAIT makes window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = window_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roiMask = getappdata(0, 'roiMask');
progressBar = waitbar(0, 'Saving... Please wait', 'WindowStyle', 'modal');
if roiMask ~= -1
    delete(roiMask);
    roiMask = -1;
    setappdata(0, 'roiMask', roiMask);
    labelledImage = getappdata(0, 'labelledImageTemp');
    newCellRegion = getappdata(0, 'newCellRegion');
    selectCellId = getappdata(0, 'cellId');
    selectedZ = getappdata(0, 'selectedZ');
    lumenImage = getappdata(0, 'lumenImageTemp');
    
    if sum(newCellRegion(:)) > 0
        %newCellRegion = imresize(double(newCellRegion), [size(labelledImage, 1)  size(labelledImage, 2)], 'nearest')>0;
        insideGland = newCellRegion>-1;
        if getappdata(0, 'canModifyOutsideGland') == 0
            insideGland = insideGland & labelledImage(:,:,selectedZ) > 0;
        end
        %% The order is important here: Because the lumen is already 0 on the labelled image
        if getappdata(0, 'canModifyInsideLumen') == 1
            insideGland(lumenImage(:,:,selectedZ) == 1) = 1;
        else
            insideGland(lumenImage(:,:,selectedZ) == 1) = 0;
        end
        
        if selectCellId > 0
            imgActualZ = labelledImage(:, :, selectedZ);
            mostCommonIndex = mode(imgActualZ(newCellRegion));
            % First condition: we do not care if you are modifying a cell
            % and moving into the background. We do care if you are
            % replacing one cell for another one.
            % Second condition: Beware if you create two non-overlapping
            % areas. If you are creating a new cell, we do not care about
            % this.
            if (mostCommonIndex ~= 0 && mostCommonIndex ~= selectCellId) || (sum(imgActualZ(newCellRegion) == selectCellId) == 0 && selectCellId <= max(labelledImage(:)))
                answer = questdlg(['You are mostly replacing ', num2str(mostCommonIndex) , ' with ', num2str(selectCellId),'. Are you sure you want to proceed?'], ...
                    'Confirm', ...
                    'Yes','No', 'No');
                if strcmp(answer, 'No')
                    close(progressBar)
                    return
                end
            end
            if selectCellId <= max(labelledImage(:))
                [x, y] = find(newCellRegion & insideGland);
            else % Add cell
                [x, y] = find(newCellRegion);
                
                colours = getappdata(0, 'colours');
                newColours = colorcube(255);
                colours(end+1, :) = newColours(randi(255), :);
                setappdata(0, 'colours', colours);
            end
            
            newIndices = sub2ind(size(labelledImage), x, y, ones(length(x), 1)*selectedZ);
            labelledImage(newIndices) = selectCellId;
            if getappdata(0, 'canModifyInsideLumen') == 1
                lumenImage(newIndices) = 0;
            end
        else
            [x, y] = find(newCellRegion);
            newIndices = sub2ind(size(labelledImage), x, y, ones(length(x), 1)*selectedZ);
            labelledImage(newIndices) = selectCellId;
            lumenImage(newIndices) = 1;
            labelledImage(lumenImage>0) = 0;
        end
        setappdata(0, 'labelledImageTemp', labelledImage);
        setappdata(0, 'lumenImageTemp', lumenImage);

        updateResizedImage();
        pause(2);
    end
end
close(progressBar)
showSelectedCell();


% --- Executes during object creation, after setting all properties.
function tbCellId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbCellId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tbCellId_Callback(hObject, eventdata, handles)
% hObject    handle to tbZFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tbZFrame as text
%        str2double(get(hObject,'String')) returns contents of tbZFrame as a double
setappdata(0,'windowListener',1);

setappdata(0,'cellId',str2double(get(hObject,'String')));
showSelectedCell();


function tbZFrame_Callback(hObject, eventdata, handles)
% hObject    handle to tbZFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tbZFrame as text
%        str2double(get(hObject,'String')) returns contents of tbZFrame as a double
setappdata(0,'windowListener',1);
labelledImage = getappdata(0, 'labelledImageTemp');
newFrameValue = str2double(get(hObject,'String'));
if newFrameValue > 0 && newFrameValue <= size(labelledImage, 3)
    setappdata(0, 'selectedZ', newFrameValue);
    showSelectedCell();
end


% --- Executes during object creation, after setting all properties.
function tbZFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbZFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in insertROI.
function insertROI_Callback(hObject, eventdata, handles)
% hObject    handle to insertROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'windowListener',0);

roiMask = getappdata(0, 'roiMask');
if roiMask ~= -1
    delete(roiMask);
end
try
    roiMask = impoly(gca);
    newCellRegion = createMask(roiMask);
    setappdata(0,'roiMask', roiMask);
    setappdata(0,'newCellRegion', newCellRegion);
catch
    disp('ROI cancelled')
end
setappdata(0,'windowListener',1);


% --- Executes on button press in increaseID.
function increaseID_Callback(hObject, eventdata, handles)
% hObject    handle to increaseID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'windowListener',1);
newValue = getappdata(0, 'cellId')+1;
labelledImage = getappdata(0, 'labelledImageTemp_Resized');

if newValue <= max(labelledImage(:))
    setappdata(0, 'cellId', newValue);
    set(handles.tbCellId,'string',num2str(newValue));
    showSelectedCell();
end

% --- Executes on button press in decreaseID.
function decreaseID_Callback(hObject, eventdata, handles)
% hObject    handle to decreaseID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'windowListener',1);
newValue = getappdata(0, 'cellId')-1;
if newValue >= 0
    setappdata(0, 'cellId', newValue);
    set(handles.tbCellId,'string',num2str(newValue));
    showSelectedCell();
end

% --- Executes on button press in increaseZ.
function increaseZ_Callback(hObject, eventdata, handles)
% hObject    handle to increaseZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'windowListener',1);
newValue = getappdata(0, 'selectedZ')+1;
labelledImage = getappdata(0, 'labelledImageTemp_Resized');

if newValue <= size(labelledImage, 3)
    setappdata(0, 'selectedZ', newValue);
    set(handles.tbZFrame,'string',num2str(newValue));
    set(handles.slider1,'Value', newValue);
    showSelectedCell();
end

% --- Executes on button press in decreaseZ.
function decreaseZ_Callback(hObject, eventdata, handles)
% hObject    handle to decreaseZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'windowListener',1);
newValue = getappdata(0, 'selectedZ')-1;
if newValue > 0
    setappdata(0, 'selectedZ', newValue);
    set(handles.tbZFrame,'string',num2str(newValue));
    set(handles.slider1,'Value', newValue);
    showSelectedCell();
end


% --- Executes on button press in modifyOutside.
function modifyOutside_Callback(hObject, eventdata, handles)
% hObject    handle to modifyOutside (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modifyOutside
toggleValue = get(hObject,'Value') == 1;

setappdata(0, 'canModifyOutsideGland', toggleValue)


% --- Executes on button press in hideLumen.
function modifyInsideLumen_Callback(hObject, eventdata, handles)
% hObject    handle to hideLumen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hideLumen
toggleValue = get(hObject,'Value') == 1;
setappdata(0, 'canModifyInsideLumen', toggleValue)


% --- Executes on button press in modifyInsideLumen.
function hideLumen_Callback(hObject, eventdata, handles)
% hObject    handle to modifyInsideLumen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modifyInsideLumen
toggleValue = get(hObject,'Value') == 1;
setappdata(0, 'hideLumen', toggleValue)
showSelectedCell();


% --- Executes on button press in btRemove.
function btRemove_Callback(hObject, eventdata, handles)
% hObject    handle to btRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellId = getappdata(0, 'cellId');

    answer = questdlg(['Are you sure to remove cell ', num2str(getappdata(0, 'cellId')) , '?'], ...
        'Remove cell', ...
        'Yes','No', 'No');
    if strcmp(answer, 'Yes')
        if cellId > 0
            labelledImage = getappdata(0, 'labelledImageTemp');
            labelledImage(labelledImage == cellId) = 0;
            setappdata(0, 'labelledImageTemp', labelledImage);
        else
            lumenImage = getappdata(0, 'lumenImageTemp');
            lumenImage(lumenImage == 1) = 0;
            setappdata(0, 'lumenImageTemp', lumenImage);
        end
        updateResizedImage();
    end

showSelectedCell();

% --- Executes on button press in btAddCell.
function btAddCell_Callback(hObject, eventdata, handles)
% hObject    handle to btAddCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
labelledImage = getappdata(0, 'labelledImageTemp');
newValue = max(labelledImage(:)) + 1;
setappdata(0, 'cellId', newValue);
set(handles.tbCellId,'string',num2str(newValue));
insertROI_Callback(hObject, eventdata, handles)


% --- Executes on button press in btMergeCells.
function btMergeCells_Callback(hObject, eventdata, handles)
% hObject    handle to btMergeCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

labelledImage = getappdata(0, 'labelledImageTemp');

prompt = {'Enter cells to be merged (comma-separated): E.g. 20,25'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer) == 0
    cellsToMergeStr = strtrim(strsplit(answer{1}, ','));
    cellsToMerge = cellfun(@str2double, cellsToMergeStr);
    if length(cellsToMerge) > 1
        labelledImageTmp = mergeLabelsOfImage(labelledImage, cellsToMerge);
        setappdata(0, 'labelledImageTemp', labelledImageTmp);
        updateResizedImage();
        showSelectedCell();
    else
        errordlg('You should add more than 1 cell label', 'MEC!');
    end
end

% --- Executes on button press in chBoxShowAll.
function chBoxShowAll_Callback(hObject, eventdata, handles)
% hObject    handle to chBoxShowAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
toggleValue = get(hObject,'Value') == 1;
setappdata(0, 'showAllCells', toggleValue)
showSelectedCell();

function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if getappdata(0,'windowListener')==1
    %% If you are not modifying ROIs
    if strcmp(eventdata.Source.SelectionType, 'normal')
        try
            pos = round(eventdata.Source.CurrentObject.Parent.CurrentPoint);
            pos = pos(1,1:2);
            
            labelledImage = getappdata(0, 'labelledImageTemp');
            labelledImageZ = labelledImage(:,:,getappdata(0, 'selectedZ'))';
            selectedCell = labelledImageZ(pos(1), pos(2));

            setappdata(0,'cellId',selectedCell);
            set(handles.tbCellId,'string',num2str(selectedCell));

            showSelectedCell();
        catch
        end
    end
end


% --- Executes on button press in btRemove2DCell.
function btRemove2DCell_Callback(hObject, eventdata, handles)
% hObject    handle to btRemove2DCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellId = getappdata(0, 'cellId');
answer = questdlg(['Are you sure to remove cell ', num2str(getappdata(0, 'cellId')) , ' on this frame?'], ...
	'Remove cell', ...
	'Yes','No', 'No');
if strcmp(answer, 'Yes')
    if cellId > 0
        labelledImage = getappdata(0, 'labelledImageTemp');
        labelledImage_selectedZ = labelledImage(:, :, getappdata(0, 'selectedZ'));
        labelledImage_selectedZ(labelledImage_selectedZ == cellId) = 0;
        labelledImage(:, :, getappdata(0, 'selectedZ')) = labelledImage_selectedZ;
        setappdata(0, 'labelledImageTemp', labelledImage);
    else
        lumenImage = getappdata(0, 'lumenImageTemp');
        lumenImage_selectedZ = lumenImage(:, :, getappdata(0, 'selectedZ'));
        lumenImage_selectedZ(lumenImage_selectedZ == 1) = 0;
        lumenImage(:, :, getappdata(0, 'selectedZ')) = lumenImage_selectedZ;
        setappdata(0, 'lumenImageTemp', lumenImage);
    end
    updateResizedImage();
end
showSelectedCell();

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
numZ = round(get(hObject,'Value'));
setappdata(0, 'selectedZ', numZ);
set(handles.tbZFrame,'string',num2str(numZ));
showSelectedCell();

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
imageSequence = getappdata(0, 'imageSequence');

set(hObject,'Max',size(imageSequence,3));
set(hObject,'Value',1);
set(hObject,'Min',1);
set(hObject,'SliderStep',[1 1]./(size(imageSequence,3)-1));
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in showBackground.
function showBackground_Callback(hObject, eventdata, handles)
% hObject    handle to showBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
toggleValue = get(hObject,'Value') == 1;
setappdata(0, 'showBackground', toggleValue)
showSelectedCell();


% --- Executes on button press in btFillHoles.
function btFillHoles_Callback(hObject, eventdata, handles)
% hObject    handle to btFillHoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
labelledImage = getappdata(0,'labelledImageTemp');
lumenImage = getappdata(0,'lumenImageTemp');
selectedZ = getappdata(0, 'selectedZ');
labelledImageZ = uint16(labelledImage(:,:,selectedZ));
lumenImageZ = lumenImage(:,:,selectedZ);

invalidRegionZ = labelledImageZ == 0 & ~lumenImageZ;
invalidRegionZ = imfill(invalidRegionZ == 0, 'holes') & lumenImageZ == 0;
 
labelledImage(:,:,selectedZ) = fill0sWithCells(labelledImageZ, labelledImageZ, invalidRegionZ == 0);
setappdata(0,'labelledImageTemp',labelledImage);
updateResizedImage();
showSelectedCell();


% --- Executes on button press in btMaintainBigObject.
function btMaintainBigObject_Callback(hObject, eventdata, handles)
% hObject    handle to btMaintainBigObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
labelledImageResized = getappdata(0,'labelledImageTemp_Resized');
labelledImage = getappdata(0,'labelledImageTemp');
progressBar = waitbar(0, 'Removing splitted objects', 'WindowStyle', 'modal');

uniqIds = unique(labelledImageResized);
uniqIds = uniqIds(2:end);

for idCell = uniqIds'
    mask3dResized = labelledImageResized == idCell;
    vol = regionprops3(mask3dResized,'Volume');
    % Update waitbar and message
    waitbar(idCell/max(uniqIds))
    if size(vol,1)>1
        mask3d = labelledImage == idCell;
        vol = regionprops3(mask3d,'Volume');
        [~,idx] = max(vol.Volume);
        labelledImage(mask3d) = 0;
        labelledImage(bwlabeln(mask3d) == idx) = idCell;
    end
end
setappdata(0,'labelledImageTemp',labelledImage)
updateResizedImage();
close(progressBar)
showSelectedCell();