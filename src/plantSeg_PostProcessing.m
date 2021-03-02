function plantSeg_PostProcessing(outputDir, fileName)
%PIPELINE Summary of this function goes here
%   Detailed explanation goes here
    mkdir(fullfile(outputDir, 'Cells', 'OutputLimeSeg'));
    mkdir(fullfile(outputDir, 'ImageSequence'));
    mkdir(fullfile(outputDir, 'Lumen', 'SegmentedLumen'));
    mkdir(fullfile(outputDir, 'Results'));
    mkdir(fullfile(outputDir, 'ApicalAndBasal_Labelled'));


    if exist(fullfile(outputDir, 'Results', 'zScaleOfGland.mat'), 'file') == 0
        zScale = inputdlg('Insert z-scale of Gland');
        zScale = str2double(zScale{1});

        save(fullfile(outputDir, 'Results', 'zScaleOfGland.mat'), 'zScale');
    else
        load(fullfile(outputDir, 'Results', 'zScaleOfGland.mat')); 
    end
    
    if exist(fullfile(outputDir, 'Results', 'pixelScaleOfGland.mat'), 'file') == 0
        pixelScale = inputdlg('Insert pixel width of Gland');
        pixelScale = str2double(pixelScale{1});

        save(fullfile(outputDir, 'Results', 'pixelScaleOfGland.mat'), 'pixelScale');
    else
        load(fullfile(outputDir, 'Results', 'pixelScaleOfGland.mat')); 
    end
    
    resizeImg = 0.25;

    tipValue = 0;

    imageSequenceFiles = [dir(fullfile(outputDir, 'ImageSequence/*.tif'));dir(fullfile(outputDir, 'ImageSequence/*.tiff'))];
    NoValidFiles = startsWith({imageSequenceFiles.name},'._','IgnoreCase',true);
    imageSequenceFiles=imageSequenceFiles(~(NoValidFiles));
      
    %% if there is only 1 image, convert to imageSequence, then load.
    if size(imageSequenceFiles,1) == 1
        fname = fullfile(imageSequenceFiles.folder, imageSequenceFiles.name);
        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:num_images
            demoImg = imread(fname, k);
            imwrite(uint16(demoImg) , [imageSequenceFiles.folder '\image' num2str(k,'%03.f') '.tif']) ;
        end
        mkdir([imageSequenceFiles.folder,'\rawImageSequence\'])
        movefile(fname, [imageSequenceFiles.folder,'\rawImageSequence\' imageSequenceFiles.name]);
    end
    
    if exist(fullfile(outputDir, 'Results', '3d_layers_info.mat'), 'file')
        load(fullfile(outputDir, 'Results', '3d_layers_info.mat'))
    else
        colours = [];
        selpath = fullfile(outputDir, fileName);
        tiff_info = imfinfo(selpath); % return tiff structure, one element per image
        tiff_stack = imread(selpath, 1) ; % read in first image
        %concatenate each successive tiff to tiff_stack
        for ii = 2 : size(tiff_info, 1)
            temp_tiff = imread(selpath, ii);
            tiff_stack = cat(3 , tiff_stack, temp_tiff);
        end
        
        %set the background to the '0' label
        if min(tiff_stack(:))==1
            labelledImage = double(tiff_stack)-1;
        else
            labelledImage = double(tiff_stack);
        end
        
        
        if size(dir(fullfile(outputDir, 'Lumen/SegmentedLumen', '*.tif')),1) > 0
            [labelledImage, lumenImage] = processLumen(fullfile(outputDir, 'Lumen', filesep), labelledImage, resizeImg, tipValue);
        else
            [indx,~] = listdlg('PromptString',{'Lumen selection'},'SelectionMode','single','ListString',{'Select the biggest cell','Draw in matlab','Lumen from VolumeSegmenter'});
            switch indx                   
                case 1
                    %%Posible idea: try catch this line and if an error occurs get
                    %%the biggest 'cell' from plantSeg
                    %[labelledImage, lumenImage] = inferLumen(labelledImage);

                    [cellsVolume] = regionprops3(labelledImage, 'Volume');
                    [~, lumenIndex] = max(table2array(cellsVolume));
                    lumenImage = labelledImage == lumenIndex;
                    labelledImage(labelledImage == lumenIndex) = 0;
                case 2
                    lumenImage = zeros(size(labelledImage));
                case 3
                    load(fullfile(outputDir, 'lumen.mat'),'labels')
                    lumenImage = double(labels)==1;
            end
            
        end
        
        % Keep only the biggest object
%         objectsImage = bwlabeln(labelledImage>0);
%         objectsVolume = regionprops3(objectsImage, 'Volume');
%         [~, biggestObject] = max(table2array(objectsVolume));
%         labelledImage(objectsImage ~= biggestObject) = 0;
%         outsideGland = getOutsideGland(labelledImage);
          outsideGland = labelledImage == 0;
%         labelledImage = fill0sWithCells(labelledImage, labelledImage, outsideGland | lumenImage);
%         labelledImage(lumenImage) = 0;
        
        idLumen = unique(lumenImage(:,:,:));
        if ismember(idLumen,1)
            %% Put both lumen and labelled image at a 90 degrees
            orientationGland = regionprops3(lumenImage>0, 'Orientation');
            glandOrientation = -orientationGland.Orientation(1);
        else
            glandOrientation = 0;
        end
        
        [basalLayer, apicalLayer, colours] = postprocessGland(labelledImage,outsideGland, lumenImage, outputDir, colours, tipValue);
    end
    
    segmentationPostProcessing(labelledImage,lumenImage,apicalLayer,basalLayer,outputDir,resizeImg,tipValue,glandOrientation,colours);
    
end

