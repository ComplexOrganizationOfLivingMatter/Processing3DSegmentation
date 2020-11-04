function limeSeg_PostProcessing(outputDir)
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

    tipValue = 5;
    
    imageSequenceFiles = dir(fullfile(outputDir, 'ImageSequence/*.tif'));
    NoValidFiles = startsWith({imageSequenceFiles.name},'._','IgnoreCase',true);
    imageSequenceFiles=imageSequenceFiles(~NoValidFiles);
    demoFile =  imageSequenceFiles(3);
    demoImg = imread(fullfile(demoFile.folder, demoFile.name));
    
    imageSequence = zeros(size(demoImg,1),size(demoImg,2),size(imageSequenceFiles, 1));
    for numImg = 1:size(imageSequenceFiles, 1)
        actualFile = imageSequenceFiles(numImg);
        actualImg = imread(fullfile(actualFile.folder, actualFile.name));
        imageSequence(:, :, numImg) = actualImg;
    end
    imgSize = size(imageSequence);
    imgSize(1:2)= imgSize(1:2).*resizeImg;

    if exist(fullfile(outputDir, 'Results', '3d_layers_info.mat'), 'file')
        load(fullfile(outputDir, 'Results', '3d_layers_info.mat'))
    else
        colours = [];
        
        [labelledImage] = processCells(fullfile(outputDir, 'Cells', 'OutputLimeSeg'), resizeImg, imgSize, zScale, tipValue);
        [indx,~] = listdlg('PromptString',{'Valid region selection'},'SelectionMode','single','ListString',{'valid Region','No'});
        switch indx
            case 1
                 load(fullfile(outputDir, 'validRegion.mat'))
                 validRegion = addTipsImg3D(tipValue,imresize3(double(labels)==1,imgSize,'nearest'));
                 labelledImage(validRegion==0)=0;
        end
        [indx,~] = listdlg('PromptString',{'Lumen selection'},'SelectionMode','single','ListString',{'Infer lumen','Draw in matlab', 'LimeSeg lumen', 'Photoshop lumen','Volume Segmenter lumen'});
        switch indx
            case 1
                [labelledImage, lumenImage] = inferLumen(labelledImage);
            case 2
                lumenImage = zeros(size(labelledImage));
            case 3
                lumenImage = processCells(fullfile(outputDir, 'Lumen', 'SegmentedLumen'), resizeImg, imgSize, zScale, tipValue)>0;
                labelledImage(lumenImage) = 0;
            case 4
                [labelledImage, lumenImage] = processLumen(fullfile(outputDir, 'Lumen', filesep), labelledImage, resizeImg, tipValue);
            case 5
                load(fullfile(outputDir, 'validRegion.mat'))
                lumenImage = logical(addTipsImg3D(tipValue,imresize3(double(labels)==2,imgSize,'nearest')));

        end
        
        try
            outsideGland = double(processCells(fullfile(outputDir, 'OutsideGland'), resizeImg, imgSize, zScale, tipValue)) == 0;
            labelledImage(outsideGland) = 0;
        catch ex
            %% Get invalid region
            [outsideGland] = getOutsideGland(labelledImage);
            %It add pixels and remove some
            validRegion = imfill(bwmorph3(labelledImage>0 | imdilate(lumenImage, strel('sphere', 5)), 'majority'), 'holes');
            %outsideGland = validRegion == 0;
            questionedRegion = imdilate(outsideGland, strel('sphere', 2));
            outsideGland(questionedRegion) = ~validRegion(questionedRegion);
        end
        
        %% Put both lumen and labelled image at a 90 degrees
        if sum(lumenImage(:))>0
            outsideGland(lumenImage) = 0;

            orientationGland = regionprops3(lumenImage>0, 'Orientation');
            glandOrientation = -orientationGland.Orientation(1);
            %labelledImage = imrotate(labelledImage, glandOrientation);
            %lumenImage = imrotate(lumenImage, glandOrientation);
        else
            glandOrientation=0;
        end      
        
        answer = questdlg('Would you fill empty space with cell labels?','Choose', 'yes', 'no', 'no');
        
        % Handle response
        if strcmp(answer, 'yes')
            labelledImage = fill0sWithCells(labelledImage, labelledImage, outsideGland | lumenImage);
        end
        
        [basalLayer, apicalLayer, colours] = postprocessGland(labelledImage, outsideGland, lumenImage, outputDir, colours, tipValue);
        
        labelledImage = addTipsImg3D(-tipValue,labelledImage);
        lumenImage = addTipsImg3D(-tipValue,lumenImage);
    end
    imgSize(1:2)= imgSize(1:2)./resizeImg;
    labelledImage = imresize3(labelledImage,imgSize,'nearest');
    lumenImage = imresize3(double(lumenImage),imgSize,'nearest');
        
    segmentationPostProcessing(labelledImage,lumenImage,apicalLayer,basalLayer,outputDir,resizeImg,tipValue,glandOrientation,colours)
end