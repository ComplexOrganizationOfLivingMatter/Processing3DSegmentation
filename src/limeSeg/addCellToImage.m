function [labelledImage] = addCellToImage(pixelLocations, labelledImage, numCell)
%ADDCELLTOIMAGE Summary of this function goes here
%   We added for the same x,y several zs, because we found that some of the
%   zs were not completed (i.e. some zs of some cells were composed by only
%   a few pixels).

%     for numPixel = 1:size(pixelLocations, 1)
% %         zPixels = pixelLocations(numPixel, 3)-numDepth:1:pixelLocations(numPixel, 3)+numDepth;
% %         zPixels(zPixels < 1) = [];
% 
%         %We add these zeros to get the final size of labelled image
%         %But we do not add the perim of the cell because we only want the
%         %locations from the alphaShape
%         labelledImage(pixelLocations(numPixel, 1)+1, pixelLocations(numPixel, 2)+1, pixelLocations(numPixel, 3)+1) = 0;
%     end

    %Real final size of labelledImage
    labelledImage(max(pixelLocations(:, 1)), max(pixelLocations(:, 2)), max(pixelLocations(:, 3))) = 0;

    %filledCell = imfill(labelledImage == numCell, [pixelLocations(numPixel, 1)+1, pixelLocations(numPixel, 2)+1, pixelLocations(numPixel, 3)+1], 4);
    
    [labelledImage] = smoothObject(labelledImage,pixelLocations, numCell);
end

