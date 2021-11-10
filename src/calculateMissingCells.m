function [answer, apical3dInfo, notFoundCellsApical, basal3dInfo, notFoundCellsBasal] = calculateMissingCells(labelledImage, lumenImage, apicalLayer, basalLayer, colours, noValidCells)
%CALCULATEMISSINGCELLS Calculate and plot missing cells in both layers
%   We calculate which cells are missing from apical and basal layer, which
%   may be a mistake. In addition, we plot the layers and where are the
%   missing cells.

    allCells = 1:max(labelledImage(:));
    allCells(allCells == 0) = [];
    
    %% Cells that does not exist
    cellsPixels = regionprops3(labelledImage, 'Volume');
    cellsWithoutPixels = find(cellsPixels.Volume == 0);
    
    [apical3dInfo] = calculateNeighbours3D(apicalLayer, 2, apicalLayer == 0);
    apical3dInfo = apical3dInfo.neighbourhood';
    
    if length(allCells) ~= length(apical3dInfo)
        addingCells = length(allCells) - length(apical3dInfo);
        apical3dInfo(end+addingCells) = {[]};
    end
    notFoundCellsApical = find(cellfun(@(x) isempty(x), apical3dInfo))';
    % Remove non-existing cells
    notFoundCellsApical = setdiff(notFoundCellsApical , cellsWithoutPixels);

    %Display missing cells in basal
    [basal3dInfo] = calculateNeighbours3D(basalLayer, 2, basalLayer == 0);
    basal3dInfo = basal3dInfo.neighbourhood';

    if length(allCells) ~= length(basal3dInfo)
        addingCells = length(allCells) - length(basal3dInfo);
        basal3dInfo(end+addingCells) = {[]};
    end
    notFoundCellsBasal = find(cellfun(@(x) isempty(x), basal3dInfo))';
    notFoundCellsBasal = setdiff(notFoundCellsBasal , cellsWithoutPixels);
    
    answer = 'Yes';

    
end

