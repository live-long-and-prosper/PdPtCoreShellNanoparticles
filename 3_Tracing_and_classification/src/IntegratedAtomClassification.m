function [type,pos_sort] = IntegratedAtomClassification(pos,res,vol,typeNum,dataTypeToClassify,classifyMethod,plotResult)
%% written by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn
arguments
    pos(3,:) single % use voxel as unit，must be positive (directly leading to positions in the volume)
    res(1,1) single % Angstrom，voxel size
    vol(:,:,:) single % 3D volume
    typeNum(1,1) double {mustBeInteger} = 2 % expected number of atom types
    dataTypeToClassify(1,:) char {mustBeMember(dataTypeToClassify,{'boxIntensity','totalIntensity','peakIntensity','meanIntensity_otsu'})} = 'boxIntensity'
    classifyMethod(1,:) char {mustBeMember(classifyMethod,{'kmeans','spectral','gaussian','totalDifference'})} = 'totalDifference'
    % 'totalDifference' refer to unbaised classification method (Yongsoo Yang et al. Nature 2017)
    plotResult(1,1) double {mustBeInteger} = 1 % 0: do not plot result, 1: plot result in new figure, 2: plot in figure handle outside
end

if sum(pos(:) < 0) > 0
    error('Please input positive positions that directly indicated atom positions in the 3Dvolume !\n')
end

if strcmp(classifyMethod,'totalDifference')
    dataTypeToClassify = 'boxIntensity';
end

if strcmp(dataTypeToClassify,'boxIntensity')
    % estimate reasonable box radius
    [~,meanBondLength] = CalculateLatticeConstant(pos,res);
    boxRadius = floor((round(meanBondLength/res)-1)/2);
    interpMethod = 'cubic';
    if strcmp(classifyMethod,'totalDifference')
        [boxINT_total,boxINT_list] = GetIntegratedIntensityInBox(vol,pos,boxRadius*2+1,interpMethod);
    else
        data = GetIntegratedIntensityInBox(vol,pos,boxRadius*2+1,interpMethod);
    end
    pos_sort = pos;
else
    mask3DPadSize = 0;
    [~, statsLists] = WaterShedSegment3D(vol, pos, 'standard', mask3DPadSize);
    pos_watershed = single(statsLists.otsuC_interp);
    [pos_sort,~,~] = PairTwoModels(pos,pos_watershed);
end
num = size(pos_sort,2);

if strcmp(dataTypeToClassify,'totalIntensity')
    data = statsLists.totalINT;
elseif strcmp(dataTypeToClassify,'peakIntensity')
    data = statsLists.maxINT;
elseif strcmp(dataTypeToClassify,'meanIntensity_otsu')
    data = statsLists.meanINT_interp;
end

if strcmp(classifyMethod,'totalDifference')
    fprintf('Execute classification using method in 2017 FePt Nature.\n')
    boxINT_total_sort = sort(boxINT_total);
    thresholds = linspace(min(boxINT_total_sort),max(boxINT_total_sort),typeNum+1);
    type = ones(1,num);
    for i = 1:typeNum-1
        type(boxINT_total > thresholds(i+1)) = i+1;
    end
    boxLength = size(boxINT_list,1);
    averageBox = zeros(boxLength,typeNum);
    for i = 1:typeNum
        averageBox(:,i) = mean(boxINT_list(:,type==i),2);
    end
    flag = 1;
    while flag
        type_matrix = zeros(typeNum,num);
        for i = 1:typeNum
            type_matrix(i,:) = vecnorm(boxINT_list-averageBox(:,i),1);
        end
        [~,type_new] = min(type_matrix);
        if vecnorm(type-type_new,1) == 0 % if classification result not changed
            flag = 0;
        else
            type = type_new;
            for i = 1:typeNum
                averageBox(:,i) = mean(boxINT_list(:,type==i),2);
            end
        end
    end
    if plotResult
        colormap = 'parula';
        BinWidth = (max(boxINT_total)-min(boxINT_total))*0.02;
        figure
        hold on
        cmaps = feval(colormap,typeNum+1);
        for i = 1:typeNum
            histogram(boxINT_total(type == i),'FaceColor',cmaps(i,:),'BinWidth',BinWidth)
        end
        legend(cellstr(string(1:typeNum)));
    end
    for i = 1:typeNum
        fprintf('cluster No.%d : %d points (%.1f%%)\n', i, sum(type==i), sum(type==i)/num*100)
    end
else
    fprintf(['Execute 1D clustering using', classifyMethod, ' method.\n'])
    excludeOutlierDuringClustering = 0;
    % pay attention to transposition
    type = IntegratedClustering(transpose(data),classifyMethod,typeNum,plotResult,excludeOutlierDuringClustering);
end

end
