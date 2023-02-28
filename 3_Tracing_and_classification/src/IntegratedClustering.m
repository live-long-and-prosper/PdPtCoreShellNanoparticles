function [idx,varargout] = IntegratedClustering(X,method,clusterNum,plotResult,excludeOutlierDuringClustering,colormap)
%% written by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn
arguments
    X(:,:) double % size:n*m，n is the length of data, m is the dimension size of each datum
    method(1,:) char {mustBeMember(method,{'kmeans','spectral','gaussian','dbscan'})} = 'kmeans'
    clusterNum(1,1) double {mustBePositive, mustBeInteger} = 2
    plotResult(1,1) double {mustBeInteger} = 1 % 0: do not plot result, 1: plot result in new figure, 2: plot in figure handle outside
    excludeOutlierDuringClustering(1,1) double = 0 % use isoutlier to remove outliers in the data
    colormap(1,:) char = 'parula'
end

[num,dimensions] = size(X);
fprintf(['Clustering %d points (%d dimensions) using ',method,' method'],num,dimensions)
% normalize
for i = 1:dimensions
    X(:,i) = rescale(X(:,i));
end

if excludeOutlierDuringClustering
    outliers = isoutlier(X,1);
    outliers = logical(sum(outliers,2));
    X_cluster = X(~outliers,:);
    X_outliers = X(outliers,:);
    OutlierNum = sum(outliers);
else
    X_cluster = X;
    X_outliers = [];
    OutlierNum = 0;
end

if strcmp(method,'kmeans')
    idx = kmeans(X_cluster,clusterNum);
elseif strcmp(method,'spectral')
    idx = spectralcluster(X_cluster,clusterNum);
elseif strcmp(method,'gaussian')
    % use kmeans first to get reasonable result
    idx = kmeans(X_cluster,clusterNum);
    gm = fitgmdist(X_cluster,clusterNum,'Options',statset('MaxIter',1000,'TolFun',1e-5),'Start',idx);
    % although MATLAB use kmeans++ to get initial points, but the result is
    % worse than use kmeans function seperately at first!
    idx = cluster(gm,X_cluster);
elseif strcmp(method,'dbscan') % dbscan needs to manually adjust two parameters...
    epsilon = 0.17;
    minpts = round(num/10);
    idx = dbscan(X_cluster,epsilon,minpts);
end

clusterTypes = unique(idx);
clusterNum_real = length(clusterTypes);
fprintf(' to %d types\n',clusterNum_real)

% cluster the outliers afterwards
clusterCenters = zeros(clusterNum_real,dimensions);
for i = 1:clusterNum_real
    clusterCenters(i,:) = mean(X_cluster(idx==clusterTypes(i),:),1);
end
if OutlierNum > 0
    distMatrix = zeros(OutlierNum,clusterNum_real,dimensions);
    for i = 1:dimensions
        distMatrix(:,:,i) = X_outliers(:,i)-clusterCenters(:,i)';
    end
    tempDistance = vecnorm(distMatrix,2,3);
    [~,indice] = min(tempDistance,[],2);
    idx_all = zeros(num,1);
    idx_all(~outliers) = idx;
    idx_all(outliers) = clusterTypes(indice);
    idx = idx_all;
end

% sort clusters (except dbscan)
if ~strcmp(method,'dbscan')
    % calculate center of mass
    tempValue = zeros(1,clusterNum_real);
    for i = 1:clusterNum_real
        tempValue(i) = vecnorm(mean(X(idx==clusterTypes(i),:),1));
    end
    % sort
    [~,sortIndice] = sort(tempValue);
    for i = 1:clusterNum_real
        idx(idx==sortIndice(i)) = i+clusterNum_real;
    end
    idx = idx-clusterNum_real;
    clusterTypes = 1:clusterNum_real;
end

for i = 1:clusterNum_real
    fprintf('cluster No.%d : %d points (%.1f%%)\n',clusterTypes(i), sum(idx == clusterTypes(i)), sum(idx == clusterTypes(i))/num*100)
end

% plot
if plotResult > 0
    pointSize = 15;
    BinWidth = 0.02;
    if excludeOutlierDuringClustering
        BinWidth = BinWidth/10;
    end
    % figure handle
    if plotResult == 1
        figure
    else
        gcf;
    end
    if dimensions == 3
        cmaps = feval(colormap,clusterNum_real+1);
        for i = 1:clusterNum_real
            plotAtoms(X(idx == clusterTypes(i),:)',clusterTypes(i),{pointSize,cmaps(i,:)})
        end
    elseif dimensions == 2
        cmaps = feval(colormap,clusterNum_real+1);
        gscatter(X(:,1),X(:,2),idx,cmaps(1:end-1,:),'.',pointSize)
    elseif dimensions == 1
        hold on
        cmaps = feval(colormap,clusterNum_real+1);
        for i = 1:clusterNum_real
            histogram(X(idx == clusterTypes(i)),'FaceColor',cmaps(i,:),'BinWidth',BinWidth)
        end
        legend(cellstr(string(clusterTypes)));
    elseif dimensions > 3
        fprintf('plot with first 3 dimensions\n')
        cmaps = feval(colormap,clusterNum_real+1);
        subplot(121)
        for i = 1:clusterNum_real
            plotAtoms(X(idx == clusterTypes(i),1:3)',clusterTypes(i),{pointSize,cmaps(i,:)})
        end
        subplot(122)
        for i = 1:clusterNum_real
            plotAtoms(X(idx == clusterTypes(i),[1 2 4])',clusterTypes(i),{pointSize,cmaps(i,:)})
        end
    end
    title(['clustering by ',method],'FontSize',16)
end

if nargout > 1
    varargout{1} = clusterNum_real;
end

end
