function [goodPos,tooClosePos] = PolynomialFit3DAtoms(vol,roughPos,Dmin,res,fitRadius)
%% modified from Jianwei Miao's previous AET work by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn
arguments
    vol(:,:,:) double
    roughPos(3,:) double % the list should been sort by integrated value from low to high
    Dmin(1,1) double = 1.5 % Angstrom
    res(1,1) double = 0.3434 % Angstrom
    fitRadius(1,1) double = 3 % voxel
end

RoughAtomPosListForTrace = round(roughPos');
num = size(roughPos,2);
MaxIter = 14;
% CritIter = 7;
Q = 0.5;
Alpha = 1;
cropHalfSize = fitRadius;
% get polynomial power array
fitCoeff = [];
for i=0:4
    for j=0:4
        for k=0:4
            if i+j+k <= 4                
                if max([i j k]) == 4
                    fitCoeff(end+1,:) = [i j k -1];
                else                
                    fitCoeff(end+1,:) = [i j k 0];
                end
            end
        end
    end
end

[X,Y,Z]      = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
SphereInd    = find(X.^2+Y.^2+Z.^2 <=(fitRadius+0.5)^2);
XYZdata.X    = X(SphereInd); XYZdata.Y = Y(SphereInd); XYZdata.Z = Z(SphereInd);

Orders       = fitCoeff(:,1:3);
PosArr       = zeros(num,3);
TotPosArr    = zeros(num,3);

exitFlagArr  = zeros(1,num);
CoeffArr     = repmat(fitCoeff(:,4),[1,num]);

% perform the main tracing loop
parfor i = 1:num
    endFlag = 0;
    iterNum = 0;
    tempRoughPos = RoughAtomPosListForTrace(i,:);
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4; % end the cycle
          endFlag = 1;
        end
        cropXind = tempRoughPos(1) + (-cropHalfSize:cropHalfSize);
        cropYind = tempRoughPos(2) + (-cropHalfSize:cropHalfSize);
        cropZind = tempRoughPos(3) + (-cropHalfSize:cropHalfSize);

        cropVol = vol(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2 ) / cropHalfSize^2 );
        
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        opts = optimset('Display','off');
        
        [p1,~] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX == -100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1; % fail
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            TotPosArr(i,:) = PosArr(i,:) + tempRoughPos;
        end
    end
end

tempIndice = find(TotPosArr(:,1)==0);
if ~isempty(tempIndice)
    TotPosArr(tempIndice,:) = [];
    fprintf('%d atoms fit failed !\n', length(tempIndice))
end

% iteratively delete too closed atoms until the atom number is not changed
flag  = 1;
goodPos = TotPosArr';
tooClosePos = zeros(size(goodPos));
tooCloseNum = 0;
while flag
    [minDistanceList,minDistanceList_index] = CalculateMinimumDistanceList(goodPos*res);
    list1 = find(minDistanceList < Dmin);
    if isempty(list1)
        flag = 0;
    else
        list2 = minDistanceList_index(list1);
        list3 = unique(min([list1;list2])); % remove too closed atoms
        tooClosePos(:,tooCloseNum+1:tooCloseNum+length(list3)) = goodPos(:,list3);
        tooCloseNum = tooCloseNum + length(list3);
        goodPos(:,list3) = [];
    end
end
if tooCloseNum < size(tooClosePos,2)
    tooClosePos(:,tooCloseNum+1:end) = [];
end

fprintf('%d atoms fitted, %d atoms deleted because there are too close, totally %d atoms\n',size(goodPos,2),size(tooClosePos,2),num-length(tempIndice))

end

%% child function
function [minDistanceList, varargout] = CalculateMinimumDistanceList(pos)
arguments
    pos(3,:) single
end

totalAtomNum = size(pos,2);
try
    distMatrix = single(zeros(totalAtomNum,totalAtomNum,3));
    for i = 1:3
        distMatrix(:,:,i) = pos(i,:)-pos(i,:)';
    end
    distanceMatrix = single(vecnorm(distMatrix,2,3));
    clear distMatrix
    distanceMatrix(distanceMatrix == 0) = nan;

    if nargout == 1
        minDistanceList = double(min(distanceMatrix));
    elseif nargout >= 2
        [minDistanceList,minDistanceList_index] = min(distanceMatrix);
        minDistanceList = double(minDistanceList);
        varargout{1} = double(minDistanceList_index);
        if nargout == 3
            DistaneMatrixSort = sort(distanceMatrix);
            if size(pos,2) >= 12
                varargout{2} = double(mean(DistaneMatrixSort(1:12,:)));
            else
                varargout{2} = double(mean(DistaneMatrixSort(1:end-1,:)));
                fprintf('number of input positions is less than 12, output the mean of all distances instead\n')
            end
        end
    end
catch % when memory is not enough...
    clear distMatrix
    fprintf('memory not enough, try doing one by one...')
    minDistanceList = single(zeros(1,totalAtomNum));
    if nargout == 1
        parfor i = 1:totalAtomNum
            minDistanceList(i) = min(vecnorm(pos(:,i) - pos));
        end
    elseif nargout == 2
        minDistanceList_index = single(zeros(1,totalAtomNum));
        parfor i = 1:totalAtomNum
            [minDistanceList(i), minDistanceList_index(i)] = min(vecnorm(pos(:,i) - pos));
        end
        varargout{1} = double(minDistanceList_index);
    elseif nargout == 3
        minDistanceList_index = single(zeros(1,totalAtomNum));
        minDistance12_List = single(zeros(1,totalAtomNum));
        parfor i = 1:totalAtomNum
            [tempList, tempIndice] = sort(vecnorm(pos(:,i) - pos));
            minDistanceList(i) = tempList(1);
            minDistanceList_index(i) = tempIndice(1);
            minDistance12_List(i) = mean(tempList(1:12));
        end
        varargout{1} = double(minDistanceList_index);
        varargout{2} = double(minDistance12_List);
    end
end
   
end
