function PTMResult = PTM01(pos, res, type)

% Colin Ophus - 2021 Jan
% modified by Zezhou Li - 2022 Oct
% 71 - This is the first of a series of scripts designed to index Philipp's
% decahedral atomic coordinates. This script constructs the neighbor
% network, in order to make NN distance calculations

arguments
    pos(3,:) double
    res(1,1) double  = 0.3434 % Angstrom
    type(1,:) double = []
end

num = length(pos);
% make all positions positive
pos = pos - mean(pos,2);
szMax = ceil(max(pos(:)))+round(5/res); % add edge of 5 Angstrom
pos = pos + szMax;

% Output original coordinates
PTMResult.xyz = pos';
PTMResult.res = res;
if nargin > 2
    if length(type) == num
        PTMResult.type = type;
    else
        error('Please input right number of atom types and atom positions !\n')
    end
end

[~, radiusNNmean, radiusNNmax] = CalculateLatticeConstant(pos,res);
PTMResult.radiusNNmean = radiusNNmean;
PTMResult.radiusNNmax = radiusNNmax;

% calculate the number of Nearest Neighbors with their distance and indice
distanceMatrix = CalculateDistanceList(PTMResult.xyz'*res);
distanceMatrix(distanceMatrix > radiusNNmax) = 0;
PTMResult.NNnum = sum(distanceMatrix>0,2);
PTMResult.NNlist = cell(num,1);
for i = 1:num
    PTMResult.NNlist{i} = find(distanceMatrix(i,:));
end

% plot # of NNs
figure(55)
clf
histogram(PTMResult.NNnum,0:max(PTMResult.NNnum))
xlabel('Number of Neighbors')
ylabel('Count')

end

%% child function
function dist = CalculateDistanceList(pos)
arguments
    pos(3,:) double
end
% tic
totalAtomNum = size(pos,2);
distMatrix = zeros(totalAtomNum,totalAtomNum,3);
for i = 1:3
    distMatrix(:,:,i) = pos(i,:)-pos(i,:)';
end
dist = vecnorm(distMatrix,2,3);

end

