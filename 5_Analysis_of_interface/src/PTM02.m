function PTMResult = PTM02(PTMResult, AdaptiveFitMode, DrawFlag, numIterICP, numTheta, scoreThreshold)
arguments
    PTMResult struct
    AdaptiveFitMode(1,1) double {mustBeNonnegative} = 1 % Variable standard bond length
    DrawFlag(1,1) double {mustBeNonnegative} = 0 % plot or not
    numIterICP(1,1) double {mustBePositive, mustBeInteger} = 20
    numTheta(1,1) double {mustBePositive, mustBeInteger} = 5
    scoreThreshold(1,1) double {mustBePositive} = 8
end
%% Polyhedral Template Matching Main Code
% Colin Ophus - 2020 Dec
% modified by Zezhou Li - 2022 Oct

% 72 - measure FCC and HCP (twin) ordering from all coordinates using a
% cost function based on polyhedral matching.

% numIterICP = 20;

numTesselations = 4;
% numTheta = 5;
% step is 24° to search 72° (for decahedron and isocahedron) and 120° (for fcc and hcp)
phiArray = linspace(0,pi*2/3,numTheta+1);
phiArray(end) = []; 

radiusNN = PTMResult.radiusNNmean / PTMResult.res; % mean bond length
scoreRadiusRange = [0 0.5] * radiusNN;  % for cost function

% Define the uvw array for primary orientation
% uAll = [ ...
%     [0 0 1];
%     [0 1 1]/sqrt(2);
%     [1 1 1]/sqrt(3)];
uvwOrient = [ ...
    [0 0 1];
    [0 1 0];
    [1 0 0]];
f = [1 2 3];
% tesselate
for a0 = 1:numTesselations
    [uvwOrient,f] = tesselate(uvwOrient,f);
end
% symmetry
uvwOrient = [ ...
    uvwOrient;
    uvwOrient .* [-1 1 1];
    uvwOrient .* [-1 -1 1];
    uvwOrient .* [1 -1 1]];
uvwOrient(:) = uvwOrient ./ sqrt(sum(uvwOrient.^2,2));
% remove duplicates
Nv = size(uvwOrient,1);
del = false(Nv,1);
r2max = 1e-4;
for a0 = 1:Nv
    if del(a0) == false
        d2a = (uvwOrient(:,1) - uvwOrient(a0,1)).^2 + ...
            (uvwOrient(:,2) - uvwOrient(a0,2)).^2 + ...
            (uvwOrient(:,3) - uvwOrient(a0,3)).^2;
        d2b = (uvwOrient(:,1) + uvwOrient(a0,1)).^2 + ...
            (uvwOrient(:,2) + uvwOrient(a0,2)).^2 + ...
            (uvwOrient(:,3) + uvwOrient(a0,3)).^2;
        inds = sort(find((d2a < r2max) | (d2b < r2max)));
        
        if length(inds) > 1
            del(inds(2:end)) = true;
        end
    end
end
uvwOrient(del,:) = [];
uvwOrient(:) = uvwOrient ./ sqrt(sum(uvwOrient.^2,2));

% Full sphere
uvwOrient = [uvwOrient; -uvwOrient];

% output
PTMResult.uvwOrient = sortrows(uvwOrient,[1 2 3]);


figure(67)
clf
scatter3(PTMResult.uvwOrient(:,2),...
    PTMResult.uvwOrient(:,1),...
    PTMResult.uvwOrient(:,3),'rs','filled')
axis equal
box on
view([0 0 1])

%% generate perfect 12-CN model for hcp、fcc、decahedron and isocahedron
% Math init
a_FCC = radiusNN * sqrt(2);
h_ClosePacking = a_FCC / sqrt(3);
a_Isocahedron = radiusNN * 4 / sqrt(10+2*sqrt(5));
h_Isocahedron = a_Isocahedron * sqrt(1-1/(4*sin(pi/5)^2));
r1 = radiusNN;
r2 = radiusNN / sqrt(3);
r3 = radiusNN * sqrt(3) / 2;
r4 = a_Isocahedron / (2*sin(pi/5));
t1 = linspace(0,2*pi,6+1)';
t1(end) = [];
t2 = linspace(0,2*pi,3+1)' + pi / 6;
t2(end) = [];
t3 = t2 + pi / 3;
t4 = linspace(0,2*pi,5+1)' + pi / 10;
t4(end) = [];
t5 = t4 + pi / 5;

% HCP basis
xyzHCP = [ ...
    [[cos(t1) sin(t1)]*r1 zeros(6,1)];      % 6 points in the plane
    [[cos(t2) sin(t2)]*r2 ones(3,1) * -h_ClosePacking];  % 3 points above
    [[cos(t2) sin(t2)]*r2 ones(3,1) *  h_ClosePacking];  % 3 points below
    ];
% FCC basis
xyzFCC = [ ...
    [[cos(t1) sin(t1)]*r1 zeros(6,1)];      % 6 points in the plane
    [[cos(t2) sin(t2)]*r2 ones(3,1) * -h_ClosePacking];  % 3 points above
    [[cos(t3) sin(t3)]*r2 ones(3,1) *  h_ClosePacking];  % 3 points below
    ];
% Decahedron basis
xyzDH = [ ...
    [[cos(t4) sin(t4)]*r3 ones(5,1)* radiusNN/2];  % 5 points in the plane above
    [[cos(t4) sin(t4)]*r3 ones(5,1)*-radiusNN/2];  % 5 points in the plane below
    [zeros(2,2)               [1;-1] * radiusNN];  % 2 points on the 5-fold axis
    ];
% Isocahedron basis
xyzIH = [ ...
    [[cos(t4) sin(t4)]*r4 ones(5,1)* h_Isocahedron];  % 5 points in the plane above
    [[cos(t5) sin(t5)]*r4 ones(5,1)*-h_Isocahedron];  % 5 points in the plane below
    [zeros(2,2)               [1;-1]*a_Isocahedron];  % 2 points on the 5-fold axis
    ];

figure;
subplot(141)
plotAtoms(xyzHCP',ones(1,12),struct('Atomsize',200,'DrawAlphaShape',1))
title('hcp model')
subplot(142)
plotAtoms(xyzFCC',ones(1,12),struct('Atomsize',200,'DrawAlphaShape',1))
title('fcc model')
subplot(143)
plotAtoms(xyzDH',ones(1,12),struct('Atomsize',200,'DrawAlphaShape',1))
title('decahedron model')
subplot(144)
plotAtoms(xyzIH',ones(1,12),struct('Atomsize',200,'DrawAlphaShape',1))
title('isocahedron model')

% output
PTMResult.numTesselations = numTesselations;
PTMResult.phiArray = phiArray;
PTMResult.xyzHCP = xyzHCP;
PTMResult.xyzFCC = xyzFCC;
PTMResult.xyzDH = xyzDH;
PTMResult.xyzIH = xyzIH;

% pre-compute all orientation arrays
N = [size(xyzHCP,1) size(xyzHCP,2) ...
    size(PTMResult.uvwOrient,1) length(phiArray)];
xyzHCPall = zeros(N);
xyzFCCall = zeros(N);
xyzDHall = zeros(N);
xyzIHall = zeros(N);
for a0 = 1:size(PTMResult.uvwOrient,1)
    vec3 = PTMResult.uvwOrient(a0,:);
    vec1 = cross([1e-8 0 1],vec3);
    vec2 = cross(vec3,vec1);
    vec1(:) = vec1 / norm(vec1);
    vec2(:) = vec2 / norm(vec2);
    vec3(:) = vec3 / norm(vec3);
    
    for a1 = 1:length(phiArray)
        phi = phiArray(a1);
        m = [cos(phi) sin(phi) 0;
            -sin(phi) cos(phi) 0;
            0 0 1];
        
        % rotate basis functions
        xyzHCPall(:,:,a0,a1) = (xyzHCP * m) * [vec1; vec2; vec3];
        xyzFCCall(:,:,a0,a1) = (xyzFCC * m) * [vec1; vec2; vec3];
        xyzDHall(:,:,a0,a1) = (xyzDH * m) * [vec1; vec2; vec3];
        xyzIHall(:,:,a0,a1) = (xyzIH * m) * [vec1; vec2; vec3];
    end
end

%% main loop
num = length(PTMResult.xyz);
fprintf('wait for about %d minutes...\n', ceil(num/35/60))
xyz = PTMResult.xyz;
NNlist = PTMResult.NNlist;
tempDataAll = zeros(4,num,4);
% parallel to get faster
parfor a2 = 1:num
    dxyz = xyz(NNlist{a2},:) - xyz(a2,:);
    localRadiusNN = mean(vecnorm(dxyz'));
    if AdaptiveFitMode == 0
        RadiusNN_ScaleRatio = 1;
    else
        RadiusNN_ScaleRatio = localRadiusNN/radiusNN;
    end
    % global search
    score_HCP = zeros(N(3:4));
    score_FCC = zeros(N(3:4));
    score_DH = zeros(N(3:4));
    score_IH = zeros(N(3:4));
    for a0 = 1:N(3)
        for a1 = 1:N(4)
            score_HCP(a0,a1) = calcScore(xyzHCPall(:,:,a0,a1)*RadiusNN_ScaleRatio,dxyz,scoreRadiusRange);
            score_FCC(a0,a1) = calcScore(xyzFCCall(:,:,a0,a1)*RadiusNN_ScaleRatio,dxyz,scoreRadiusRange);
            score_DH(a0,a1) = calcScore(xyzDHall(:,:,a0,a1)*RadiusNN_ScaleRatio,dxyz,scoreRadiusRange);
            score_IH(a0,a1) = calcScore(xyzIHall(:,:,a0,a1)*RadiusNN_ScaleRatio,dxyz,scoreRadiusRange);
        end
    end
    % 4 dimensions: local structure, a0 for best initial similarity score, a1 for best initial score, best score after icp fit 
    tempData = zeros(4,4);
    % hcp
    [tempData(1,1), Ind] = max(score_HCP,[],'all');
    [tempData(1,2), tempData(1,3)] = ind2sub(N(3:4),Ind);
    [tempData(2,1), Ind] = max(score_FCC,[],'all');
    [tempData(2,2), tempData(2,3)] = ind2sub(N(3:4),Ind);
    [tempData(3,1), Ind] = max(score_DH,[],'all');
    [tempData(3,2), tempData(3,3)] = ind2sub(N(3:4),Ind);
    [tempData(4,1), Ind] = max(score_IH,[],'all');
    [tempData(4,2), tempData(4,3)] = ind2sub(N(3:4),Ind);

    % icp (iterative cloud point algorithm)
    if size(dxyz,1)>=4
    % ICP for HCP
        xyzAlign = xyzHCPall(:,:,tempData(1,2),tempData(1,3))*RadiusNN_ScaleRatio;
        [TR, TT] = icp(dxyz', xyzAlign', numIterICP, 'Matching','kDtree');
        xyzAlign = (TR * xyzAlign' + TT)';
        score = calcScore(xyzAlign,dxyz,scoreRadiusRange);
        tempData(1,4) = score;
    % ICP for FCC
        xyzAlign = xyzFCCall(:,:,tempData(2,2), tempData(2,3))*RadiusNN_ScaleRatio;
        [TR, TT] = icp(dxyz', xyzAlign', numIterICP, 'Matching','kDtree');
        xyzAlign = (TR * xyzAlign' + TT)';
        score = calcScore(xyzAlign,dxyz,scoreRadiusRange);
        tempData(2,4) = score;
    % ICP for DH
        xyzAlign = xyzDHall(:,:,tempData(3,2), tempData(3,3))*RadiusNN_ScaleRatio;
        [TR, TT] = icp(dxyz', xyzAlign', numIterICP, 'Matching','kDtree');
        xyzAlign = (TR * xyzAlign' + TT)';
        score = calcScore(xyzAlign,dxyz,scoreRadiusRange);
        tempData(3,4) = score;
    % ICP for IH
        xyzAlign = xyzIHall(:,:,tempData(4,2), tempData(4,3))*RadiusNN_ScaleRatio;
        [TR, TT] = icp(dxyz', xyzAlign', numIterICP, 'Matching','kDtree');
        xyzAlign = (TR * xyzAlign' + TT)';
        score = calcScore(xyzAlign,dxyz,scoreRadiusRange);
        tempData(4,4) = score;
    end
    tempDataAll(:,a2,:) = tempData;
end

PTMResult.hcpData = zeros(num,7);
PTMResult.fccData = zeros(num,7);
PTMResult.dhData = zeros(num,7);
PTMResult.ihData = zeros(num,7);
PTMResult.hcpData(:,[1:3,7]) = squeeze(tempDataAll(1,:,:));
PTMResult.fccData(:,[1:3,7]) = squeeze(tempDataAll(2,:,:));
PTMResult.dhData(:,[1:3,7]) = squeeze(tempDataAll(3,:,:));
PTMResult.ihData(:,[1:3,7]) = squeeze(tempDataAll(4,:,:));

% output the final uvw orientations
sub = PTMResult.hcpData(:,1) > 0;
PTMResult.hcpData(sub,4:6) = PTMResult.uvwOrient(PTMResult.hcpData(sub,2),:);
sub = PTMResult.fccData(:,1) > 0;
PTMResult.fccData(sub,4:6) = PTMResult.uvwOrient(PTMResult.fccData(sub,2),:);
sub = PTMResult.dhData(:,1) > 0;
PTMResult.dhData(sub,4:6) = PTMResult.uvwOrient(PTMResult.dhData(sub,2),:);
sub = PTMResult.ihData(:,1) > 0;
PTMResult.ihData(sub,4:6) = PTMResult.uvwOrient(PTMResult.ihData(sub,2),:);

% calculate the similarity score of each structure
PTMResult.FCCscore = max(PTMResult.fccData(:,1),PTMResult.fccData(:,7));
PTMResult.HCPscore = max(PTMResult.hcpData(:,1),PTMResult.hcpData(:,7));
PTMResult.DHscore = max(PTMResult.dhData(:,1),PTMResult.dhData(:,7));
PTMResult.IHscore = max(PTMResult.ihData(:,1),PTMResult.ihData(:,7));

% compare the score to determine the most likely structure of each atom
[~,PTMResult.classify] = max([PTMResult.FCCscore,PTMResult.HCPscore,PTMResult.DHscore,PTMResult.IHscore],[],2);

% distinguish atoms with different local structures by threshold
% thres = 6; % score threshold is set as 6
PTMResult.scoreThreshold = scoreThreshold;
PTMResult.ind_fcc = PTMResult.classify==1 & PTMResult.FCCscore>scoreThreshold;
PTMResult.ind_hcp = PTMResult.classify==2 & PTMResult.HCPscore>scoreThreshold;
PTMResult.ind_dh = PTMResult.classify==3 & PTMResult.DHscore>scoreThreshold;
PTMResult.ind_ih = PTMResult.classify==4 & PTMResult.IHscore>scoreThreshold;

PTMResult.ind_chaos = ~(PTMResult.ind_fcc | PTMResult.ind_hcp | PTMResult.ind_dh | PTMResult.ind_ih);

%% plot
if DrawFlag
    figure;
    subplot(241)
    plotAtoms(PTMResult.xyz',PTMResult.FCCscore,struct('DrawColorbar',1))
    title('fcc score')
    subplot(242)
    plotAtoms(PTMResult.xyz',PTMResult.HCPscore,struct('DrawColorbar',1))
    title('hcp score')
    subplot(243)
    plotAtoms(PTMResult.xyz',PTMResult.DHscore,struct('DrawColorbar',1))
    title('dh score')
    subplot(244)
    plotAtoms(PTMResult.xyz',PTMResult.IHscore,struct('DrawColorbar',1))
    title('isoca score')
    subplot(245)
    histogram(PTMResult.FCCscore)
    subplot(246)
    histogram(PTMResult.HCPscore)
    subplot(247)
    histogram(PTMResult.DHscore)
    subplot(248)
    histogram(PTMResult.IHscore)
    
    figure;
    plotAtoms(PTMResult.xyz',PTMResult.classify,struct('DrawColorbar',1))
    title('total classification')
    
    figure;
    subplot(221)
    plotAtoms(PTMResult.xyz(PTMResult.classify==1,:)',ones(1,sum(PTMResult.classify==1)),struct('DrawColorbar',0,'Atomsize',5,'Colormap','gray'))
    plotAtoms(PTMResult.xyz(PTMResult.ind_fcc,:)',round(PTMResult.FCCscore(PTMResult.ind_fcc),1),struct('DrawColorbar',1,'Atomsize',30,'Colormap','parula'))
    title('fcc with fcc score')
    subplot(222)
    plotAtoms(PTMResult.xyz(PTMResult.classify==2,:)',ones(1,sum(PTMResult.classify==2)),struct('DrawColorbar',0,'Atomsize',5,'Colormap','gray'))
    plotAtoms(PTMResult.xyz(PTMResult.ind_hcp,:)',round(PTMResult.HCPscore(PTMResult.ind_hcp),1),struct('DrawColorbar',1,'Atomsize',30,'Colormap','parula'))
    title('hcp with hcp score')
    subplot(223)
    plotAtoms(PTMResult.xyz(PTMResult.classify==3,:)',ones(1,sum(PTMResult.classify==3)),struct('DrawColorbar',0,'Atomsize',5,'Colormap','gray'))
    plotAtoms(PTMResult.xyz(PTMResult.ind_dh,:)',round(PTMResult.DHscore(PTMResult.ind_dh),1),struct('DrawColorbar',1,'Atomsize',30,'Colormap','parula'))
    title('decahedron with dh score')
    subplot(224)
    plotAtoms(PTMResult.xyz(PTMResult.classify==4,:)',ones(1,sum(PTMResult.classify==4)),struct('DrawColorbar',0,'Atomsize',5,'Colormap','gray'))
    plotAtoms(PTMResult.xyz(PTMResult.ind_ih,:)',round(PTMResult.IHscore(PTMResult.ind_ih),1),struct('DrawColorbar',1,'Atomsize',30,'Colormap','parula','TickNumber',10))
    title('isocahedron with ih score')
        
    figure;
    plotAtoms(PTMResult.xyz(PTMResult.ind_chaos,:)',PTMResult.classify(PTMResult.ind_chaos),struct('DrawColorbar',1,'Atomsize',50,'Colormap','parula'))
    title('not ideally classified atoms')

end
end


%% child functions
function [vNew,fNew] = tesselate(v,f)

% tesselate all faces f from vertices v
% rMax = 1e-3;
% init
vNew = [v; zeros(size(v,1),3)];
fNew = zeros(size(f,1)*4,3);
vCount = size(v,1);
fCount = 0;

% main
for a0 = 1:size(f,1)
    inds = f(a0,:);
    
    vSub = v(inds,:);
    vAdd = [ ...
        vSub(1,:)+vSub(2,:);
        vSub(1,:)+vSub(3,:);
        vSub(2,:)+vSub(3,:);
        ] / 2;
    vInds = vCount + (1:3);
    
    fAdd = [ vInds;
        inds(1) vInds(1) vInds(2);
        inds(2) vInds(1) vInds(3);
        inds(3) vInds(2) vInds(3)];
    fInds = fCount + (1:4);
    
    vNew(vInds,:) = vAdd;
    fNew(fInds,:) = fAdd;

    vCount = vCount + 3;
    fCount = fCount + 4;
end
end


function [score] = calcScore(dxyzRef,dxyz,scoreRadiusRange)
r2max = scoreRadiusRange(2)^2;
score = 0;
for a0 = 1:size(dxyzRef,1)
    [d2min,ind] = min(sum((dxyzRef(a0,:) - dxyz).^2,2));
    if d2min < r2max % smaller than atom radius ( scoreRadiusRange(2) )
        score = score + ...
            min(max((scoreRadiusRange(2) - sqrt(d2min)) / (scoreRadiusRange(2) - scoreRadiusRange(1)),0),1);
    end
    dxyz(ind,:) = [];
end
end


function plotAtoms(pos,property,para)
arguments
    pos (3,:) double % YXZ
    property (1,:) double = []
    para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',0,'DrawAlphaShape',0,'FaceAlpha',0.5);
end

if isstruct(para)
    if ~isfield(para,'Colormap')
        para.Colormap = 'parula';
    end
    if ~isfield(para,'Atomsize')
        para.Atomsize = 50;
    end
    if ~isfield(para,'TickNumber')
        para.TickNumber = 'auto';
    end
    if ~isfield(para,'DrawColorbar')
        para.DrawColorbar = 0;
    end
    if ~isfield(para,'DrawAlphaShape')
        para.DrawAlphaShape = 0;
    end
    if ~isfield(para,'FaceAlpha')
        para.FaceAlpha = 0.5;
    end
elseif ischar(para)
    if strcmp(para,'triangle')
        para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',0,'DrawAlphaShape',1,'FaceAlpha',0.5);
    else
        try colormap(para)
            para = struct('Colormap',para,'Atomsize',50,'TickNumber','auto','DrawColorbar',0,'DrawAlphaShape',0,'FaceAlpha',0.5);
        catch
            para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',0,'DrawAlphaShape',0,'FaceAlpha',0.5);
        end
    end
elseif isnumeric(para)
    para = struct('Colormap','parula','Atomsize',round(para(1)),'TickNumber','auto','DrawColorbar',0,'DrawAlphaShape',0,'FaceAlpha',0.5);
else
    para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',0,'DrawAlphaShape',0,'FaceAlpha',0.5);
end

num=size(pos,2);

if nargin == 1
    property = ones(1,num);
end

if num~=length(property)
    error('Not equal number of atom properties inputed !')
end
if sum(isnan(property))>0
    error('There exists nan value in your input properties !')
end

maxproperty=max(property);
minproperty=min(property);
colormapNumber = min(101,length(unique(property)));
if ~strcmp(para.TickNumber,'auto')
    if para.TickNumber > colormapNumber
        colormapNumber = para.TickNumber;
    end
end
colorlist=feval(para.Colormap,colormapNumber);
property_rescaled=round(rescale(property)*(size(colorlist,1)-1)+1);
colorIndex=colorlist(property_rescaled,:);
if strcmp(para.TickNumber,'auto')
    tickLabelList = unique([round(property),floor(minproperty),ceil(maxproperty)]);
    para.TickNumber = length(tickLabelList);
    tickLabels = mat2cell(num2str(tickLabelList'),ones(1,para.TickNumber));
elseif isnumeric(para.TickNumber)
    tickLabels=cell(1,para.TickNumber);
    tickLabels{1} = num2str(minproperty);
    for i=2:para.TickNumber
        tickLabels{i}=num2str(minproperty+(maxproperty-minproperty)/(para.TickNumber-1)*(i-1));
    end
end

scatter3(pos(2,:)',pos(1,:)',pos(3,:)',para.Atomsize,colorIndex,"MarkerEdgeColor",'black',"MarkerFaceColor",'flat','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0.5);
colormap(colorlist)
if para.DrawColorbar
    tickPositions = 1/para.TickNumber/2 : 1/para.TickNumber : 1-1/para.TickNumber/2;
    colorbar('Ticks',tickPositions,'TickLabels',tickLabels,'Limits',[0,1])
end
if para.DrawAlphaShape
    hold on
    if num>3
        shp = alphaShape(pos(2,:)',pos(1,:)',pos(3,:)');
        a = criticalAlpha(shp,'one-region');
        shp2 = alphaShape(pos(2,:)',pos(1,:)',pos(3,:)',a*2);
        plot(shp2,'FaceAlpha',para.FaceAlpha,'EdgeAlpha',0.5)
    elseif num==3
        patch(pos(2,:)',pos(1,:)',pos(3,:)','red','FaceAlpha',0.5)
    elseif num<=2
        fprintf('Not enough atoms for Triangle mesh drawing!\n')
    end
end

axis equal image vis3d

xlabel 'X'
ylabel 'Y'
zlabel 'Z'
hold on

end