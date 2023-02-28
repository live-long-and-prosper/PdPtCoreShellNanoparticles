function [roughPos,boxIntensity] = FindLocalMaximums(vol,localMaximumBoxSize,intensityBoxRadius)
%% modified from Jianwei Miao's previous AET work by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn
arguments
    vol(:,:,:) double
    localMaximumBoxSize(1,1) double = 5 % default diameter as 5 (5 = 2*2+1)
    intensityBoxRadius(1,1) double {mustBeInteger,mustBeNonnegative} = 1 % default diameter as 3 (3*3*3 box)
end

%% get list of local maximum, then sort the integrated intensity in box from low to high
[~, ~, ~, sphere] = create_box(localMaximumBoxSize); % 5为BoxSize1
M = find_maxima3D(vol,sphere);
% M(RecVol<10)=0;
M=M.*(My_paddzero(ones(size(vol)-10),size(vol))); % set local maximums which near boundaries (defined as 5 voxels) all to zero
PossiblePosNum=sum(M(:));
fprintf('There are possiblely %d atoms in this particle.\n',PossiblePosNum);
PeaksBox3=zeros(size(vol));
% get integrated intensity of the box centered at the local maximum
ind_find=find(M);
[yy,xx,zz] = ind2sub(size(M),ind_find);
tempPeakList=zeros(1,PossiblePosNum);
for kkk=1:PossiblePosNum
    IntensityBox3 = sum(vol(yy(kkk)-intensityBoxRadius:yy(kkk)+intensityBoxRadius, ...
        xx(kkk)-intensityBoxRadius:xx(kkk)+intensityBoxRadius, ...
        zz(kkk)-intensityBoxRadius:zz(kkk)+intensityBoxRadius),'all');
    PeaksBox3(yy(kkk),xx(kkk),zz(kkk)) = IntensityBox3;
    tempPeakList(kkk) = IntensityBox3;
end
[PeakList, ind] = sort(tempPeakList,'descend');
PeakPosList=[yy(ind),xx(ind),zz(ind)]';

%% manual select a threshold and delete non-atoms (weak value)
NumOfBins=100;
[Count,BinEdges] = histcounts(PeakList,NumOfBins);
BinLocation=(BinEdges(1:end-1)+BinEdges(2:end))/2;
TF = islocalmin(smooth(Count,5)); % find local min after smooth
K=find(TF);
figure;plot(BinLocation,Count,BinLocation(TF),Count(TF),'r*');title(['histogram of Local Maximum in Reconstruction, with ',num2str(length(K)),' local minimums exist'])
flg=input('\n choose which local minimum as threshold? \n(input interger 1,2,3,etc, default as 1; if not satisfy with this result, input 0)');
if isempty(flg)
    fprintf('Use the first maximum By Default\n');
    flg=1; %  select the first one as default
    tempLocation=K(flg);
    threshold=BinLocation(tempLocation);
elseif flg==0
    figure;histogram(PeakList,NumOfBins)
    threshold=round(input('Manual choose a threshold?\n'));
% % elseif flg<1 || flg>length(K)
    fprintf('Warning! Use the first minimum in default\n');
    flg=1;
    tempLocation=K(flg);
    threshold=BinLocation(tempLocation);
else
    flg=round(flg);
    tempLocation=K(flg);
    threshold=BinLocation(tempLocation);
end
Th = threshold/2;
figure;histogram(PeakList,NumOfBins);hold on;plot([threshold,threshold],[1,max(Count)],"LineWidth",1);plot([Th,Th],[1,max(Count)],"LineWidth",1);hold off
figure;plot(PeakList);hold on;plot([1,PossiblePosNum],[threshold,threshold],"LineWidth",1);plot([1,PossiblePosNum],[Th,Th],"LineWidth",1);hold off

NumofRoughPos = sum(PeakList>threshold);
fprintf('There are approximately %d atoms in this particle.\n',NumofRoughPos);
RoughAtomPosList = PeakPosList(:,1:NumofRoughPos);
RoughAtomintensityList = PeakList(1:NumofRoughPos);

%% calculate 3D mask and delete all maximums outside the mask
paddingSize = 0;
mask3D = getMask3D(vol, paddingSize);

posIND = sub2ind(size(vol),RoughAtomPosList(1,:),RoughAtomPosList(2,:),RoughAtomPosList(3,:));
RemoveList = logical(mask3D(posIND));
% RemoveList = ones(1,size(RoughAtomPosList,2));
% for i = 1:NumofRoughPos
%     if mask3D(RoughAtomPosList(1,i),RoughAtomPosList(2,i),RoughAtomPosList(3,i)) == 0
%         RemoveList(i) = 0;
%     end
% end
fprintf('%d atom postions are deleted because they are far away from the bulk particle\n',sum(~RemoveList))
roughPos = RoughAtomPosList(:,RemoveList);
boxIntensity = RoughAtomintensityList(RemoveList);
fprintf('%d local maximum postions found finally\n',size(roughPos,2))

end

%% children functions
% create box
function [boxCoordinates, BoxCenter, BoxRadius, sphere] = create_box(BoxSize)

%M. Bartels, UCLA, 2014
%small function to create box coordinates and correpsonding spherical mask
%use odd BoxSize

BoxCenter = (BoxSize+1)/2;
BoxRadius = (BoxSize-1)/2;

%box coordinate systems
[boxX,boxY,boxZ]=ndgrid(-BoxRadius:BoxRadius,-BoxRadius:BoxRadius,-BoxRadius:BoxRadius);

%calculate spherical mask
sphere = sqrt(boxX.^2+boxY.^2+boxZ.^2)<=BoxSize/2;

boxCoordinates.x = boxX;
boxCoordinates.y = boxY;
boxCoordinates.z = boxZ;

end

% find local maxima
function M = find_maxima3D(mat3d,mask)
%M. Bartels, UCLA, 2014
%find all maxima of a 3D matrix
%with a local neighbourhood defined by mask

s = size(mask,1);
c = (s+1)/2;
mask=mask>0;
mask(c,c,c) = false; %exclude the pixel itself from the neighbourhood

%for each voxel, assign the maximum in the neighbourhood
mat3d_neighbours = imdilate(mat3d,mask);

%a maximum is where the matrix is larger than all neighbors
M = mat3d > mat3d_neighbours; 

end

% get 3D mask
function mask3D = getMask3D(vol, paddingSize, otsuParameter, gaussfiltSigma)
arguments
    vol(:,:,:) double
    paddingSize(1,1) double {mustBeInteger} = 0
    otsuParameter(1,1) double {mustBePositive} = 0.5
    gaussfiltSigma(1,1) double {mustBePositive} = 3
end

vol_filt = rescale(imgaussfilt3(vol,gaussfiltSigma));
if paddingSize > 0
    vol_filt = imdilate(vol_filt,strel('sphere',paddingSize));
elseif paddingSize < 0
    vol_filt = imerode(vol_filt,strel('sphere',-paddingSize));
end

[Count,~] = histcounts(vol_filt);
otsu_th=otsuthresh(Count);
BW = imbinarize(vol_filt,otsu_th*otsuParameter);
BW=bwareaopen(BW,round(length(vol_filt(:))/100));
BW=imclose(BW,strel('sphere',5));
BW=imcomplement(bwareaopen(imcomplement(BW),round(length(vol_filt(:))/100)));
while length(table2array(regionprops3(bwconncomp(BW),"Volume")))>1
    BW=imopen(BW,strel('sphere',10));
    BW=bwareaopen(BW,round(length(vol_filt(:))/100));
end
mask3D = BW;

end
