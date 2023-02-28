function Properties = CalculatePropertiesInCoreShell(AtomPos,AtomType,volumeSize,Res,exfoliationStartPosition,stepSize,doPlot,alphaShapePara)
%% written by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn
arguments
    AtomPos(3,:) double
    AtomType(1,:) double
    volumeSize(1,:) double % size of the reconstruction volume
    Res(1,1) double  = 0.3434 % voxel size
    exfoliationStartPosition(1,:) char {mustBeMember(exfoliationStartPosition,{'clean Pd','clean Pt','surface'})} = 'clean Pd'
    stepSize(1,1) double {mustBeInteger,mustBePositive} = 1
    doPlot(1,1) double {mustBeNonnegative} = 1 % 0:do not plot, 1:plot and auto copy to clipboard, 2: plot and auto save
    alphaShapePara(1,1) double = 2
end

%% initialize
% the type label of core and shell atoms are set defaultly as 1 and 2,respectively
stepType = 'voxels';
num = length(AtomType);
[~,~,threshold] = CalculateLatticeConstant(AtomPos,Res,'fcc'); % calculate Nearest Neighbor threshold of Pd@Pt nanoparticles
coor = CalculateCoordinationNumber(AtomPos,threshold,Res,AtomType);
close all
Properties = struct();
Properties.AtomPos = AtomPos;
Properties.AtomType = AtomType;
Properties.Res = Res;
Properties.stepSize = stepSize;
Properties.stepType = stepType;
Properties.exfoliationStartPosition = exfoliationStartPosition;
Properties.alphaShapePara = alphaShapePara;
Properties.labeledVolume = zeros(volumeSize); % labeled volume
Properties.Depth = 10;
Properties.AtomIndexByDepth = zeros(1,num);

% AtomPos = AtomPos+(1+volumeSize')/2;
shp = getSuitableAlphaShape(AtomPos([2,1,3],:)',alphaShapePara);
[qx,qy,qz] = meshgrid(1:volumeSize(1),1:volumeSize(2),1:volumeSize(3));
tf = inShape(shp,qx,qy,qz);
tf_total = imdilate(tf,strel('cube',2*3+1));

%% calculate core atom concentration in each layer
AtomPos_positive_round=round(AtomPos);
if strcmp(stepType,'voxels')
    if strcmp(exfoliationStartPosition,'clean Pd')
        % use CN == 12 as criterion
        purePdlist = (coor(1,:)==12) & (AtomType==1);
        pos_purePd_positive = AtomPos(:,purePdlist);
        shp = getSuitableAlphaShape(pos_purePd_positive([2,1,3],:)',alphaShapePara);
        [qx,qy,qz] = meshgrid(1:volumeSize(1),1:volumeSize(2),1:volumeSize(3));
        tf = inShape(shp,qx,qy,qz);
        tf_dilate = imdilate(tf,strel('cube',2*2+1)).*tf_total; % dilate 2 voxels
        Properties.labeledVolume = tf_dilate;
        
        relativedepth = 1;
        DoErode = 1;
        tf_temp = tf_dilate;
        while DoErode
            tf_temp = imerode(tf_temp,strel('cube',2*stepSize+1)); % strel 2r+1
            if sum(tf_temp(:))==0 || sum(tf_temp(sub2ind(volumeSize,AtomPos_positive_round(1,:),AtomPos_positive_round(2,:),AtomPos_positive_round(3,:))))==0
                DoErode=0;
            else
                relativedepth = relativedepth+1;
                Properties.labeledVolume = Properties.labeledVolume+tf_temp;
            end
        end  
        DoDilate = 1;
        tf_temp = tf_dilate;
        while DoDilate
            currentV = tf_temp;
            tf_temp = imdilate(tf_temp,strel('cube',2*stepSize+1)).*tf_total; % strel 2r+1
            dilatedVolume = tf_temp-currentV;
            if sum(dilatedVolume(:))==0 || sum(dilatedVolume(sub2ind(volumeSize,AtomPos_positive_round(1,:),AtomPos_positive_round(2,:),AtomPos_positive_round(3,:))))==0
                DoDilate = 0;
            else
                relativedepth = relativedepth+1;
                Properties.labeledVolume = tf_temp+Properties.labeledVolume;
            end
        end
        
        Properties.Depth = relativedepth;
        Properties.labeledVolume = (relativedepth+1-Properties.labeledVolume).*tf_total; % define where depth=1
        figure;sliceViewer(Properties.labeledVolume);
        
    elseif strcmp(exfoliationStartPosition,'clean Pt')
        % use CN == 12 as criterion
        purePtlist = (coor(2,:)==12) & (AtomType==2);
        pos_temp = AtomPos(:,~purePtlist);
        shp = getSuitableAlphaShape(pos_temp([2,1,3],:)',alphaShapePara);
        [qx,qy,qz] = meshgrid(1:volumeSize(1),1:volumeSize(2),1:volumeSize(3));
        tf = inShape(shp,qx,qy,qz);
        tf_dilate = imerode(tf,strel('cube',2*2+1)).*tf_total;
        Properties.labeledVolume = tf_dilate;
        
        relativedepth = 1;
        DoErode = 1;
        tf_temp = tf_dilate;
        while DoErode
            tf_temp=imerode(tf_temp,strel('cube',2*stepSize+1));
            if sum(tf_temp(:))==0 || sum(tf_temp(sub2ind(volumeSize,AtomPos_positive_round(1,:),AtomPos_positive_round(2,:),AtomPos_positive_round(3,:))))==0
                DoErode=0;
            else
                relativedepth = relativedepth + 1;
                Properties.labeledVolume=Properties.labeledVolume+tf_temp;
            end
        end  
        DoDilate = 1;
        tf_temp = tf_dilate;
        while DoDilate
            currentV = tf_temp;
            tf_temp = imdilate(tf_temp,strel('cube',2*stepSize+1)).*tf_total;
            dilatedVolume = tf_temp-currentV;
            if sum(dilatedVolume(:))==0 || sum(dilatedVolume(sub2ind(volumeSize,AtomPos_positive_round(1,:),AtomPos_positive_round(2,:),AtomPos_positive_round(3,:))))==0
                DoDilate=0;
            else
                relativedepth=relativedepth+1;
                Properties.labeledVolume = tf_temp+Properties.labeledVolume;
            end
        end
        
        Properties.Depth = relativedepth;
        Properties.labeledVolume = (relativedepth+1-Properties.labeledVolume).*tf_total;
    
    elseif strcmp(exfoliationStartPosition,'surface')
        Properties.labeledVolume = tf_total;
        relativedepth = 1;
        DoErode = 1;
        tf_temp = tf_total;
        while DoErode
            tf_temp = imerode(tf_temp,strel('cube',2*stepSize+1));
            if sum(tf_temp(:))==0 || sum(tf_temp(sub2ind(volumeSize,AtomPos_positive_round(1,:),AtomPos_positive_round(2,:),AtomPos_positive_round(3,:))))==0
                DoErode=0;
            else
                relativedepth=relativedepth+1;
                Properties.labeledVolume=Properties.labeledVolume+tf_temp;
            end
        end
        
        Properties.Depth = relativedepth;
        Properties.labeledVolume = (relativedepth+1-Properties.labeledVolume).*tf_total;
        figure;sliceViewer(Properties.labeledVolume);
    
    else
        error('wrong input parameters !\n')
    end
end

%% calculate other data...
Properties.meanCoreAtomDensityByDepth = zeros(1,Properties.Depth);
Properties.meanCoorNumByDepth = zeros(4,Properties.Depth); % four type of CN in each layer（1-1，2-2，2 around 1，1 around 2) (core-core，shell-shell，shell atoms around core atom，core atoms around shell atom）
Properties.atomNumberByDepth = zeros(3,Properties.Depth);
for i = 1:Properties.Depth
    V_temp = (Properties.labeledVolume==i).*tf_total;
    posIndex = logical(V_temp(sub2ind(volumeSize,AtomPos_positive_round(1,:),AtomPos_positive_round(2,:),AtomPos_positive_round(3,:))));
    Properties.AtomIndexByDepth(posIndex) = i;
    Properties.meanCoreAtomDensityByDepth(i) = 1/(1+sum(posIndex & AtomType==2)/sum(posIndex & AtomType==1));
    Properties.meanCoorNumByDepth(1,i) = mean(coor(1,(posIndex & AtomType==1)));
    Properties.meanCoorNumByDepth(2,i) = mean(coor(2,(posIndex & AtomType==2)));
    Properties.meanCoorNumByDepth(3,i) = mean(coor(2,(posIndex & AtomType==1)));
    Properties.meanCoorNumByDepth(4,i) = mean(coor(1,(posIndex & AtomType==2)));
    Properties.atomNumberByDepth(1,i) = sum(Properties.AtomIndexByDepth==i);
    Properties.atomNumberByDepth(2,i) = sum(Properties.AtomIndexByDepth==i & Properties.AtomType==1);
    Properties.atomNumberByDepth(3,i) = sum(Properties.AtomIndexByDepth==i & Properties.AtomType==2);
end

[~,Properties.locationOfInterface] = min(abs(Properties.meanCoreAtomDensityByDepth-0.5));

Properties.Depth = sum(~isnan(Properties.meanCoreAtomDensityByDepth));
Properties.meanCoreAtomDensityByDepth = Properties.meanCoreAtomDensityByDepth(1:Properties.Depth);
Properties.meanCoorNumByDepth = Properties.meanCoorNumByDepth(:,1:Properties.Depth);
Properties.atomNumberByDepth = Properties.atomNumberByDepth(:,1:Properties.Depth);

% interpolate and define interface as 5% ~ 95%
x = 1:Properties.Depth;
y = Properties.meanCoreAtomDensityByDepth;
step = 0.01;
xq = 1:step:Properties.Depth;
meanPdDensity_interp = interp1(x,y,xq,'makima');
% find the depth position of Pd conc. == 95% 5% 99% 1%
Properties.DepthLocation_interface = xq([find(meanPdDensity_interp <= 0.95,1,'first'),find(meanPdDensity_interp >= 0.05,1,'last'),find(meanPdDensity_interp <= 0.99,1,'first'),find(meanPdDensity_interp >= 0.01,1,'last')]);
[~,id] = min(abs(meanPdDensity_interp-0.5));
Properties.DepthLocation_medium = xq(id);
if strcmp(stepType,'voxels')
    Properties.Width_interface = (Properties.DepthLocation_interface(2)-Properties.DepthLocation_interface(1))*Properties.Res*Properties.stepSize;
else
    Properties.Width_interface = nan;
end

%% plot
if doPlot
    close all
    figure;tiledlayout(2,1,'TileSpacing','none')
    set(gcf,'Position',[100,100,700,600])
    step = Properties.Res*Properties.stepSize;
    nexttile
%     x = [Properties.DepthLocation_interface(1)-0.5,Properties.DepthLocation_interface(1)-0.5,Properties.DepthLocation_interface(2)+0.5,Properties.DepthLocation_interface(2)+0.5];
    x1 = [Properties.DepthLocation_interface(3),Properties.DepthLocation_interface(3),Properties.DepthLocation_interface(1),Properties.DepthLocation_interface(1)];
    x2 = [Properties.DepthLocation_interface(1),Properties.DepthLocation_interface(1),Properties.DepthLocation_interface(2),Properties.DepthLocation_interface(2)];
    x3 = [Properties.DepthLocation_interface(2),Properties.DepthLocation_interface(2),Properties.DepthLocation_interface(4),Properties.DepthLocation_interface(4)];
    y = [0,1,1,0];
    patch('XData',x1*step,'YData',y,'FaceColor','#f9eda6','FaceAlpha',1,'EdgeColor','none')
    hold on
    patch('XData',x2*step,'YData',y,'FaceColor','#bbd597','FaceAlpha',1,'EdgeColor','none')
    patch('XData',x3*step,'YData',y,'FaceColor','#c7e6ea','FaceAlpha',1,'EdgeColor','none')
    plot((1:Properties.Depth)*step,Properties.meanCoreAtomDensityByDepth,'LineWidth',2)
    legend({'99%~95%','95%~5%','5%~1%','Pd'},'AutoUpdate','off','Location','northeast')
    xlim([1,Properties.Depth]*step)
%     xlabel('depth (Å)')
    ylabel('concentration')
    set(gca,'FontSize',14,'FontName','Arial','fontweight','bold','lineWidth',1.5,'BoxStyle','full','Box','on');
    title(['start from ',Properties.exfoliationStartPosition, ', layer thickness: ',num2str(Properties.stepSize),' ',Properties.stepType],'FontSize',16,'FontName','Arial','fontweight','bold')
    
    nexttile
    patch('XData',x1*step,'YData',y*max(Properties.atomNumberByDepth(1,:))*1.2,'FaceColor','#f9eda6','FaceAlpha',1,'EdgeColor','none')
    hold on
    patch('XData',x2*step,'YData',y*max(Properties.atomNumberByDepth(1,:))*1.2,'FaceColor','#bbd597','FaceAlpha',1,'EdgeColor','none')
    patch('XData',x3*step,'YData',y*max(Properties.atomNumberByDepth(1,:))*1.2,'FaceColor','#c7e6ea','FaceAlpha',1,'EdgeColor','none')
    plot(1:Properties.Depth*step,Properties.atomNumberByDepth(1,:),'k-','LineWidth',2)
    plot(1:Properties.Depth*step,Properties.atomNumberByDepth(2,:),'Color','#0072BD','LineWidth',2)
    plot(1:Properties.Depth*step,Properties.atomNumberByDepth(3,:),'Color','#D95319','LineWidth',2)
    legend({'99%~95%','95%~5%','5%~1%','total','Pd','Pt'},'AutoUpdate','off','Location','northeast')
    xlim([1,Properties.Depth]*step)
    ylim([0,max(Properties.atomNumberByDepth(1,:))*1.1])
    xlabel('depth (Å)')
    ylabel('atom number')
    set(gca,'FontSize',14,'FontName', 'Arial','fontweight','bold','lineWidth',1.5,'BoxStyle','full','Box','on');
    
    % copygraphics(gcf,'ContentType','vector','BackgroundColor','none','Resolution',600)
%     exportgraphics(gcf,['D:\start from',Properties.exfoliationStartPosition,'-step ',num2str(Properties.stepSize),' ',Properties.stepType,'.emf'],'ContentType','vector','BackgroundColor','none','Resolution',600)
end

end

%% child functions
function coor=CalculateCoordinationNumber(pos,threshold,varargin)
    num=size(pos,2);
    
    if numel(varargin)>0
        Res=varargin{1};
    else
        Res=1;
    end
    
    if numel(varargin)>1
        type=varargin{2};
        if num~=length(type)
            error('Please Input equal number of atom types !\n')
        end
        if ~isvector(type)
            error('Please Input right format of atom types (i.e. Input a vector !)\n')
        end
    else
        type=ones(1,num);
    end
    
    typeList=unique(type);
    typeNum=length(typeList);
    coor=zeros(typeNum+1,num);
    for i=1:num
        difference=pos-pos(:,i);
        distance=sqrt(difference(1,:).^2+difference(2,:).^2+difference(3,:).^2);
        selectlist=(distance<=threshold/Res) & (distance>0.1);
        for kkk=1:typeNum
            coor(kkk,i)=sum(selectlist & (type==kkk));
        end
    end
    coor(end,:)=sum(coor(1:typeNum,:),1);
    
    figure;
    for kkk=1:typeNum
        subplot(1,typeNum+1,kkk);histogram(coor(kkk,:));title(['CN count to atoms type ',num2str(kkk)])
    end
    subplot(1,typeNum+1,typeNum+1);histogram(coor(typeNum+1,:));title('CN count total')
end

function shp = getSuitableAlphaShape(posXYZ,alphaPara)
arguments
    posXYZ(:,:) % n*2 or n*3
    alphaPara(1,1) double = 2
end

if size(posXYZ,2) == 2
    shp = alphaShape(posXYZ(:,1),posXYZ(:,2));
    a = criticalAlpha(shp,'one-region');
    shp = alphaShape(posXYZ(:,1),posXYZ(:,2),a*alphaPara);
elseif size(posXYZ,2) == 3
    shp = alphaShape(posXYZ(:,1),posXYZ(:,2),posXYZ(:,3));
    a = criticalAlpha(shp,'one-region');
    shp = alphaShape(posXYZ(:,1),posXYZ(:,2),posXYZ(:,3),a*alphaPara);
else
    error('please input point positions in right forum !')
end

end