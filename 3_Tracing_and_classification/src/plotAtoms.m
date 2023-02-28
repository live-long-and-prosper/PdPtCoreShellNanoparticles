function plotAtoms(posYXZ,property,para)
%% written by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn

%% examples
% figure;plotAtoms(randn(3,100),1.5,200)
% figure;plotAtoms(randn(3,100),'C_N',200)
% figure;plotAtoms(randn(3,100),1:100,'parula')
% figure;plotAtoms(randn(3,100),1.5,'mesh')
% figure;plotAtoms(randn(3,100),1.5,{50,'cyan'});plotAtoms(randn(3,50),3.1,{50,'red'})
% figure;plotAtoms(randn(3,100),1:100,'jet')
% figure;plotAtoms(randn(3,10),1:10,{50,'jet','mesh'})
% figure;plotAtoms(randn(3,20),2.5,{100,'black','mesh'})
% figure;plotAtoms(randn(3,20),1:20,{100,'black','summer','mesh'})
% figure;plotAtoms(randn(3,100),2,{0,'red','mesh'})
% figure;plotAtoms(randn(3,100),'type 1',{20,'red'});plotAtoms(randn(3,20),'type 2',{50,'blue'})

arguments
    posYXZ (3,:) double % default as YXZ, not XYZ
    property (1,:) = 1
    para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',0,'FaceAlpha',0.5,'Color','#EDB120');
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
        para.DrawColorbar = 1;
    end
    if ~isfield(para,'DrawAlphaShape')
        para.DrawAlphaShape = 0;
    end
    if ~isfield(para,'FaceAlpha')
        para.FaceAlpha = 0.5;
    end
    if ~isfield(para,'Color')
        para.Color = '#EDB120'; % orange
    end
elseif ischar(para) || isstring(para)
    if strcmp(para,'triangle') || strcmp(para,'mesh') % DrawAlphaShape
        para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',1,'FaceAlpha',0.5,'Color','#EDB120');
    else
        try colormap(para) % draw colormap
            para = struct('Colormap',para,'Atomsize',50,'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',0,'FaceAlpha',0.5,'Color','#EDB120');
        catch % select atom Color 
            para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',0,'FaceAlpha',0.5,'Color',para);
        end
    end
elseif isnumeric(para)
    if length(para) == 1
        para = struct('Colormap','parula','Atomsize',round(para(1)),'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',0,'FaceAlpha',0.5,'Color','#EDB120');
    elseif length(para) == 3
        para = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',0,'FaceAlpha',0.5,'Color',para);
    end
else
    para_temp = struct('Colormap','parula','Atomsize',50,'TickNumber','auto','DrawColorbar',1,'DrawAlphaShape',0,'FaceAlpha',0.5,'Color','#EDB120');
    if iscell(para)
        for i = 1:length(para)
            a = para{i};
            if ischar(a) || isstring(a)
                try colormap(a)
                    para_temp.Colormap = a;
                catch
                    if strcmp(a,'triangle') || strcmp(a,'mesh')
                        para_temp.DrawAlphaShape = 1;
                    else
                        para_temp.Color = a;
                    end
                end
            elseif isnumeric(a)
                if length(a) == 1
                    para_temp.Atomsize = a;
                elseif length(a) == 3
                    para_temp.Color = a;
                end
            end
        end
    end
    para = para_temp;
end

num=size(posYXZ,2);

if nargin == 1
    property = ones(1,num);
end

if para.Atomsize ~= 0
    if ~isnumeric(property) || length(property) == 1 || length(unique(property)) == 1
        scatter3(posYXZ(2,:)',posYXZ(1,:)',posYXZ(3,:)',para.Atomsize,"MarkerEdgeColor",'#4c4c4c',"MarkerFaceColor",para.Color,'MarkerFaceAlpha',0.75);
        if ~isempty(gca().Legend)
            legendCells = gca().Legend.String;
        else
            legendCells = cell(0);
        end
        num_legend = length(legendCells);
        if isnumeric(property)
            legendCells{num_legend+1} = num2str(property(1));
        else
            legendCells{num_legend+1} = property;
        end
        legend(legendCells,'AutoUpdate','off','Location','northwest')
    else
        if num~=length(property)
            error('Not equal number of atom properties inputed !')
        end
        if sum(isnan(property))>0
            error('There exists nan value in your input properties !')
        end
        
        maxproperty=max(property);
        minproperty=min(property);
        if minproperty<0 && maxproperty>0
            FlagSymmetricScaleBar = 1;
            property2 = [property, -property];
            colormapNumber = min(101,length(unique(property2)));
            if ~strcmp(para.TickNumber,'auto')
                if para.TickNumber > colormapNumber
                    colormapNumber = para.TickNumber;
                end
            end
            colorlist=feval(para.Colormap,colormapNumber);
            property_rescaled=round(rescale(property2)*(size(colorlist,1)-1)+1);
        else
            FlagSymmetricScaleBar = 0;
            colormapNumber = min(101,length(unique(property)));
            if ~strcmp(para.TickNumber,'auto')
                if para.TickNumber > colormapNumber
                    colormapNumber = para.TickNumber;
                end
            end
            colorlist=feval(para.Colormap,colormapNumber);
            property_rescaled=round(rescale(property)*(size(colorlist,1)-1)+1);
        end
        
        colorIndex=colorlist(property_rescaled(1:num),:);
        % plot scale bar symmetrically
        if FlagSymmetricScaleBar
            if strcmp(para.TickNumber,'auto')
                tickLabelList = unique([round(property2),floor(-maxproperty),ceil(maxproperty)]);
                para.TickNumber = min(length(tickLabelList),20);
                tickLabels = mat2cell(num2str(round(linspace(min(tickLabelList),max(tickLabelList),para.TickNumber),1)'),ones(1,para.TickNumber));
            elseif isnumeric(para.TickNumber)
                tickLabels=cell(1,para.TickNumber);
                tickLabels{1} = num2str(-maxproperty);
                for i=2:para.TickNumber
                    tickLabels{i}=num2str(-maxproperty+2*maxproperty/(para.TickNumber-1)*(i-1));
                end
            end
        else
            if strcmp(para.TickNumber,'auto')
                tickLabelList = unique([round(property),floor(minproperty),ceil(maxproperty)]);
                para.TickNumber = min(length(tickLabelList),20);
                tickLabels = mat2cell(num2str(round(linspace(min(tickLabelList),max(tickLabelList),para.TickNumber),1)'),ones(1,para.TickNumber));
            elseif isnumeric(para.TickNumber)
                tickLabels=cell(1,para.TickNumber);
                tickLabels{1} = num2str(minproperty);
                for i=2:para.TickNumber
                    tickLabels{i}=num2str(minproperty+(maxproperty-minproperty)/(para.TickNumber-1)*(i-1));
                end
            end
        end
        
        scatter3(posYXZ(2,:)',posYXZ(1,:)',posYXZ(3,:)',para.Atomsize,colorIndex,"MarkerEdgeColor",'#4c4c4c',"MarkerFaceColor",'flat','MarkerFaceAlpha',0.75);
        colormap(colorlist)
        if ~isempty(gca().Legend)
            legendCells = gca().Legend.String;
        else
            legendCells = cell(0);
        end
        num_legend = length(legendCells);
        legendCells{num_legend+1} = 'colormap atoms';
        legend(legendCells,'AutoUpdate','off','Location','northwest')
        if para.DrawColorbar
            tickPositions = 1/para.TickNumber/2 : 1/para.TickNumber : 1-1/para.TickNumber/2;
            colorbar('Ticks',tickPositions,'TickLabels',tickLabels,'Limits',[0,1])
        end
    end
end

if para.DrawAlphaShape
    hold on
    if num>3
        shp = alphaShape(posYXZ(2,:)',posYXZ(1,:)',posYXZ(3,:)');
        a = criticalAlpha(shp,'one-region');
        shp2 = alphaShape(posYXZ(2,:)',posYXZ(1,:)',posYXZ(3,:)',a*2);
        plot(shp2,'FaceColor',rand(1,3),'FaceAlpha',para.FaceAlpha,'EdgeAlpha',0.5)
        if ~isempty(gca().Legend)
            legendCells = gca().Legend.String;
        else
            legendCells = cell(0);
        end
        num_legend = length(legendCells);
        legendCells{num_legend+1} = 'surface';
        legend(legendCells,'AutoUpdate','off','Location','northwest')
    elseif num==3
        patch(posYXZ(2,:)',posYXZ(1,:)',posYXZ(3,:)','red','FaceAlpha',0.5)
    else % num<=2
        fprintf('Not enough atoms for Triangle mesh drawing!\n')
    end
end


axis equal image vis3d

xlabel 'X'
ylabel 'Y'
zlabel 'Z'

hold on

end
