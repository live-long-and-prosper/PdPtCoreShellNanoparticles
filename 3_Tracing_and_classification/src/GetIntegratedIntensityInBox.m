function [boxIntensity_total,varargout] = GetIntegratedIntensityInBox(vol,pos,boxDiameter,interpMethod)
%% written by Zezhou Li, Peking University, Feb. 2023
%% Advanced Imaging Group, College of Chemistry and Molecular Engineering, Peking University, Beijing, China
%% Contact:  Website: https://www.chem.pku.edu.cn/jhzhou/  &  Email: jhzhou@pku.edu.cn
arguments
    vol(:,:,:) single
    pos(3,:) single
    boxDiameter(1,1) single {mustBePositive, mustBeInteger} = 7
    interpMethod(1,:) char {mustBeMember(interpMethod,{'linear','cubic'})} = 'cubic'
end

num = size(pos,2);
volSize = size(vol);
F = griddedInterpolant({1:volSize(1),1:volSize(2),1:volSize(3)},vol,interpMethod);
boxIntensity_total = zeros(1,num);
if nargout == 1
    for i = 1:num
        tempPos = pos(:,i);
        vec1 = linspace(tempPos(1)-boxDiameter/2+0.5,tempPos(1)+boxDiameter/2-0.5,boxDiameter);
        vec2 = linspace(tempPos(2)-boxDiameter/2+0.5,tempPos(2)+boxDiameter/2-0.5,boxDiameter);
        vec3 = linspace(tempPos(3)-boxDiameter/2+0.5,tempPos(3)+boxDiameter/2-0.5,boxDiameter);
        tempBox = F({vec1,vec2,vec3});
        boxIntensity_total(i) = sum(tempBox,'all');
    end
else
    boxIntensity_list = zeros(boxDiameter^3,num);
    for i = 1:num
        tempPos = pos(:,i);
        vec1 = linspace(tempPos(1)-boxDiameter/2+0.5,tempPos(1)+boxDiameter/2-0.5,boxDiameter);
        vec2 = linspace(tempPos(2)-boxDiameter/2+0.5,tempPos(2)+boxDiameter/2-0.5,boxDiameter);
        vec3 = linspace(tempPos(3)-boxDiameter/2+0.5,tempPos(3)+boxDiameter/2-0.5,boxDiameter);
        tempBox = F({vec1,vec2,vec3});
        boxIntensity_list(:,i) = tempBox(:);
    end
    boxIntensity_total = sum(boxIntensity_list,1);
    varargout{1} = boxIntensity_list;
end


end
