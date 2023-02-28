%% calculate PDF, mean bond length and Nearest Neighbor threshold (first-nearest-neighbor shell distance)
clear
% load atom position and type first, use PB as example
load('\4_Final_coordinates\coordinates_PB.mat')
doPlot = 1;
[LatticeConst, MeanBondLength, Threshold, PDF] = CalculateLatticeConstant(pos,res,'fcc',doPlot);

%% analysis of core-shell interface
clear
% load atom position and type first, use PB as example
load('\4_Final_coordinates\coordinates_PB.mat')
res = 0.3434;
% load 3D reconstruction to get the size of volume, use PB as example
load('\2_Final_reconstruction_volume\reconstruction_PB.mat')
volumeSize = size(vol);
% exfoliation start from...
exfoliationStartPosition = 'clean Pd'; % 'clean Pd','clean Pt','surface'
stepSize = 3; % thickness of each layer
doPlot = 1; % plot the result and auto copy to clipboard
Properties = CalculatePropertiesInCoreShell(pos,type,volumeSize,res,exfoliationStartPosition,stepSize,doPlot);

%% polyhedral template matching
clear
% load atom position and type first, use PB as example
load('\4_Final_coordinates\coordinates_PB.mat')
res = 0.3434;
% polyhedral template matching (need a few minites)
PTMResult = PTM01(pos, res);
PTMResult = PTM02(PTMResult, 1, 1, 5, 5);
% draw results clearly
figure;
subplot(131)
plotAtoms(PTMResult.xyz(PTMResult.ind_fcc,:)','fcc',{'#EDB120',10})
plotAtoms(PTMResult.xyz(PTMResult.ind_hcp,:)','hcp',{'#0072BD',30})
plotAtoms(PTMResult.xyz(PTMResult.ind_dh,:)','dh',{'#D95319',30})
subplot(132)
plotAtoms(PTMResult.xyz(PTMResult.ind_hcp,:)','hcp',{'#0072BD',30})
plotAtoms(PTMResult.xyz(PTMResult.ind_dh,:)','dh',{'#D95319',50})
subplot(133)
plotAtoms(PTMResult.xyz(PTMResult.ind_chaos,:)','others',10)


%% strain calculation
clear
% load atom position and type first, use PB as example
load('\4_Final_coordinates\coordinates_PB.mat')
res = 0.3434;
% rotate coordinates to zone axis direction
rotateEulerAngle = [-2.5,-8,6]; % for PB
% rotateEulerAngle = [35.5,7.5,-4.5]; % for EPB
Phi = rotateEulerAngle(1);
Theta = rotateEulerAngle(3);
Psi = rotateEulerAngle(2);
R1=[cosd(Phi) -sind(Phi) 0;
        sind(Phi)  cosd(Phi) 0;
           0         0     1];
R2=[cosd(Theta)   0   sind(Theta);
         0       1        0    ;
    -sind(Theta)  0   cosd(Theta)];
R3=[1      0         0     ;
    0  cosd(Psi) -sind(Psi);
    0  sind(Psi)  cosd(Psi)];
RotationMatrix = (R2*R3*R1);
pos_rotated = RotationMatrix*pos;
figure;plotAtoms(pos_rotated,type,10)
% find fivefold axis atoms by polyhedral template matching
PTMResult = PTM01(pos_rotated, res);
PTMResult = PTM02(PTMResult, 1, 1, 5, 5);
pos_fivefold = pos_rotated(:,PTMResult.ind_dh);
% move fivefold axis to center
center_fivefold = mean(pos_fivefold,2);
pos_center = pos_rotated - center_fivefold;
figure;plotAtoms(pos_center,type,10)

% fit fcc lattice
sStrain = strain_01(pos_center*res,type);
% calculate displacements
sStrain = strain_02(sStrain);
% calculate strain
sStrain = strain_03(sStrain);
% draw global strain
strain_04(sStrain, 0)
% draw local strain for 5 grains
strain_04(sStrain, 1)

