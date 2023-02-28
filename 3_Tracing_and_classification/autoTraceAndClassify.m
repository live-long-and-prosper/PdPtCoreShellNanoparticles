%% 0. load 3D reconstructed volume
load('\2_Final_reconstruction_volume\reconstruction_EPB.mat') % use EPB as example
res = 0.3434; % Angstrom, resolution of voxel (voxel size)

%% 1. compute local maximum and remove non-atoms
volSize = size(vol);
roughPos = FindLocalMaximums(vol);

%% 2. polynomial fit and exclude too closed atoms
Dmin = 2.21; % Angstrom, minimum distance chose for Pd@Pt core-shell nanoparticles
% fit sub-voxel resolutions of atomic coordinates, parallel computing is
% auto set to save time, which may need rather large memory...
[pos0,tooClosePos] = PolynomialFit3DAtoms(vol,roughPos,Dmin,res); % a few miniutes to hours to wait...
% plot result
figure;plotAtoms(pos0,'good',{'red',20})
hold on;plotAtoms(tooClosePos,'too close',{'blue',50})

%% 3. after manual check, classify integrated intensity in the box by k-means
[type,pos] = IntegratedAtomClassification(pos0,res,vol,2,'boxIntensity','kmeans',1);
% plot result
figure;plotAtoms(pos,type)