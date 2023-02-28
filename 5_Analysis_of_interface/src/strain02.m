function [sStrain] = strain02(sStrain)
tic

% Colin Ophus - 2022 April 2
% This set of scripts calculates 3D strain tensors for pentagonal bipyramid
% (PB) structures.  Uses shared lattice vectors between twin planes.

% 03 - Calculate the 3 displacement tensors (x,y,z) for each of the 5
% lattice sectors, using a best fit FCC lattice.  Assuming vertical lattice
% vector direction (110) is the same for all sectors.

% Also calculate the displacement arrays for the global best fit lattice.

sStrain.voxel_size = 0.5;
sStrain.vol_size = [200 200 140];  % in voxels!
% sStrain.vol_sigma = (2.75*1.0) / sStrain.voxel_size;  % scale KDE sigma with voxel size
sStrain.vol_sigma = (3.89*1.0) / sStrain.voxel_size;  % scale KDE sigma with voxel size

% Lattice constant for sample
sStrain.a_lat = sqrt(2)*2.75;
sStrain.dens_mean = (4/sStrain.a_lat^3) * sStrain.voxel_size^3;

% Lattice vectors
dist_out_of_plane = sStrain.a_lat/sqrt(2);
dist_in_plane =  sStrain.a_lat*sqrt(3/2);


% Find the 5 reference lattices.
% Store each [4x3] lattice separately.
sStrain.lat_ref_fcc = zeros(4,3,5);
sStrain.lat_ref_inds = zeros(5,3);
theta = acos(1/3);
for a0 = 1:5
    % indices
    ind1 = 2 + a0;
    ind2 = 2 + mod(a0,5) + 1;
    ind3 = 2;
    sStrain.lat_ref_inds(a0,:) = [ind1 ind2 ind3];
    
    u0 = sStrain.lat(ind1,:);
    v0 = sStrain.lat(ind2,:);
    w0 = sStrain.lat(ind3,:);
    
    % find mid plane vector between u0 and v0
    uv = u0/norm(u0) + v0/norm(v0);
    uv = uv / norm(uv);
    
    % Generate coordinate system perpendicular to w0
    x = cross(uv,w0);
    y = cross(w0,x);
    x = x / norm(x);
    y = y / norm(y);
    
    % New FCC reference lattice
    u = y*cos(theta/2) + x*sin(theta/2);
    v = y*cos(theta/2) - x*sin(theta/2);
    w = w0 / norm(w0);
    
    % Scale to reference vector lengths
    u = u * dist_in_plane;
    v = v * dist_in_plane;
    w = w * dist_out_of_plane;

    % output
    sStrain.lat_ref_fcc(:,:,a0) = [sStrain.lat(1,:); u; v; w];
end


% % plot lattice sectors
% figure(56)
% clf
% hold on
% for a0 = 1:5
%     lat = sStrain.lat_ref_fcc(:,:,a0);
%     
%     p = [0 0 0;
%         lat(2,:);
%         lat(2,:)+ lat(3,:);
%         lat(3,:);
%         0 0 0];
%     line(p(:,2),p(:,1),p(:,3),...
%         'linewidth',2,'color','r')
% end
% hold off
% axis equal
% box on


% KDE - make kernel
x = mod((0:(sStrain.vol_size(1)-1)) + sStrain.vol_size(1)/2, ...
    sStrain.vol_size(1)) - sStrain.vol_size(1)/2;
y = mod((0:(sStrain.vol_size(2)-1)) + sStrain.vol_size(2)/2, ...
    sStrain.vol_size(2)) - sStrain.vol_size(2)/2;
z = mod((0:(sStrain.vol_size(3)-1)) + sStrain.vol_size(3)/2, ...
    sStrain.vol_size(3)) - sStrain.vol_size(3)/2;
x = reshape(x,[sStrain.vol_size(1) 1 1]);
y = reshape(y,[1 sStrain.vol_size(2) 1]);
z = reshape(z,[1 1 sStrain.vol_size(3)]);
k = exp((x.^2 + y.^2 + z.^2)/(-2*sStrain.vol_sigma^2));
k(:) = k / sum(k(:));
k(:) = k / sStrain.dens_mean;
k(:) = fftn(k);


% r = ceil(2*sStrain.vol_sigma);
% v = -r:r;


% init
% best fit fcc lattices
sStrain.lat_fcc_fit = zeros(4,3,5);
% 5 displacement arrays for 5 deca sectors
% dim 4 - [dx,dy,dz,count]
% dim 5 - 5 sectors
vol = zeros(sStrain.vol_size);
vol_count = zeros(sStrain.vol_size);
sStrain.vol_disp = zeros( ...
    sStrain.vol_size(1),...
    sStrain.vol_size(2),...
    sStrain.vol_size(3),...
    4,...
    5);
sStrain.vol_disp_global = zeros( ...
    sStrain.vol_size(1),...
    sStrain.vol_size(2),...
    sStrain.vol_size(3),...
    4);

% main loop for fcc grains
for a0 = 1:5
    % Determine which sites are in this sector
    inds = sStrain.lat_ref_inds(a0,:);
    basis_test = sStrain.basis;
    basis_test(:,1) = 0;
    basis_test(:,inds) = 0;
    keep = sStrain.atoms_indexed & sum(abs(basis_test),2) == 0;

    % Get sites
    xyz = sStrain.atoms(keep,1:3);
    basis = sStrain.basis(keep,[1 inds]);
    
    % Fit mean lattice for this segment
    sStrain.lat_fcc_fit(:,:,a0) = basis \ xyz;
    % update origin for fcc ref lattice
    sStrain.lat_ref_fcc(1,:,a0) = sStrain.lat_fcc_fit(1,:,a0);
    
     % Reference lattice
    lat_ref = sStrain.lat_ref_fcc(:,:,a0);
        
    % Calculate displacements
    xyz_fit = basis * lat_ref;
    dxyz = xyz - xyz_fit;
    
    % volume site locations
    x = xyz(:,1) / sStrain.voxel_size + 0.5 + sStrain.vol_size(1)/2;
    y = xyz(:,2) / sStrain.voxel_size + 0.5 + sStrain.vol_size(2)/2;
    z = xyz(:,3) / sStrain.voxel_size + 0.5 + sStrain.vol_size(3)/2;
    x = min(max(x,1),sStrain.vol_size(1)-1);
    y = min(max(y,2),sStrain.vol_size(2)-1);
    z = min(max(z,3),sStrain.vol_size(3)-1);
    
    % trilinear interpolation
    xF = floor(x);
    yF = floor(y);
    zF = floor(z);
    dx = x - xF;
    dy = y - yF;
    dz = z - zF;
    inds3D = [ ...
        xF   yF   zF  ;
        xF+1 yF   zF  ;
        xF   yF+1 zF  ;
        xF+1 yF+1 zF  ;
        xF   yF   zF+1;
        xF+1 yF   zF+1;
        xF   yF+1 zF+1;
        xF+1 yF+1 zF+1];
    weights = [ ...
        (1-dx).*(1-dy).*(1-dz);
        (  dx).*(1-dy).*(1-dz);
        (1-dx).*(  dy).*(1-dz);
        (  dx).*(  dy).*(1-dz);
        (1-dx).*(1-dy).*(  dz);
        (  dx).*(1-dy).*(  dz);
        (1-dx).*(  dy).*(  dz);
        (  dx).*(  dy).*(  dz);
        ];
    inds = sub2ind( ...
        sStrain.vol_size, ...
        inds3D(:,1),...
        inds3D(:,2),...
        inds3D(:,3));
    
    % Local density
    vol_count(:) = 0;
    vol_count(inds) = weights;
    vol_count(:) = ifftn(fftn(vol_count).*k, 'symmetric');
    sub = vol_count > (max(vol_count(:))*0.01);
    vol_counts_sub = vol_count(sub);
    sStrain.vol_disp(:,:,:,4,a0) = vol_count;
    
    % xyz displacements
    for a1 = 1:3
        vol(:) = 0;
        vol(inds) = weights .* repmat(dxyz(:,a1),[8 1]);
        vol(:) = ifftn(fftn(vol).*k, 'symmetric');
        vol(sub) = vol(sub) ./ vol_counts_sub;
        vol(~sub) = 0;
        sStrain.vol_disp(:,:,:,a1,a0) = vol;
    end
end




% Calculate global displacements
keep = sStrain.atoms_indexed ;

% Get sites
xyz = sStrain.atoms(keep,1:3);
basis = sStrain.basis(keep,:);

% Fit mean reference lattice (just in case)
lat_ref = basis \ xyz;

% Calculate displacements
xyz_fit = basis * lat_ref;
dxyz = xyz - xyz_fit;

% volume site locations
x = xyz(:,1) / sStrain.voxel_size + 0.5 + sStrain.vol_size(1)/2;
y = xyz(:,2) / sStrain.voxel_size + 0.5 + sStrain.vol_size(2)/2;
z = xyz(:,3) / sStrain.voxel_size + 0.5 + sStrain.vol_size(3)/2;
x = min(max(x,1),sStrain.vol_size(1)-1);
y = min(max(y,2),sStrain.vol_size(2)-1);
z = min(max(z,3),sStrain.vol_size(3)-1);

% trilinear interpolation
xF = floor(x);
yF = floor(y);
zF = floor(z);
dx = x - xF;
dy = y - yF;
dz = z - zF;
inds3D = [ ...
    xF   yF   zF  ;
    xF+1 yF   zF  ;
    xF   yF+1 zF  ;
    xF+1 yF+1 zF  ;
    xF   yF   zF+1;
    xF+1 yF   zF+1;
    xF   yF+1 zF+1;
    xF+1 yF+1 zF+1];
weights = [ ...
    (1-dx).*(1-dy).*(1-dz);
    (  dx).*(1-dy).*(1-dz);
    (1-dx).*(  dy).*(1-dz);
    (  dx).*(  dy).*(1-dz);
    (1-dx).*(1-dy).*(  dz);
    (  dx).*(1-dy).*(  dz);
    (1-dx).*(  dy).*(  dz);
    (  dx).*(  dy).*(  dz);
    ];
inds = sub2ind( ...
    sStrain.vol_size, ...
    inds3D(:,1),...
    inds3D(:,2),...
    inds3D(:,3));

% Local density
vol_count(:) = 0;
vol_count(inds) = weights;
vol_count(:) = ifftn(fftn(vol_count).*k, 'symmetric');
sub = vol_count > (max(vol_count(:))*0.01);
vol_counts_sub = vol_count(sub);
sStrain.vol_disp_global(:,:,:,4) = vol_count;

% xyz displacements
for a1 = 1:3
    vol(:) = 0;
    vol(inds) = weights .* repmat(dxyz(:,a1),[8 1]);
    vol(:) = ifftn(fftn(vol).*k, 'symmetric');
    vol(sub) = vol(sub) ./ vol_counts_sub;
    vol(~sub) = 0;
    sStrain.vol_disp_global(:,:,:,a1) = vol;
end


% quick plotting
figure(111)
clf

% size(vol_count)
% imagesc(squeeze(sum(vol(:,:,71),3)))
% imagesc(squeeze(sum(vol_count,2)))
z = 71;
imagesc([ ...
    sum(sStrain.vol_disp(:,:,z,1,:),5) ...
    sum(sStrain.vol_disp(:,:,z,2,:),5) ...
    sum(sStrain.vol_disp(:,:,z,3,:),5) ...
    ])
% imagesc([ ...
%     sStrain.vol_disp(:,:,z,1,1) ...
%     sStrain.vol_disp(:,:,z,2,1) ...
%     sStrain.vol_disp(:,:,z,3,1) ...
%     ])
axis equal off
colormap(gray(256))
set(gca,'position',[0 0 1 1])

% hold on
% 
% scatter3( ...
%     xyz(:,2),...
%     xyz(:,1),...
%     xyz(:,3),...
%     'marker','o',...
%     'sizedata',10,...
%     'markeredgecolor',[0 0 0],...
%     'markerfacecolor',[1 1 1]);
% 
% 
% scatter3( ...
%     xyz_fit(:,2),...
%     xyz_fit(:,1),...
%     xyz_fit(:,3),...
%     'marker','s',...
%     'sizedata',100,...
%     'linewidth',0.5,...
%     'markeredgecolor',[1 0 0],...
%     'markerfacecolor','none');
% 
% hold off
% 
% axis equal
% box on
% 
% view([0 0 1])

toc
end