function [sStrain] = strain01(pos,type,fivefoldCenter)
arguments
    pos(3,:) double
    type(1,:) double
    fivefoldCenter(1,3) double = [0,0,0]
end
% Colin Ophus - 2022 April 2
% modified by Zezhou Li - 2022 Oct
% This set of scripts calculates 3D strain tensors for pentagonal bipyramid
% (PB) and EPB structures.  Uses shared lattice vectors between twin planes.

% 01 - load the atomic positions and species, fit the lattices.  Run this
% script without an output to plot the initial lattice guess.
% Running this script with an output will use marching lattice vectors to
% assign a lattice sit ID to all atoms (if possible, within the tolerance).

% Inputs:
% pos - [3xN] array containing atoms positions in Angstroms.
% ID - [1xN] array containing species index.
report_progress = true;
plot_vectors = 8;  % number of vectors to plot for initial guess
sStrain.site_tolerance = 0.6875;  % tolerance for site assignment
% sStrain.site_tolerance = 0.5;  % tolerance for site assignment
% sStrain.site_tolerance = 2.75/6;  % tolerance for site assignment
sStrain.max_in_plane_inds = 2;  % how many in-plane non-zero indices allowed
% sStrain.num_fit_loops = 2;   % how many fitting iterations to run

sStrain.allow_negative_indices = false;

% Lattice initial guess:  [7x2] array
%     [or; u; v1; v2; v3; v4; v5]
%     or - origin
%     u  - vertical lattice vector

% PB from Jihan and Zezhou
% for PB
sStrain.lat = [ ...
    fivefoldCenter;
    0 0 2.75;
    [cos(pi*0/5) sin(pi*0/5)]*4.7829   0;
    [cos(pi*2/5) sin(pi*2/5)]*4.6873   0;
    [cos(pi*4/5) sin(pi*4/5)]*4.7527   0;
    [cos(pi*6/5) sin(pi*6/5)]*4.7277   0;
    [cos(pi*8/5) sin(pi*8/5)]*4.6856   0;
    ];

% for PdPt210605
% 拟合后不断更新
% sStrain.lat = [ ...
%     5.76 10 -0.7; % y x z，不是x y z
%     0 0 2.75;
%     [cos(pi*0/5) sin(pi*0/5)]*4.80   0;
%     [cos(pi*2/5) sin(pi*2/5)]*4.64   0;
%     [cos(pi*4/5) sin(pi*4/5)]*4.72   0;
%     [cos(pi*6/5) sin(pi*6/5)]*4.68   0;
%     [cos(pi*8/5) sin(pi*8/5)]*4.70   0;
%     ];
% sStrain.lat = [
%     5.8872   10.0050   -0.7812;
%     0.0113   -0.0018    2.7388;
%     4.4204    1.5385    0.0409;
%    -0.1641    4.8057    0.0000;
%    -4.5779    1.2421   -0.0497;
%    -2.6585   -3.9853   -0.0147;
%     2.8756   -3.8218    0.0511];
% sStrain.lat = [
%     5.8787   10.0185   -0.7803;
%     0.0107   -0.0007    2.7386;
%     4.4334    1.5449    0.0408;
%    -0.1727    4.7941    0.0006;
%    -4.5776    1.2555   -0.0487;
%    -2.6570   -3.9915   -0.0150;
%     2.8748   -3.8263    0.0511];


% Outputs:
% sStrain - struct containing all data, for example:
% sStrain.atoms - [Nx4] array containing atom data [x y z ID].
%               - Note that I use column vectors for atom data!

% init atoms
sStrain.atoms = [pos' type'];
sStrain.atoms_range = [ ...
    min(sStrain.atoms(:,1:3),[],1);
    max(sStrain.atoms(:,1:3),[],1);
    ];
sStrain.atoms_num = size(sStrain.atoms,1);


% Main script
if nargout == 0
    % plot the initial lattice guess
    lat = sStrain.lat;
    
    % [001] projection
    figure(101)
    clf
    hold on
    
    s = abs(sStrain.atoms(:,3) - lat(1,3)) < 1.5;
    scatter3( ...
        sStrain.atoms(s,2),...
        sStrain.atoms(s,1),...
        sStrain.atoms(s,3),...
        'marker','o',...
        'linewidth',0.5,...
        'sizedata',10,...
        'markerfacecolor',[1 1 1],...
        'markeredgecolor',[0 0 0])
    
    c = [0 0.7 1];
    a = (1:plot_vectors)';
    for a0 = 3:7
        line( ...
            [0 lat(a0,2)]*plot_vectors+lat(1,2),...
            [0 lat(a0,1)]*plot_vectors+lat(1,1),...
            [0 lat(a0,3)]*plot_vectors+lat(1,3),...
            'linewidth',2,...
            'color',c)
        p = lat(a0,:).*a;
        scatter3( ...
            p(:,2)+lat(1,2),...
            p(:,1)+lat(1,1),...
            p(:,3)+lat(1,3),...
            'marker','s',...
            'linewidth',1,...
            'sizedata',200,...
            'markerfacecolor','none',...
            'markeredgecolor',c)
    end
    hold off
    axis equal
    box on
    set(gca,'ydir','reverse')
    view([0 0 1])
    xlim(sStrain.atoms_range(:,2)+[-1; 1])
    ylim(sStrain.atoms_range(:,1)+[-1; 1])
    zlim(sStrain.atoms_range(:,3)+[-1; 1])
    
    
    
    % [100] projection
    figure(102)
    clf
    hold on
    
    s = abs(sStrain.atoms(:,2) - lat(1,2)) < 1;
    scatter3( ...
        sStrain.atoms(s,2),...
        sStrain.atoms(s,1),...
        sStrain.atoms(s,3),...
        'marker','o',...para.Property
        'linewidth',0.5,...
        'sizedata',10,...
        'markerfacecolor',[1 1 1],...
        'markeredgecolor',[0 0 0])
    
    c = [0 0.7 1;
        0 0.7 0];
    lat = sStrain.lat;
    
    scatter3( ...
        lat(1,2),...
        lat(1,1),...
        lat(1,3),...
        'marker','s',...
        'linewidth',1,...
        'sizedata',200,...
        'markerfacecolor','none',...
        'markeredgecolor',[0 0 0])
    
    a = (1:round(plot_vectors/2))';
    line( ...
        [0 1]*lat(2,2)*round(plot_vectors/2)+lat(1,2),...
        [0 1]*lat(2,1)*round(plot_vectors/2)+lat(1,1),...
        [0 1]*lat(2,3)*round(plot_vectors/2)+lat(1,3),...
        'linewidth',2,...
        'color',c(1,:))
    p = lat(2,:).*a;
    scatter3( ...
        p(:,2)+lat(1,2),...
        p(:,1)+lat(1,1),...
        p(:,3)+lat(1,3),...
        'marker','s',...
        'linewidth',1,...
        'sizedata',200,...
        'markerfacecolor','none',...
        'markeredgecolor',c(1,:))
    
    a = (-round(plot_vectors/2):-1)';
    line( ...
        [0 1]*lat(2,2)*-round(plot_vectors/2)+lat(1,2),...
        [0 1]*lat(2,1)*-round(plot_vectors/2)+lat(1,1),...
        [0 1]*lat(2,3)*-round(plot_vectors/2)+lat(1,3),...
        'linewidth',2,...
        'color',c(2,:))
    p = lat(2,:).*a;
    scatter3( ...
        p(:,2)+lat(1,2),...
        p(:,1)+lat(1,1),...
        p(:,3)+lat(1,3),...
        'marker','s',...
        'linewidth',1,...
        'sizedata',200,...
        'markerfacecolor','none',...
        'markeredgecolor',c(2,:))
    
    %     for a0 = 3:7
    %         line( ...
    %             [0 lat(a0,2)]*plot_vectors+lat(1,2),...
    %             [0 lat(a0,1)]*plot_vectors+lat(1,1),...
    %             [0 lat(a0,3)]*plot_vectors+lat(1,3),...
    %             'linewidth',2,...
    %             'color',c)
    %         p = lat(a0,:).*a;
    %         scatter3( ...
    %             p(:,2)+lat(1,2),...
    %             p(:,1)+lat(1,1),...
    %             p(:,3)+lat(1,3),...
    %             'marker','s',...
    %             'linewidth',1,...
    %             'sizedata',200,...
    %             'markerfacecolor','none',...
    %             'markeredgecolor',c)
    %     end
    hold off
    axis equal
    box on
    set(gca,'ydir','reverse')
    view([1 0 0])
    xlim(sStrain.atoms_range(:,2)+[-1; 1])
    ylim(sStrain.atoms_range(:,1)+[-1; 1])
    zlim(sStrain.atoms_range(:,3)+[-1; 1])
    
else
    % Fit the lattice indices for all atomic sites
    
    %%% Fitting 01 %%%%
    % First, only search for the cardinal lattice directions
    
    
    % basis vectors to check for NNs
    d_basis = [ ...
        0   1 0 0 0 0 0;
        %         0   0 1 0 0 0 0;
        %         0   0 0 1 0 0 0;
        %         0   0 0 0 1 0 0;
        %         0   0 0 0 0 1 0;
        %         0   0 0 0 0 0 1;
        ...
        0   -1 0 0 0 0 0;
        0   0 -1 0 0 0 0;
        0   0 0 -1 0 0 0;
        0   0 0 0 -1 0 0;
        0   0 0 0 0 -1 0;
        0   0 0 0 0 0 -1;
        ...
        0   0.5 0.5 0.0 0.0 0.0 0.0;
        0   0.5 0.0 0.5 0.0 0.0 0.0;
        0   0.5 0.0 0.0 0.5 0.0 0.0;
        0   0.5 0.0 0.0 0.0 0.5 0.0;
        0   0.5 0.0 0.0 0.0 0.0 0.5;
        ...
        0   -0.5 0.5 0.0 0.0 0.0 0.0;
        0   -0.5 0.0 0.5 0.0 0.0 0.0;
        0   -0.5 0.0 0.0 0.5 0.0 0.0;
        0   -0.5 0.0 0.0 0.0 0.5 0.0;
        0   -0.5 0.0 0.0 0.0 0.0 0.5;
        ...
        0   0.5 -0.5 0.0 0.0 0.0 0.0;
        0   0.5 0.0 -0.5 0.0 0.0 0.0;
        0   0.5 0.0 0.0 -0.5 0.0 0.0;
        0   0.5 0.0 0.0 0.0 -0.5 0.0;
        0   0.5 0.0 0.0 0.0 0.0 -0.5;
        ...
        0   -0.5 -0.5 0.0 0.0 0.0 0.0;
        0   -0.5 0.0 -0.5 0.0 0.0 0.0;
        0   -0.5 0.0 0.0 -0.5 0.0 0.0;
        0   -0.5 0.0 0.0 0.0 -0.5 0.0;
        0   -0.5 0.0 0.0 0.0 0.0 -0.5;
        ...
        0   0.0 0.5 0.5 0.0 0.0 0.0;
        0   0.0 0.0 0.5 0.5 0.0 0.0;
        0   0.0 0.0 0.0 0.5 0.5 0.0;
        0   0.0 0.0 0.0 0.0 0.5 0.5;
        0   0.0 0.5 0.0 0.0 0.0 0.5;
        ...
        0   0.0 -0.5 0.5 0.0 0.0 0.0;
        0   0.0 0.0 -0.5 0.5 0.0 0.0;
        0   0.0 0.0 0.0 -0.5 0.5 0.0;
        0   0.0 0.0 0.0 0.0 -0.5 0.5;
        0   0.0 -0.5 0.0 0.0 0.0 0.5;
        ...
        0   0.0 0.5 -0.5 0.0 0.0 0.0;
        0   0.0 0.0 0.5 -0.5 0.0 0.0;
        0   0.0 0.0 0.0 0.5 -0.5 0.0;
        0   0.0 0.0 0.0 0.0 0.5 -0.5;
        0   0.0 0.5 0.0 0.0 0.0 -0.5;
        ...
        0   0.0 -0.5 -0.5 0.0 0.0 0.0;
        0   0.0 0.0 -0.5 -0.5 0.0 0.0;
        0   0.0 0.0 0.0 -0.5 -0.5 0.0;
        0   0.0 0.0 0.0 0.0 -0.5 -0.5;
        0   0.0 -0.5 0.0 0.0 0.0 -0.5;
        ];
    
    % Loops
    r2max = sStrain.site_tolerance^2;
    comp = 0;
    %     for a0 = 1:sStrain.num_fit_loops
    dxyz = d_basis * sStrain.lat;
    
    % init
    sStrain.basis = zeros(sStrain.atoms_num,7);
    sStrain.basis(:,1) = 1;
    mark = false(size(sStrain.atoms,1),1);
    [~,ind] = min(sum((sStrain.atoms(:,1:3) - sStrain.lat(1,:)).^2,2));
    mark(ind) = true;
    test = [sStrain.atoms(ind,1:3) + dxyz ...
        sStrain.basis(ind,:) + d_basis];
    
    % March over all lattice sites
    count = 0;
    march = true;
    while march == true
        if report_progress == true
            % progress
            comp = sum(mark)/size(sStrain.atoms,1);
            progressbar(comp,2);
        end
        
        % Find all sites closest to test locations
        dist2_ind = zeros(size(test,1),2);
        for a1 = 1:size(test,1)
            dist2 = sum((sStrain.atoms(:,1:3) - test(a1,1:3)).^2,2);
            dist2(mark) = inf;  % remove sites already identified
            [dist2_ind(a1,1),dist2_ind(a1,2)] = min(dist2);
        end
        
        % Add new sites
        sub = dist2_ind(:,1) < r2max;
        inds_add = dist2_ind(sub,2);
        basis_add = test(sub,3+(1:7));
        mark(inds_add) = true;
        sStrain.basis(inds_add,:) = basis_add;
        
        % New candidate sites
        test_init = zeros(size(d_basis,1)*length(inds_add),10);
        for a1 = 1:length(inds_add)
            inds = (1:size(d_basis,1)) + (a1-1)*size(d_basis,1);
            test_init(inds,1:3) = sStrain.atoms(inds_add(a1),1:3) + dxyz;
            test_init(inds,4:10) = sStrain.basis(inds_add(a1),:) + d_basis;
        end
        
        % Remove duplicates, already existing entries, and any rows
        % with more than 4 non-zero basis components
        basis_unique = unique(test_init(:,4:10),'rows');
        already_indexed = ismember(basis_unique,sStrain.basis,'rows');
        basis_unique(already_indexed,:) = [];
        too_many_inds = sum(basis_unique(:,3:end)~=0,2) ...
            > sStrain.max_in_plane_inds;
        basis_unique(too_many_inds,:) = [];
        
        % remove negative indices if required
        if sStrain.allow_negative_indices == false
            inds_neg = any(basis_unique(:,3:end) < 0,2);
            basis_unique(inds_neg,:) = [];
        end
        
        % Final list of candidate peaks
        if size(basis_unique,1) > 0
            test = [zeros(size(basis_unique,1),3) basis_unique];
            for a1 = 1:size(basis_unique,1)
                sub = all(test_init(:,4:10)==basis_unique(a1,:),2);
                test(a1,1:3) = mean(test_init(sub,1:3),1);
            end
        else
            march = false;
        end
        
        count = count + 1;
        if count > 6
            % update lattice
            sStrain.lat = sStrain.basis(mark,:) \ sStrain.atoms(mark,1:3);
        end
    end
    
    % update lattice
    sStrain.lat = sStrain.basis(mark,:) \ sStrain.atoms(mark,1:3);
    
    % Report statistics
    disp([ ...
        num2str(sum(mark)) '/' ...
        num2str(size(sStrain.atoms,1)) ' sites indexed'])
    %     disp(['Iter ' num2str(a0) '/' ...
    %         num2str(sStrain.num_fit_loops) ', ' ...
    %         num2str(sum(mark)) '/' ...
    %         num2str(size(sStrain.atoms,1)) ' sites indexed'])
    %     end
    
    if report_progress == true && comp < 1
        progressbar(1,2);
    end
    
    % fit lattice positions
    sStrain.atoms_indexed = mark;
    sStrain.xyz_fit = sStrain.basis * sStrain.lat;
    
    
    % plot results
    figure(101)
    clf
    hold on
    
    
    scatter3( ...
        sStrain.atoms(~mark,2),...
        sStrain.atoms(~mark,1),...
        sStrain.atoms(~mark,3),...
        'marker','o',...
        'linewidth',0.5,...
        'sizedata',20,...
        'markerfacecolor',[0 0.7 1],...
        'markeredgecolor','none')
    scatter3( ...
        sStrain.atoms(mark,2),...
        sStrain.atoms(mark,1),...
        sStrain.atoms(mark,3),...
        'marker','o',...
        'linewidth',0.5,...
        'sizedata',20,...
        'markerfacecolor',[0 0 0],...
        'markeredgecolor','none')
    
    
    %     s = abs(sStrain.atoms(:,3) - sStrain.lat(1,3)) < 1.5;
    %     scatter3( ...
    %         sStrain.atoms(:,2),...
    %         sStrain.atoms(:,1),...
    %         sStrain.atoms(:,3),...
    %         'marker','o',...
    %         'linewidth',0.5,...
    %         'sizedata',5,...
    %         'markerfacecolor',[0 0 0],...
    %         'markeredgecolor','none')
    
    
    %     scatter3( ...
    %         sStrain.xyz_fit(:,2),...
    %         sStrain.xyz_fit(:,1),...
    %         sStrain.xyz_fit(:,3),...
    %         'marker','o',...
    %         'linewidth',0.5,...
    %         'sizedata',20,...
    %         'markerfacecolor','none',...
    %         'markeredgecolor',[0 1 0])
    
    v = [sStrain.atoms(mark,1:3);
        sStrain.xyz_fit(mark,1:3)];
    f = [(1:sum(mark))' (1:sum(mark))' + sum(mark)];
    patch('vertices',v(:,[2 1 3]),...
        'faces',f,...
        'edgecolor',[1 0 0],...
        'linewidth',1);
    
    
    
    
    hold off
    axis equal
    box on
    set(gca,'ydir','reverse')
    view([0 0 1])
    xlim(sStrain.atoms_range(:,2)+[-1; 1])
    ylim(sStrain.atoms_range(:,1)+[-1; 1])
    zlim(sStrain.atoms_range(:,3)+[-1; 1])
    
    
end




end

