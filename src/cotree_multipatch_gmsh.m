function [tree_dofs,cotree_dofs] = ...
  cotree_multipatch_gmsh (sp_h1, sp_hcurl, msh, drchlt_sides, ...
                          intrfc_dofs, geometry, isslave, ...
                          filename, exe_path)

if (sp_h1.npatch ~= sp_hcurl.npatch)
  error ('Spaces not defined on the same geometry')  
end

if nargin < 9
  exe_path = './topoprocessor';
  if nargin < 8
    filename = 'input.msh';
  end
end

ndim = sp_hcurl.ncomp;

% Compute the coordinates of the mapped Greville points
% And the numbering of the elements
greville_points = zeros (ndim, sp_h1.ndof);
nel_per_patch = zeros (sp_h1.npatch, 1);
for iptc = 1:sp_h1.npatch
  sph1_loc = sp_h1.sp_patch{iptc};
  pts = aveknt(sph1_loc.knots, sph1_loc.degree+1);
  greville_points(:,sp_h1.gnum{iptc}) = msh.msh_patch{iptc}.map(pts);
  
  ndof_dir_h1 = sph1_loc.ndof_dir;
  nelem_dir = ndof_dir_h1 - 1;
  nel_per_patch(iptc) = prod (nelem_dir);
end

nel = sum(nel_per_patch);
elements = zeros (2^ndim, nel);
inds = cell (ndim, 1);

prev_el = [0; cumsum(nel_per_patch)];
for iptc = 1:sp_h1.npatch
  sph1_loc = sp_h1.sp_patch{iptc};
  ndof_dir_h1 = sph1_loc.ndof_dir;
  nel_patch = nel_per_patch(iptc);
  nelem_dir = ndof_dir_h1 - 1;
  gnum = sp_h1.gnum{iptc};
  
  [inds{:}] = ind2sub (nelem_dir, 1:nel_patch);
  elem_indices = prev_el(iptc)+(1:nel_patch);
% Almost dimension independent
  elements(1,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}, inds{2}, inds{3}));
  elements(2,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}+1, inds{2}, inds{3}));
  elements(3,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}, inds{2}+1, inds{3}));
  elements(4,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}+1, inds{2}+1, inds{3}));
  elements(5,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}, inds{2}, inds{3}+1));
  elements(6,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}+1, inds{2}, inds{3}+1));
  elements(7,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}, inds{2}+1, inds{3}+1));
  elements(8,elem_indices) = gnum(sub2ind(ndof_dir_h1, inds{1}+1, inds{2}+1, inds{3}+1));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sp_loc = sp_hcurl.sp_patch{iptc};
  
  S_loc = zeros(sp_loc.ndof, 1);
  T_loc = zeros(sp_loc.ndof, 1);
  
  ndof_dir_h1 = sph1_loc.ndof_dir;
  cumsum_ndof = sp_loc.cumsum_ndof;
  for icomp = 1:sp_loc.ncomp_param
    inds = cell (sp_loc.ncomp_param, 1);
    ndof = sp_loc.scalar_spaces{icomp}.ndof;
    ndof_dir = sp_loc.scalar_spaces{icomp}.ndof_dir;
    [inds{:}] = ind2sub (ndof_dir, 1:ndof);
    S_loc(cumsum_ndof(icomp)+(1:ndof)) = sub2ind(ndof_dir_h1, inds{:});
    inds{icomp} = inds{icomp} + 1;
    T_loc(cumsum_ndof(icomp)+(1:ndof)) = sub2ind(ndof_dir_h1, inds{:});
  end
  
  % Without orientation
  S(sp_hcurl.gnum{iptc}) = gnum(S_loc);
  T(sp_hcurl.gnum{iptc}) = gnum(T_loc);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

switch_indices = [1 2 4 3 5 6 8 7];
elements = elements(switch_indices,:);

ptc_nsides = zeros(1,numel(msh.boundaries));
for ii=1:numel(msh.boundaries)
    ptc_nsides(ii) = msh.boundaries(ii).nsides;
end

Nbnd = cumsum ([0, ptc_nsides]);
bnd_dofs = [];
for iref = drchlt_sides
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    
  boundary_gnum = sp_hcurl.boundary.gnum;
  bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
end
drchlt_dofs = sp_hcurl.boundary.dofs(bnd_dofs);

weights = ones (size(S));
weights(drchlt_dofs) = 0;
%weights(intrfc_dofs) = 2;
if isslave
    weights(drchlt_dofs) = 0.1;
    weights(intrfc_dofs) = 0;
end
% if isslave == 2
%     weights(intrfc_dofs) = 0.1;
% end

% find corresponding nodes for edge dofs
nodes_intrfc = unique([S(intrfc_dofs), T(intrfc_dofs)]);

nodes_drchlt = unique([S(drchlt_dofs), T(drchlt_dofs)]);
nodes_drchlt = setdiff(nodes_drchlt, nodes_intrfc); % remove interface dofs from dirichlet dofss


G = graph(S,T,weights);
export_to_gmsh(filename, sp_h1.ndof, greville_points, nel, elements, nodes_intrfc);
system([exe_path, ' ', filename]); % calls the topoprocessor executable
g2_edges = load('tree.txt');

dofs = G.findedge(g2_edges(:,1),g2_edges(:,2));

% From the numbering in the graph to the numbering in GeoPDEs
space2graph = G.findedge(S,T);
[~,graph2space] = ismember (1:sp_hcurl.ndof, space2graph);

tree_dofs = sort (graph2space(dofs));
cotree_dofs = setdiff(1:sp_hcurl.ndof, tree_dofs);

if isslave == 1
   tree_dofs = setdiff(tree_dofs, intrfc_dofs);
   cotree_dofs = union(cotree_dofs, intrfc_dofs);
end

figure
plot_tree (geometry, sp_h1, sp_hcurl, tree_dofs);
end

function plot_tree (geometry, sp_h1, sp_hcurl, tree_dofs)

ndim = numel(sp_h1.sp_patch{1}.degree);
for iptc = 1:sp_h1.npatch
  sp_patch = sp_h1.sp_patch{iptc};
  for idim = 1:ndim
    grev_pts{idim} = aveknt(sp_patch.knots{idim}, sp_patch.degree(idim)+1);
  end
  mapped_pts = geometry(iptc).map(grev_pts);
  rdim = size(mapped_pts, 1);
  
  mapped_pts = reshape (mapped_pts, [rdim, sp_patch.ndof_dir]);

  for idim = 1:ndim
  for jdim = idim+1:ndim
    vv = arrayfun(@(x) 1:x(:), sp_patch.ndof_dir, 'UniformOutput', false);
    for idof = 1:sp_patch.ndof_dir(idim)
      vv{idim} = idof;
      for jdof = 1:sp_patch.ndof_dir(jdim)
        vv{jdim} = jdof;
        xpts = mapped_pts(:,vv{:});
        plot3 (xpts(1,:),xpts(2,:),xpts(3,:), 'color', 'k');
        hold on
      end
    end
  end
  end

  [~,~,tree_patch] = intersect (tree_dofs, sp_hcurl.gnum{iptc});
  sp_patch = sp_hcurl.sp_patch{iptc};
  indices = cell(ndim, 1);
  for icomp = 1:sp_patch.ncomp
    dofs = sp_patch.cumsum_ndof(icomp)+1:sp_patch.cumsum_ndof(icomp+1);
    [~,~,tree_comp] = intersect(tree_patch, dofs);
    [indices{:}] = ind2sub(sp_patch.scalar_spaces{icomp}.ndof_dir,tree_comp);

% Plot edges of the tree
    for ii = 1:numel(indices{1})
      vv = cellfun(@(x) x(ii), indices, 'UniformOutput', false);
      vv{icomp} = [indices{icomp}(ii), indices{icomp}(ii)+1];
      xpts = mapped_pts(:,vv{:});
      plot3(xpts(1,:),xpts(2,:),xpts(3,:), 'color', 'r', 'LineWidth', 2);
    end
  end
  
end

end
