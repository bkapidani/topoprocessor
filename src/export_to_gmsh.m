function export_to_gmsh (filename, ndof, vertices, nel, elements, ...
                         drchlt_vertices)

if nargin < 6
  drchlt_vertices = [];
  num_intrfc_dofs = 0;
else
  offset = size(elements,2);
  num_intrfc_dofs = size(drchlt_vertices,2);
  drchlt_vertices = [offset+1:offset+num_intrfc_dofs; 
                     15*ones(1,num_intrfc_dofs); ones(1,num_intrfc_dofs); 2*ones(1,num_intrfc_dofs);
                     drchlt_vertices];
end

vertices = [(1:ndof); vertices];

elements = [1:nel; 
            5*ones(1,nel); ones(2,nel);
            elements];

fid = fopen (filename, 'w');
fprintf(fid, '$MeshFormat\n');
fprintf(fid, '2.2 0 8\n');
fprintf(fid, '$EndMeshFormat\n');

fprintf(fid, '$Nodes\n');
fprintf(fid, '%i\n', ndof); % ndof nodes
fprintf(fid, '%i %3.12f %3.12f %3.12f\n', vertices);
fprintf(fid, '$EndNodes\n');

fprintf(fid, '$Elements\n');
fprintf(fid, '%i\n', nel+num_intrfc_dofs); % nel elements
% El. number, el. type (5=8-node hexahedron), unclear, unclear, vertices
fprintf(fid, '%i %i %i %i %i %i %i %i %i %i %i %i\n', elements);
if numel(drchlt_vertices)>0
  fprintf(fid, '%i %i %i %i %i\n', drchlt_vertices);
end
fprintf(fid, '$EndElements\n');
fclose(fid);

end
