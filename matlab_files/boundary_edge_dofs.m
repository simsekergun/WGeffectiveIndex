function bnd_dofs = boundary_edge_dofs(edges, elems, Nn)
% DOF indices of edges on the mesh boundary (for PEC condition).
    all_e = [elems(:,[1,2]); elems(:,[2,3]); elems(:,[3,1])];
    sorted_e = sort(all_e, 2);
    [~, ~, ic] = unique(sorted_e, 'rows');
    cnt = accumarray(ic, 1);
    [~, ia] = unique(sort(edges,2), 'rows');
    bnd_edge_mask = (cnt(ia) == 1);
    bnd_dofs = find(bnd_edge_mask)';
end


