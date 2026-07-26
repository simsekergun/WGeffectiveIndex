function bnd = boundary_nodes(elems, Nn)
% Find nodes on the boundary (appear in only one triangle shared edge).
    all_edges = [elems(:,[1,2]); elems(:,[2,3]); elems(:,[3,1])];
    sorted    = sort(all_edges, 2);
    [uniq, ~, ic] = unique(sorted, 'rows');
    cnt = accumarray(ic, 1);
    bnd_edges = uniq(cnt == 1, :);
    bnd = unique(bnd_edges(:));
end


