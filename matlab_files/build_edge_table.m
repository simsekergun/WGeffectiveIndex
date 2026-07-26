function [edges, elem2edge, edge_sign] = build_edge_table(elems)
% Build global edge list and element-to-edge connectivity with orientation.
    Ne   = size(elems,1);
    % All directed half-edges (element, local index)
    local_pairs = [2 3; 3 1; 1 2];
    all_edges   = zeros(3*Ne, 2);
    for ie = 1:Ne
        for k = 1:3
            all_edges(3*(ie-1)+k,:) = elems(ie, local_pairs(k,:));
        end
    end
    % Canonical edge = sorted pair
    sorted_edges = sort(all_edges, 2);
    [edges, ~, ic] = unique(sorted_edges, 'rows');
    Nedge = size(edges,1);

    elem2edge = reshape(ic, 3, Ne)';   % [Ne x 3]
    edge_sign = zeros(Ne,3);
    for ie = 1:Ne
        for k = 1:3
            raw = all_edges(3*(ie-1)+k,:);
            if raw(1) < raw(2)
                edge_sign(ie,k) = 1;
            else
                edge_sign(ie,k) = -1;
            end
        end
    end
end


