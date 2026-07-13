%% convertToQuadraticMesh function
% This function converts a linear mesh into a quadratic mesh.
% 
%% Inputs
% * *NODE*: A matrix of size Nx2 containing the x and y coordinates of 
%           the nodes. Each row represents a node, with the first column 
%           as the x-coordinate and the second column as the y-coordinate.
% * *ELEM*: A connectivity matrix of size Mx4 (for Q4 elements) or Mx3 
%           (for T3 elements). Each row represents an element, with 
%           columns specifying the indices of the nodes that form the 
%           element.
% 
%% Outputs
% * *NODE_quad*: An updated NODE matrix that includes the original nodes 
%                and the newly created mid-edge nodes.
% * *ELEM_quad*: An updated connectivity matrix that includes the 
%                quadratic connectivity for each element.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function [NODE_quad, ELEM_quad] = convertToQuadraticMesh(NODE, ELEM)
    % Initialize
    num_nodes = size(NODE, 1);
    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    new_nodes = [];

    num_elems = numel(ELEM);
    ELEM_quad = cell(num_elems, 1);  % Initialize cell for quadratic elements

    for elem_idx = 1:num_elems
        elem = ELEM{elem_idx}(:)';  % Ensure row vector
        num_nodes_per_elem = numel(elem);

        % Determine edge pairs based on element type
        if num_nodes_per_elem == 3
            % Triangular element (T3 -> T6)
            edge_pairs = [1, 2; 2, 3; 3, 1];
        elseif num_nodes_per_elem == 4
            % Quadrilateral element (Q4 -> Q8)
            edge_pairs = [1, 2; 2, 3; 3, 4; 4, 1];
        else
            error('Unsupported element type. ELEM should have 3 or 4 nodes.');
        end

        new_elem = elem;

        for edge_idx = 1:size(edge_pairs, 1)
            n1 = elem(edge_pairs(edge_idx, 1));
            n2 = elem(edge_pairs(edge_idx, 2));

            edge_key = sprintf('%d-%d', min(n1, n2), max(n1, n2));

            if isKey(edge_map, edge_key)
                mid_node = edge_map(edge_key);
            else
                mid_coord = (NODE(n1, :) + NODE(n2, :)) / 2;
                new_nodes = [new_nodes; mid_coord];
                mid_node = num_nodes + size(new_nodes, 1);
                edge_map(edge_key) = mid_node;
            end

            new_elem(end + 1) = mid_node;
        end

        ELEM_quad{elem_idx} = new_elem;
    end

    NODE_quad = [NODE; new_nodes];
end
