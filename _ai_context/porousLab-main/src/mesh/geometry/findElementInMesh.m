%% findElementInMesh function
% This function determines the element in a 2D mesh that contains a given 
% point P. It supports both triangular and quadrilateral elements. The 
% function calculates the area of the element and compares it with the 
% sum of the areas formed by the point P and the vertices of the element. 
% If the areas match within a tolerance, the point is considered to be 
% inside the element.
% 
%% Inputs
% * *NODE*: A matrix containing the coordinates of the nodes. Each row 
%           represents a node, with the first column as the x-coordinate 
%           and the second column as the y-coordinate.
% * *ELEM*: A matrix containing the indices of the nodes that form each 
%           element. Each row represents an element, and the number of 
%           columns corresponds to the number of vertices.
% * *P*: A vector representing the coordinates of the point to locate 
%        in the mesh.
% 
%% Outputs
% * *elemId*: The index of the element in which the point P is located. If 
%             the point is not inside any element, elemId is set to 0.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function elemId = findElementInMesh(NODE, ELEM, P)
% Define in which convex quadrilateral element a point is inside
    nelem = size(ELEM, 1);

    for i = 1:nelem

        % Number of vertices
        nv = length(ELEM{i});

        % Get the corner nodes
        if (nv == 6) || (nv ==8)
            nv = nv / 2;
        end
        
        % Get the coordinates of the corner nodes
        vertices = NODE(ELEM{i}(1:nv), :); 
        
        % Define vertices
        if nv == 3
            pv1 = P - vertices(1, :);
            pv2 = P - vertices(2, :);
            pv3 = P - vertices(3, :);
        elseif nv == 4
            pv1 = P - vertices(1, :);
            pv2 = P - vertices(2, :);
            pv3 = P - vertices(3, :);
            pv4 = P - vertices(4, :);
        end

        % Calculate the cross product
        % Define vertices
        if nv == 3
            vA1 = cross([pv1, 0], [pv2, 0]);
            vA2 = cross([pv2, 0], [pv3, 0]);
            vA3 = cross([pv3, 0], [pv1, 0]);
        elseif nv == 4
            vA1 = cross([pv1, 0], [pv2, 0]);
            vA2 = cross([pv2, 0], [pv3, 0]);
            vA3 = cross([pv3, 0], [pv4, 0]);
            vA4 = cross([pv4, 0], [pv1, 0]);
        end
        
        
        % Calculate the area of each triangle
        if nv == 3
            A1 = 0.5 * norm(vA1);
            A2 = 0.5 * norm(vA2);
            A3 = 0.5 * norm(vA3);
            Ai = A1 + A2 + A3;
        elseif nv == 4
            A1 = 0.5 * norm(vA1);
            A2 = 0.5 * norm(vA2);
            A3 = 0.5 * norm(vA3);
            A4 = 0.5 * norm(vA4);
            Ai = A1 + A2 + A3 + A4;
        end
        
        % Area of the element
        Aelem = calculateQuadrilateralArea(vertices);
        
        if abs(Aelem - Ai ) < 1e-9
            elemId = i;
            return;
        end
    end

    % If P is outside all elements
    elemId = 0;

end

function area = calculateQuadrilateralArea(vertices)
    % Ensure the vertices are in a counterclockwise order
    if isClockwise(vertices)
        vertices = flipud(vertices);
    end
    
    % Calculate the area using the shoelace formula
    n = size(vertices, 1);
    area = 0.5 * abs(sum(vertices(1:end-1, 1) .* vertices(2:end, 2)) ...
                   + vertices(end, 1) * vertices(1, 2) ...
                   - sum(vertices(1:end-1, 2) .* vertices(2:end, 1)) ...
                   - vertices(end, 2) * vertices(1, 1));
end

function result = isClockwise(vertices)
    % Check if the vertices are listed in clockwise order
    n = size(vertices, 1);
    sum_angle = 0;
    for i = 1:n
        v1 = vertices(i, :);
        v2 = vertices(mod(i, n) + 1, :);
        sum_angle = sum_angle + atan2(v2(2) - v1(2), v2(1) - v1(1));
    end
    result = sum_angle < 0;
end
