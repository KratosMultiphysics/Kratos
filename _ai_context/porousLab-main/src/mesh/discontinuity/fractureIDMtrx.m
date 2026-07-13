%% fractureIDMtrx Function
% This function identifies the elements in a mesh that are intersected by 
% fractures and computes an enrichment matrix to indicate these 
% intersections. It also computes an adjacency matrix to represent 
% connectivity between elements that share fracture nodes. The function 
% uses geometric intersection tests to determine whether a fracture 
% segment intersects an element edge.
% 
%% Inputs
% * *NODE*: Matrix containing the coordinates of the mesh nodes. Each row 
%           represents a node, with columns for x and y coordinates.
% * *ELEM*: Matrix defining the connectivity of the mesh elements. Each 
%           row represents an element, with columns indicating the node 
%           indices.
% * *NODE_D*: Matrix containing the coordinates of the nodes defining the 
%             fractures. Each row represents a fracture node, with columns 
%             for x and y coordinates.
% * *FRACT*: Matrix defining the fractures in the mesh. Each row 
%            represents a fracture segment, with columns indicating the 
%            indices of the nodes that define the segment.
%
%% Outputs
% * *IDenr*: Matrix indicating which elements are enriched by which 
%            fractures. A value of 1 indicates that the corresponding 
%            element is crossed by the corresponding fracture.
% * *IDAdj*: Sparse adjacency matrix indicating connectivity between 
%            elements based on shared fracture nodes. A value of 1 
%            indicates adjacency.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function [IDenr, IDAdj] = fractureIDMtrx(NODE,ELEM,NODE_D,FRACT)

% Initialize the matrix
IDenr = zeros(size(ELEM,1),size(FRACT,1)); 

nFract  = size(FRACT,1);
nElem   = size(ELEM,1);
counter = 0;

PINT = [];

IDFracNd = zeros(size(NODE_D,1),size(ELEM,1));

% Loop through the fracture segments --------------------------------------
for i = 1:nFract

    % Points that the define the fracture segment
    ndFrac1 = FRACT(i,1);
    ndFrac2 = FRACT(i,2);
    p1 = NODE_D(ndFrac1,:);
    p2 = NODE_D(ndFrac2,:);

    % Loop through the continuum elements ---------------------------------
    for el = 1:nElem 

        % Get the number of edges of the element
        nEdges = size(ELEM,2);

        % Get the coordinates of the element (repeat the first one to close
        % the polygon)
        cX = [NODE(ELEM(el,:),1); NODE(ELEM(el,1),1)];
        cY = [NODE(ELEM(el,:),2); NODE(ELEM(el,1),2)];

        % Loop through the edges of the element ---------------------------
        for j = 1:nEdges

            % Points that defined the element edge
            p3 = [cX(j)  , cY(j)];
            p4 = [cX(j+1), cY(j+1)];

            % Evaluate if the segments p1-p2 and p3-p4 intersect
            [flagInt,pint] = intersectionSegment([p1; p2],[p3; p4]);

            % Update the intersection vector points
            if flagInt == true
                PINT = [PINT;pint];
                counter = counter + 1;

                if norm(p1-pint)<1.0e-10
                    IDFracNd(ndFrac1,el) = 1;
                end
                if norm(p2-pint)<1.0e-10
                    IDFracNd(ndFrac2,el) = 1;
                end
            end

        end
        
        % If the fracture is crossing the element domain, it will have two
        % intersection points
        if size(unique(PINT,'rows')) == 2
            IDenr(el,i) = 1;
        end
        PINT = [];
        counter = 0;
    end
end

% Adjacency matrix
IDAdj = sparse(size(ELEM,1),size(ELEM,1));
for el = 1:nElem
    fracNds = IDFracNd(:,el);
    nFracNds = sum(fracNds);
    if nFracNds > 0
        id = find(fracNds);
        elAdj = setdiff(find(sum(IDFracNd(id,:),1)),el);
        IDAdj(el,elAdj) = 1;
    end
end


end