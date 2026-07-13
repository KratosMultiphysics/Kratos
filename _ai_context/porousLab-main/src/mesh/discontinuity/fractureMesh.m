%% fractureMesh Function
% This function creates a fracture mesh by dividing a fracture segment 
% based on the continuum finite element mesh. It does not assume that the 
% continuum mesh is structured. The function also perturbs nodes to avoid 
% coinciding nodes between the fracture and continuum meshes.
%
%% Inputs
% * *NODE*: Nodal coordinates of the continuum mesh (Nx2 array).
% * *ELEM*: Element connectivity of the continuum mesh (MxK array).
% * *XD*: Coordinates of the fracture segment endpoints (2x2 array).
% * *SEGD*: Indices of the fracture segment endpoints in XD (1x2 array).
% * *aperture*: Aperture of the fracture.
% * *leakoff*: Leakoff coefficient for the fracture.
% * *TIP*: Logical flag indicating if the fracture has a tip.
% * *FixedPressureJump*: Logical flag to fix pressure jump (default: true).
% * *FixedPf*: Logical flag to fix pressure field (default: false).
% * *FixedDisplJump*: Logical flag to fix displacement jump (default: true).
% * *ptol*: Perturbation tolerance for node adjustments (default: 1.0e-2).
%
%% Outputs
% * *NODE_D*: Nodal coordinates of the fracture mesh.
% * *FRACT*: Element connectivity of the fracture mesh.
% * *NODE_D_TIPS*: Coordinates of the fracture tip nodes.
% * *NODE*: Updated nodal coordinates of the continuum mesh.
% * *W*: Fracture aperture values.
% * *LEAKOFF*: Leakoff values for the fracture.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function [NODE_D, FRACT, NODE_D_TIPS, NODE, W, LEAKOFF] = fractureMesh(NODE, ELEM, XD, ...
    SEGD,aperture,leakoff, TIP, FixedPressureJump, FixedPf, FixedDisplJump, ptol)

if nargin < 10
    FixedPf           = false;
    FixedDisplJump    = true;
    FixedPressureJump = true; 
    ptol              = 1.0e-2;
end
if nargin < 11
    ptol              = 1.0e-2;
end

% -------------------------------------------------------------------------
% Attract node to the tip
% for i = 1:size(XD,1)
%     DIST = sqrt((NODE(:,1) - XD(i,1)).*(NODE(:,1) - XD(i,1)) + (NODE(:,2) - XD(i,2)).*(NODE(:,2) - XD(i,2)));
%     [minDist,id] = sort(DIST);
%     NODE(id(1),:) = XD(i,:);
% end

% -------------------------------------------------------------------------
% First trial 
[NODE_D, FRACT, NODE_D_TIPS, W, LEAKOFF] = fractureNodesAndElements(NODE, ELEM, XD, ...
    SEGD,aperture,leakoff, TIP, FixedPressureJump, FixedPf, FixedDisplJump);

% -------------------------------------------------------------------------
% Compute the mean characteristic lenght of the mesh

% Compute the mean characteristic length of each element in the mesh
Lce=zeros(size(ELEM,1),1);
for el = 1:size(ELEM,1)

    % Vertices of the element el coordinates
    vx = NODE(ELEM(el,:),1); 
    vy = NODE(ELEM(el,:),2);

    % Number of vertices 
    nv = length(ELEM(el,:)); 

    % Shifted vertices
    vxS = vx([2:nv 1]);
    vyS = vy([2:nv 1]); 

    % Compute the area of the element (trapezoidal rule)
    temp = vx.*vyS - vy.*vxS;
    Ae   = 0.5*sum(temp);
    
    % Characteristic lenght (quadrilateral elements)
    Lce(el) = sqrt(Ae);
end

% Compute the mean characteristic length of the elements associated with
% each node
lcmNode = zeros(size(NODE,1),1);

for i = 1:size(NODE,1)

    % Get the elements associated with this node
    idElem = any(ELEM == i, 2);

    % Compute the mean characteristic lenght of these nodes
    lcmNode(i) = mean(Lce(idElem));

end

% -------------------------------------------------------------------------
% Check if there is a node of the fracture mesh that coincides with the 
% nodes from the continuum mesh
LCMNode = repmat(lcmNode',size(NODE_D,1),1);
NDIST = zeros(size(NODE_D,1),size(NODE,1));
for i = 1:size(NODE_D,1)
    for j = 1:size(NODE,1)
        NDIST(i,j) = norm(NODE_D(i,:) - NODE(j,:));
    end
end

% Get the from the continuum mesh that are the same from the fracture mesh
[~,idNode] = find(NDIST < ptol*LCMNode);
if isempty(idNode), return, end

% -------------------------------------------------------------------------
% Borders and corners id
% * Valid only for rectangular domains

% Get the id of the nodes on each border of the domain
xmin = min(NODE(:,1)); idLeft   = find(abs(NODE(:,1) - xmin) < 1.0e-10);
xmax = max(NODE(:,1)); idRight  = find(abs(NODE(:,1) - xmax) < 1.0e-10);
ymin = min(NODE(:,2)); idBottom = find(abs(NODE(:,2) - ymin) < 1.0e-10);
ymax = max(NODE(:,2)); idTop    = find(abs(NODE(:,2) - ymax) < 1.0e-10);

% Get the id of the corner nodes
idLeftBottom  = intersect(idLeft,  idBottom);
idLeftTop     = intersect(idLeft,  idTop);
idRightBottom = intersect(idRight, idBottom);
idRightTop    = intersect(idRight, idTop);

% Remove the corner nodes from the borders array
idLeft   = setdiff(idLeft,   [idLeftBottom,   idLeftTop   ]);
idRight  = setdiff(idRight,  [idRightBottom, idRightTop   ]);
idBottom = setdiff(idBottom, [idLeftBottom,  idRightBottom]);
idTop    = setdiff(idTop,    [idRightTop,    idLeftTop    ]);

% Store array with the two sets
idCorners = [idLeftBottom, idLeftTop, idRightBottom, idRightTop];
idBorders = [idLeft; idRight; idBottom; idTop];

% Get the normal and tangential vectors of each border (RECTANGULAR DOMAIN)
nleft   = [-1.0, 0.0]; mleft   = [ 0.0,-1.0];
nright  = [ 1.0, 0.0]; mright  = [ 0.0, 1.0];
ntop    = [ 0.0, 1.0]; mtop    = [-1.0, 0.0];
nbottom = [ 0.0,-1.0]; mbottom = [ 1.0, 0.0];

% Check if there is a node to the perturbed at on of the corners
idCorners     = intersect(idCorners,idNode);
idLeftBottom  = intersect(idLeftBottom,  idNode);
idLeftTop     = intersect(idLeftTop,  idNode);
idRightBottom = intersect(idRightBottom, idNode);
idRightTop    = intersect(idRightTop, idNode);

% Check if there is a node to the perturbed at on of the borders
idBorders = intersect(idBorders,idNode);
idLeft    = intersect(idLeft,idNode);
idRight   = intersect(idRight,idNode);
idBottom  = intersect(idBottom,idNode);
idTop     = intersect(idTop,idNode);

% Get the id the center nodes
idCenter = setdiff(idNode,[idCorners,idBorders]);

% -------------------------------------------------------------------------
% Apply a perturbation to the center nodes

% Compute the normal vector of the discontinuity
% Defined considering nf = ez x mf, where ez = [0 0 1]
mf = (XD(SEGD(2),:) - XD(SEGD(1),:))/ norm(XD(SEGD(2),:) - XD(SEGD(1),:));
nf = [-mf(2) , mf(1)];

% Discontinuity center
xDCenter = (XD(SEGD(2),:) + XD(SEGD(1),:))/2.0;

if isempty(idCenter) == false
    for i = 1:length(idCenter)
        
        % Heaviside function
        dX = NODE(idCenter(i),:) - xDCenter;

        % Perturbation sign
        pertSign = sign(nf * dX');
        if pertSign >= 0, pertSign = 1; end

        NODE(idCenter(i),:) = NODE(idCenter(i),:) + ptol * lcmNode(idCenter(i)) * pertSign * nf;
    end
end

% -------------------------------------------------------------------------
% Apply a perturbation to the corner nodes
% ptol = 0.0;
if isempty(idCorners) == false
    if isempty(idLeftBottom) == false
        if dot(nf,nleft) > dot(nf,nbottom) % Projects on the left edge
            mlb = mleft;
        else % Projects on the bottom edge
            mlb = mbottom;
        end
        NODE(idLeftBottom,:) = NODE(idLeftBottom,:) + ptol * lcmNode(idLeftBottom) * mlb;
    end
    if isempty(idLeftTop) == false
        if dot(nf,nleft) > dot(nf,ntop) % Projects on the left edge
            mlt = mleft;
        else % Projects on the top edge
            mlt = mtop;
        end
        NODE(idLeftTop,:) = NODE(idLeftTop,:) + ptol * lcmNode(idLeftTop) * mlt;
    end
    if isempty(idRightBottom) == false
        if dot(nf,nright) > dot(nf,nbottom) % Projects on the right edge
            mrb = mright;
        else % Projects on the bottom edge
            mrb = mbottom;
        end
        NODE(idRightBottom,:) = NODE(idRightBottom,:) + ptol * lcmNode(idRightBottom) * mrb;
    end
    if isempty(idRightTop) == false
        if dot(nf,nright) > dot(nf,ntop) % Projects on the right edge
            mrt = mright;
        else % Projects on the top edge
            mrt = mtop;
        end
        NODE(idRightTop,:) = NODE(idRightTop,:) + ptol * lcmNode(idRightTop) * mrt;
    end
end


% -------------------------------------------------------------------------
% Apply a perturbation to the border nodes

if isempty(idBorders) == false
    if isempty(idLeft) == false
        NODE(idLeft,:) = NODE(idLeft,:) + ptol * lcmNode(idLeft) * mleft;
    end
    if isempty(idRight) == false
        NODE(idRight,:) = NODE(idRight,:) + ptol * lcmNode(idRight) * mright;
    end
    if isempty(idBottom) == false
        NODE(idBottom,:) = NODE(idBottom,:) + ptol * lcmNode(idBottom) * mbottom;
    end
    if isempty(idTop) == false
        NODE(idTop,:) = NODE(idTop,:) + ptol * lcmNode(idTop) * mtop;
    end
end

% -------------------------------------------------------------------------
% Regenerate the fracture mesh based on the perturbed mesh
[NODE_D, FRACT, NODE_D_TIPS, W, LEAKOFF] = fractureNodesAndElements(NODE, ELEM, XD, ...
    SEGD,aperture,leakoff, TIP, FixedPressureJump, FixedPf, FixedDisplJump);
end
