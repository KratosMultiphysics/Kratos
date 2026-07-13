%% regularMesh
% Generate a 2D structured mesh on a rectangle with controllable nonuniform
% spacing and configurable element type.
%
% Features
% - ISOQ4 (4-noded quads) or CST (3-noded triangles; two per quad cell).
% - Optional fixed grid lines in x and y.
% - Monotone node redistribution in each direction via an inverted-Gaussian
% PDF, integrated to a CDF (“CIDF warp”) for smooth clustering.
% - Node numbering chooses the shorter side first to reduce bandwidth.
%
%% Inputs
% * Lx : Domain length in x (scalar > 0).
% * Ly : Domain length in y (scalar > 0).
% * Nx : Target number of cells in x (integer >= 1). Used to build the
% initial uniform parameterization and to set the warp width.
% * Ny : Target number of cells in y (integer >= 1). Same role as Nx.
% * xo : Vector of fixed x coordinates to force into the mesh
% (optional, default = []). Must lie in [0, Lx].
% * yo : Vector of fixed y coordinates to force into the mesh
% (optional, default = []). Must lie in [0, Ly].
% * type: Element type (char, optional, default = 'ISOQ4').
% Accepted: 'ISOQ4' or 'CST'.
% * cx : X-warp center in parametric space s∈[0,1] (optional, default = 0).
% * Ax : X-warp amplitude in [0,1] (optional, default = 0). Ax=0 ⇒ uniform.
% * cy : Y-warp center in s∈[0,1] (optional, default = 0).
% * Ay : Y-warp amplitude in [0,1] (optional, default = 0). Ay=0 ⇒ uniform.
%
% Warp model (per direction)
% - Start with uniform s=linspace(0,1,N+1).
% - PDF(s) = 1 − A * exp(−(s−c)^2 / (2 σ^2)), then normalize ∫PDF=1.
% - CDF(s) = ∫ PDF ds; physical coords = L * CDF(s).
% - Width σ is set from N: σ = clamp(k/N, 0, 1) with k≈8–10 in code.
% Larger A increases clustering away from c (inverted Gaussian “dip”).
%
%% Outputs
% * Node: (Nn × 2) array of coordinates with Nn = (Nx_eff+1)(Ny_eff+1).
% Nx_eff and Ny_eff are the post-merge counts after inserting
% fixed lines (may exceed input Nx, Ny). Columns are [x y].
% * ELEM: Cell array of element connectivities.
% - If 'ISOQ4': length(ELEM) = Nx_eff * Ny_eff, each cell has 4 node ids.
% - If 'CST' : length(ELEM) = 2 * Nx_eff * Ny_eff, each cell has 3 node ids.
%
%% Node numbering and connectivity
% * If Nx_eff < Ny_eff, nodes are numbered first along x then y; otherwise
% first along y then x. This choice targets lower matrix bandwidth.
% * Quad local ordering is consistent with the chosen numbering branch.
% * Triangulation splits each quad into two CSTs; the diagonal orientation
% follows the numbering branch used in the code.
%
%% Fixed lines handling
% * Fixed coords in xo or yo are merged with generated grids and sorted.
% * Duplicates within a tolerance tol = 0.2*L/N (per direction) are removed.
% * Inserting fixed lines increases Nx_eff or Ny_eff accordingly.
%
%% Defaults
% xo = [], yo = [], type = 'ISOQ4', cx = 0, Ax = 0, cy = 0, Ay = 0.
%
%% Notes
% - Inputs Nx, Ny drive the initial spacing and σ only; effective counts
% may grow after merging fixed lines.
% - All indices are 1-based (MATLAB). Coordinates are in the given units.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% v1.1 Updated docs for nonuniform warping, fixed-line merge, adaptive numbering.
% v1.0 Initial version (uniform spacing, ISOQ4/CST).

% -------------------------------------------------------------------------
% getUniquePoints
% Merge generated coordinates with user-fixed coordinates, drop near-duplicates,
% and return a sorted vector.
%
% Inputs
% * x0 : Candidate coordinates (row or column vector).
% * xfix: Fixed coordinates to enforce (vector; [] allowed).
% * tol : Minimum separation; points within tol of any xfix are skipped.
%
% Output
% * x : Sorted unique coordinates containing all xfix and the subset of
% x0 that are at least tol away from every xfix.
%% Function definition
function [Node,ELEM] = regularMesh(Lx,Ly,Nx,Ny,xo,yo,type,cx,Ax,cy,Ay)

if nargin < 5, xo = []; yo = []; end
if nargin < 7, type = 'ISOQ4'; end
if nargin < 9, cx = 0; Ax = 0; end
if nargin < 11, cy = 0; Ay = 0; end

% Uniform parameter in [0,1]
s = linspace(0,1,Nx+1);  

% width of the dip
sigma = max(min(8/Nx,1.0),0.0);

% 2) inverted‐Gaussian “dip” PDF
pdf = 1 - Ax*exp( -((s - cx).^2) ./ (2*sigma^2) );

% 3) normalize PDF so area=1 (optional but cleaner)
pdf = pdf / trapz(s,pdf);

% 4) CIDF warp
ux = cumtrapz(s, pdf);
ux = ux / ux(end);     % ensure u(1)=1 exactly
xcoord = getUniquePoints(Lx*ux, xo, 0.2*Lx/Nx);

% Uniform parameter in [0,1]
s = linspace(0,1,Ny+1);   

% width of the dip
sigma = max(min(10/Ny,1.0),0.0);

% 2) inverted‐Gaussian “dip” PDF
pdf = 1 - Ay*exp( -((s - cy).^2) ./ (2*sigma^2) );

% 3) normalize PDF so area=1 (optional but cleaner)
pdf = pdf / trapz(s,pdf);

% 4) CIDF warp
uy = cumtrapz(s, pdf);
uy = uy / uy(end);     % ensure u(1)=1 exactly
ycoord = getUniquePoints(Ly*uy, yo, 0.2*Ly/Ny);

% Number of nodes in each direction
Nx = length(xcoord) - 1;
Ny = length(ycoord) - 1;

% Coordinates of the nodes
if Nx < Ny
    % Number the nodes first along x and then along y
    [Y,X]= meshgrid(ycoord,xcoord);
else
    % Number the nodes first along y and then along x
    [X,Y]= meshgrid(xcoord,ycoord);
end
Node= [reshape(X,numel(X),1) reshape(Y,numel(Y),1)];

% Initialize the element matrix
if strcmp(type,'CST')
    ELEM= cell(2*Nx*Ny,1);
else
    ELEM= cell(Nx*Ny,1);
end

% Numbering elements first along x and then along y
k= 1;
if Nx < Ny 
    for j=1:Ny
        for i=1:Nx  
            n1 = (j-1)*(Nx+1)+i; n2 = j*(Nx+1)+i;
            if strcmp(type,'CST')
                ELEM{k}   = [n1, n2+1, n2];
                ELEM{k+1} = [n1, n1+1, n2+1];
                k = k+2;
            else
                ELEM{k} = [n1, n1+1, n2+1, n2];
                k = k+1;
            end
        end
    end
else
    for j=1:Ny
        for i=1:Nx  
            n1 = (i-1)*(Ny+1)+j; n2 = i*(Ny+1)+j;
            if strcmp(type,'CST')
                ELEM{k}   = [n1 n2 n2+1];
                ELEM{k+1} = [n2+1 n1+1 n1];
                k = k+2;
            else
                ELEM{k} = [n1 n2 n2+1 n1+1];
                k = k+1;
            end
        end
    end
end

end

% -------------------------------------------------------------------------
% Filters and combines unique points from two sets of points
function x = getUniquePoints(x0,xfix,tol)
if isempty(xfix) == true
    x = x0;
    return
end
x = xfix;
for i = 1:length(x0)
    dx = abs(x0(i) - xfix);
    if any(dx < tol) == false
        x = [x, x0(i)];
    end
end
x = sort(x);
end