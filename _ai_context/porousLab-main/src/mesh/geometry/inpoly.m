function [cn, on] = inpoly ( p, node, edge, reltol )

%*****************************************************************************80
%
%  INPOLY: Point-in-polygon testing.
%
% Determine whether a series of points lie within the bounds of a polygon
% in the 2D plane. General non-convex, multiply-connected polygonal
% regions can be handled.
%
% SHORT SYNTAX:
%
%   in = inpoly(p, node);
%
%   p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
%   node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
%         The standard syntax assumes that the vertices are specified in
%         consecutive order.
%
%   in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
%         region.
%
% LONG SYNTAX:
%
%  [in, on] = inpoly(p, node, edge, tol);
%
%  edge: An Mx2 array of polygon edges, specified as connections between
%        the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
%        do not need to be specified in connsecutive order when using the
%        extended syntax.
%
%  on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
%        polygon edge. (A tolerance is used to deal with numerical
%        precision, so that points within a distance of
%        reltol*min(bbox(node)) from a polygon edge are considered "on" the 
%        edge.
%
% EXAMPLE:
%
%   polydemo;       % Will run a few examples
%
% See also INPOLYGON

% The algorithm is based on the crossing number test, which counts the
% number of times a line that extends from each point past the right-most
% region of the polygon intersects with a polygon edge. Points with odd
% counts are inside. A simple implementation of this method requires each
% wall intersection be checked for each point, resulting in an O(N*M)
% operation count.
%
% This implementation does better in 2 ways:
%
%   1. The test points are sorted by y-value and a binary search is used to
%      find the first point in the list that has a chance of intersecting
%      with a given wall. The sorted list is also used to determine when we
%      have reached the last point in the list that has a chance of
%      intersection. This means that in general only a small portion of
%      points are checked for each wall, rather than the whole set.
%
%   2. The intersection test is simplified by first checking against the
%      bounding box for a given wall segment. Checking against the bbox is
%      an inexpensive alternative to the full intersection test and allows
%      us to take a number of shortcuts, minimising the number of times the
%      full test needs to be done.
%
%  Author:
%
%    Darren Engwirda
%

%% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
   reltol = 1.0e-12;
   if nargin<3
      edge = [];
      if nargin<2
         error('Insufficient inputs');
      end
   end
end
nnode = size(node,1);
if isempty(edge)                                                           % Build edge if not passed
   edge = [(1:nnode-1)' (2:nnode)'; nnode 1];
end
if size(p,2)~=2
   error('P must be an Nx2 array.');
end
if size(node,2)~=2
   error('NODE must be an Mx2 array.');
end
if size(edge,2)~=2
   error('EDGE must be an Mx2 array.');
end
if max(edge(:))>nnode || any(edge(:)<1)
   error('Invalid EDGE.');
end

%% PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n  = size(p,1);
nc = size(edge,1);

% Choose the direction with the biggest range as the "y-coordinate" for the
% test. This should ensure that the sorting is done along the best
% direction for long and skinny problems wrt either the x or y axes.
dxy = max(p,[],1)-min(p,[],1);
if dxy(1)>dxy(2)
   % Flip co-ords if x range is bigger
   p = p(:,[2,1]);
   node = node(:,[2,1]);
end

% Polygon bounding-box
dxy = max(node,[],1)-min(node,[],1);
tol = reltol*min(dxy);
if tol==0.0
   tol = reltol;
end

% Sort test points by y-value
[y,i] = sort(p(:,2));
x = p(i,1);

%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cn = false(n,1);     % Because we're dealing with mod(cn,2) we don't have
                     % to actually increment the crossing number, we can
                     % just flip a logical at each intersection (faster!)
on = cn;
for k = 1:nc         % Loop through edges

   % Nodes in current edge
   n1 = edge(k,1);
   n2 = edge(k,2);

   % Endpoints - sorted so that [x1,y1] and [x2,y2] has y1<=y2
   %           - also get xmin = min(x1,x2), xmax = max(x1,x2)
   y1 = node(n1,2);
   y2 = node(n2,2);
   if y1<y2
      x1 = node(n1,1);
      x2 = node(n2,1);
   else
      yt = y1;
      y1 = y2;
      y2 = yt;
      x1 = node(n2,1);
      x2 = node(n1,1);
   end
   if x1>x2
      xmin = x2;
      xmax = x1;
   else
      xmin = x1;
      xmax = x2;
   end

   % Binary search to find first point with y<=y1 for current edge
   if y(1)>=y1
      start = 1;
   elseif y(n)<y1
      start = n+1;       
   else
      lower = 1;
      upper = n;
      for j = 1:n
         start = round(0.5*(lower+upper));
         if y(start)<y1
            lower = start+1;
         elseif y(start-1)<y1
            break;
         else
            upper = start-1;
         end
      end
   end

   % Loop through points
   for j = start:n
      % Check the bounding-box for the edge before doing the intersection
      % test. Take shortcuts wherever possible!

      Y = y(j);   % Do the array look-up once and make a temp scalar
      if Y<=y2
         X = x(j);   % Do the array look-up once and make a temp scalar
         if X>=xmin
            if X<=xmax

               % Check if we're "on" the edge
               on(j) = on(j) || (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<=tol);

               % Do the actual intersection test
               if (Y<y2) && ((y2-y1)*(X-x1)<(Y-y1)*(x2-x1))
                  cn(j) = ~cn(j);
               end

            end
         elseif Y<y2   % Deal with points exactly at vertices
            % Has to cross edge
            cn(j) = ~cn(j);
         end
      else
         % Due to the sorting, no points with >y
         % value need to be checked
         break
      end
   end

end

% Re-index to undo the sorting
cn(i) = cn|on;
on(i) = on;

end      % inpoly()
