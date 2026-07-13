%% isInsideRectangle Function
% This function checks whether a set of points lies inside a rectangle defined by two opposite corners, X1 and X2.
% The rectangle is assumed to be axis-aligned.
% 
%% Inputs
% * P* A matrix of size Nx2, where each row represents a point (x, y)..
% * X1: A 1x2 vector representing the coordinates of one corner of the rectangle.
% * X2: A 1x2 vector representing the coordinates of opposite corner of the rectangle.
% 
%% Outputs
% * d: A column vector of size Nx1, where each element is 1 if the corresponding point in P is inside the rectangle, and 0 otherwise.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Function definition
function d = isInsideRectangle(P,X1,X2)
d = (P(:,1)>=X1(1)) .* (P(:,2)>=X1(2)) .* (P(:,1)<=X2(1)) .* (P(:,2)<=X2(2));
