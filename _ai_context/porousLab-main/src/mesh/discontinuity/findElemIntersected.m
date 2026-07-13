%% findElemIntersected Function
% This function checks which bars (defined by their endpoints in the 
% NODE array) intersect a given rectangle defined by its minimum and 
% maximum coordinates. The algorithm uses parametric line equations and 
% checks for overlap in the parameter range for each bar.
%
%% Inputs
% * *Amin*: A vector specifying the minimum coordinates of the rectangle.
% * *Amax*: A vector specifying the maximum coordinates of the rectangle.
% * *NODE*: An matrix where each row represents the coordinates of a node.
% * *BARS*: An matrix where each row contains the indices of two nodes 
%           defining a bar.
%
%% Outputs
% * *elemID*: A logical array where each element is true if the 
%             corresponding bar intersects the rectangle, and false 
%             otherwise.
%
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function elemID = findElemIntersected(Amin,Amax,NODE,BARS)

% Amin and Amax are the rectangle's limit coords: minimum and maximum
Nb= size(BARS,1);
Tmin = zeros(Nb,1); Tmax = ones(Nb,1);
D = NODE(BARS(:,2),:) - NODE(BARS(:,1),:);
for i=1:2 % Check on X (i=1) and Y (i=2)
    T1 = ( Amin(i) - NODE(BARS(:,1),i) ) ./ D(:,i);
    T2 = ( Amax(i) - NODE(BARS(:,1),i) ) ./ D(:,i);
    ind = find(T1>T2); % We require T1<T2, swap if not
    [T1(ind),T2(ind)] = deal(T2(ind),T1(ind)); % Swap operation
    Tmin = max(Tmin,T1); Tmax = min(Tmax,T2);
end
% No intersection with rectangle if Tmin>=Tmax
elemID = (Tmin < Tmax)';