%% generateFract Function
% This function creates a set of nodes and connectivity information for a 
% fracture line defined between two points in 2D space. The fracture is 
% discretized into a specified number of divisions.
% 
%% Inputs
% * *X0*: A vector specifying the starting point of the fracture [x, y].
% * *X1*: A vector specifying the ending point of the fracture [x, y].
% * *nDiv*: An integer specifying the number of divisions for the fracture.
% 
%% Outputs
% * *NODE_D*: A matrix containing the coordinates of the nodes along the 
%             fracture. Each row represents a node, with the first column 
%             being the x-coordinate and the second column being the 
%             y-coordinate.
% * *FRACT*: A matrix representing the connectivity of the fracture 
%            segments. Each row defines a segment by specifying the 
%            indices of its two endpoints in NODE_D.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function [NODE_D,FRACT] = generateFract(X0, X1, nDiv)

NODE_D = [linspace(X0(1),X1(1),nDiv+1)',linspace(X0(2),X1(2),nDiv+1)'];

FRACT = zeros(nDiv,2);
FRACT(1,:) = [1 2];
for i = 2:nDiv
    FRACT(i,:) = FRACT(i-1,:) + [1 1];
end

end
