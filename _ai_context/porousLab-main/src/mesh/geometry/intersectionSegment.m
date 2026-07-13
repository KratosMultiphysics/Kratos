%% intersectionSegment Function
% This function determines whether two line segments intersect and, if so,
% computes the intersection point and the parametric coordinates of the
% intersection along each segment.

%% Inputs
% * *segment1*: Endpoints of the first segment, specified as 1x2 vectors.
% * *segment2*: Endpoints of the second segment, specified as 1x2 vectors.
% 
%% Outputs
% * *flagInt*: Flag indicating whether the segments intersect or not.
% * *pint*: Intersection point as a 1x2 vector. Empty if no intersection.
% * *t12*: Parametric coordinate of the intersection along the first 
%          segment. Empty if no intersection.
% * *t34*: Parametric coordinate of the intersection along the second 
%          segment. Empty if no intersection.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function [flagInt,pint,t12,t34] = intersectionSegment(segment1,segment2)

    % Extract the points
    p1 = segment1(1, :);
    p2 = segment1(2, :);
    p3 = segment2(1, :);
    p4 = segment2(2, :);

    % Coordinates of the fracture segment
    x12_l = min(p1(1),p2(1));
    x12_r = max(p1(1),p2(1));
    y12_b = min(p1(2),p2(2));
    y12_t = max(p1(2),p2(2));

    % Coordinates of the element edge
    x34_l = min(p3(1),p4(1));
    x34_r = max(p3(1),p4(1));
    y34_b = min(p3(2),p4(2));
    y34_t = max(p3(2),p4(2));

    % Tolerance of the bounding box
    tol = 1e-9;

    % Discard intersection if the continuum edge is located to
    % the left or right of the horizontal bounding box of the first
    % segment
    if ((x12_r + tol) < x34_l) || (x34_r < (x12_l - tol))
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Discard intersection if the continuum edge is located bellow
    % or above of the bounding box of the first
    % segment
    if ((y12_t + tol) < y34_b) || (y34_t < (y12_b - tol))
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Get signs of oriented twice area for points p1-p2-p3 and for
    % points p1-p2-p4
    sign123 = signArea2d(p1, p2, p3);
    sign124 = signArea2d(p1, p2, p4);

    % Check if the segments are collinear
    if (sign123 == 0.0) && (sign124 == 0.0)
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Check if the second segment is above or on the right side of
    % the first segment
    if (sign123 > 0.0) && (sign124 > 0.0)
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Check if the second segment is bellow or on the left side of
    % the first segment
    if (sign123 < 0.0) && (sign124 < 0.0)
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Get signs of oriented twice area for points p1-p2-p3 and for
    % points p1-p2-p4
    sign341 = signArea2d(p3, p4, p1);
    sign342 = signArea2d(p3, p4, p2);

    % Check if the second segment is above or on the right side of
    % the first segment
    if (sign341 > 0.0) && (sign342 > 0.0)
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Check if the second segment is bellow or on the left side of
    % the first segment
    if (sign341 < 0.0) && (sign342 < 0.0)
        flagInt = false;
        pint    = [];
        t12     = [];
        t34     = [];
        return
    end

    % Check for one point of the second segment touching the first segment
    area341 =  area2d(p3, p4, p1);
    area342 =  area2d(p3, p4, p2);
    if sign123 == 0.0
        flagInt = true;
        t34  = 0.0;
        t12  = area341 / (area341 - area342);
        pint = p3;
        return
    end
    if sign124 == 0.0
        flagInt = true;
        t34  = 1.0;
        t12  = area341 / (area341 - area342);
        pint = p4;
        return
    end

    % Check for one point of the first segment touching the second segment
    area123 =  area2d(p1, p2, p3);
    area124 =  area2d(p1, p2, p4);
    if sign341 == 0.0
        flagInt = true;
        t12 = 0.0;
        t34 = area123 / (area123 - area124);
        pint = p1;
        return
    end
    if sign342 == 0.0
        flagInt = true;
        t12 = 1.0;
        t34 = area123 / (area123 - area124);
        pint = p2;
        return
    end

    % Compute the intersection
    flagInt = true;
    t12  = area341 / (area341 - area342);
    t34  = area123 / (area123 - area124);
    v34  = p4 - p3;
    pint = p3 + v34*t34;

end

% Sign of the oriented area
function signArea = signArea2d(p1,p2,p3)
    det = area2d(p1,p2,p3);
    if abs(det) < 1e-9
        signArea = 0.0;
    else
        if det > 0.0
            signArea = 1.0;
        else
            signArea = -1.0;
        end
    end
end

% Twice the area defined by three points p1, p2 and p3
function area = area2d(p1, p2, p3)
    P2P1 = p2 - p1;
    P3P1 = p3 - p1;
    area = crossProd2D(P2P1,P3P1);
end

% Cross-product between two vectors
function cprod = crossProd2D(v1, v2)
    cprod = v1(1)*v2(2) - v2(1)*v1(2);
end