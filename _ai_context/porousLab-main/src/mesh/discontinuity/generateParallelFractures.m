function FractureData = generateParallelFractures(Lx, Ly, theta, deltaX, varargin)
%GENERATEPARALLELFRACTURES Parallel fractures with HORIZONTAL spacing deltaX.
% Inputs:
%   Lx, Ly   : domain size
%   theta    : angle from y-axis (CCW)
%   deltaX   : spacing in x between consecutive lines (constant Δx)
% Name-Value:
%   'OffsetX': horizontal shift of the first line’s x-intercept at y=0 (default 0)
%
% Output:
%   FractureData : cell(N,1), each entry [x1 y1; x2 y2] clipped to [0,Lx]x[0,Ly]

    if any([Lx<=0, Ly<=0, deltaX<=0])
        error('Require Lx>0, Ly>0, deltaX>0.');
    end

    ip = inputParser;
    addParameter(ip,'OffsetX',0,@(x)isnumeric(x)&&isscalar(x));
    parse(ip,varargin{:});
    x0off = ip.Results.OffsetX;

    % Direction (along fracture) and normal
    u = [sin(theta),  cos(theta)];
    n = [cos(theta), -sin(theta)];

    c = cos(theta);
    tol = 1e-12 * max(Lx, Ly);
    if abs(c) < 1e-14
        error('Horizontal spacing undefined for |cos(theta)|≈0 (nearly horizontal fractures).');
    end

    % Signed-distance range over rectangle -> map to x-intercept range at y=0
    C = [0 0; Lx 0; 0 Ly; Lx Ly];
    s_vals = C*n.';                     % distances of corners
    smin = min(s_vals); smax = max(s_vals);
    x0min = min(smin/c, smax/c);
    x0max = max(smin/c, smax/c);

    % Exclude tangential extremes to avoid zero-length segments
    kmin = ceil( (x0min - x0off)/deltaX + tol/deltaX );
    kmax = floor((x0max - x0off)/deltaX - tol/deltaX);
    if kmax < kmin
        FractureData = cell(0,1); return
    end

    T = 2*hypot(Lx,Ly);
    FractureData = cell(kmax - kmin + 1, 1);
    idx = 0;

    for k = kmin:kmax
        x0 = x0off + k*deltaX;      % x-intercept at y=0 kept deltaX apart
        p0 = [x0, 0];               % point on the infinite line
        p1 = p0 - T*u;              % long segment for clipping
        p2 = p0 + T*u;

        [q1, q2, ok] = clipSegmentToRect(p1, p2, Lx, Ly);
        if ok
            if hypot(q2(1)-q1(1), q2(2)-q1(2)) > tol
                idx = idx + 1;
                FractureData{idx} = [q1; q2];
            end
        end
    end

    FractureData = FractureData(1:idx);
end

% --------------- helper ----------------
function [q1, q2, ok] = clipSegmentToRect(p1, p2, Lx, Ly)
    d = p2 - p1;
    p = [-d(1),  d(1), -d(2),  d(2)];
    q = [ p1(1), Lx - p1(1),  p1(2), Ly - p1(2)];
    u1 = 0; u2 = 1; ok = true;

    for i = 1:4
        if p(i) == 0
            if q(i) < 0, ok=false; q1=[NaN,NaN]; q2=[NaN,NaN]; return; end
        else
            r = q(i)/p(i);
            if p(i) < 0
                if r > u2, ok=false; q1=[NaN,NaN]; q2=[NaN,NaN]; return; end
                if r > u1, u1 = r; end
            else
                if r < u1, ok=false; q1=[NaN,NaN]; q2=[NaN,NaN]; return; end
                if r < u2, u2 = r; end
            end
        end
    end

    if u2 <= u1 + eps, ok=false; q1=[NaN,NaN]; q2=[NaN,NaN]; return; end
    q1 = [min(max(p1(1)+u1*d(1),0.0),Lx), min(max(p1(2)+u1*d(2),0.0),Ly)];
    q2 = [min(max(p1(1)+u2*d(1),0.0),Lx), min(max(p1(2)+u2*d(2),0.0),Ly)];
end
