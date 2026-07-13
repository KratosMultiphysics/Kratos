function FractureData = generateRandomFractures(Lx, Ly, nfract, Lmin, varargin)
%GENERATEFRACTURES Uniformly distributed fractures in [0,Lx]x[0,Ly].
%   FractureData = generateFractures(Lx,Ly,nfract,Lmin,Name,Value,...)
%   Output: cell(nfract,1), each entry is [x1 y1; x2 y2].
%   Centers ~ Uniform box. Orientation theta ~ U[0,pi), measured from y-axis CCW.
%   Enforces post-clipping length >= Lmin.
%
%   Name-Value:
%     'Lmax'        double, pre-clip length upper bound. Default 0.5*hypot(Lx,Ly).
%     'Seed'        double or []. rng seed. Default [] (leave RNG).
%     'MaxAttempts' double. Default 1000*nfract.

    if nargin < 4
        error('Provide Lx, Ly, nfract, Lmin.');
    end
    if Lmin <= 0
        error('Lmin must be > 0.');
    end
    if Lmin > hypot(Lx,Ly)
        error('Lmin exceeds domain diagonal.');
    end

    ip = inputParser;
    addParameter(ip,'Lmax',0.5*hypot(Lx,Ly), @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'Seed',[], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
    addParameter(ip,'MaxAttempts',1000*nfract, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip,varargin{:});
    Lmax        = ip.Results.Lmax;
    seed        = ip.Results.Seed;
    maxAttempts = ip.Results.MaxAttempts;

    if Lmax < Lmin
        error('Lmax must be >= Lmin.');
    end
    if ~isempty(seed)
        rng(seed);
    end

    FractureData = cell(nfract,1);
    accepted = 0; attempts = 0;

    while accepted < nfract && attempts < maxAttempts
        attempts = attempts + 1;

        % Uniform center
        cx = Lx * rand();
        cy = Ly * rand();

        % Orientation from y-axis
        theta = pi * rand();                 % U[0,pi)
        u     = [sin(theta), cos(theta)];    % unit direction

        % Pre-clip half-length
        ell = Lmin + (Lmax - Lmin) * rand();
        h   = 0.5 * ell;

        % Endpoints before clipping
        p1 = [cx, cy] - h * u;
        p2 = [cx, cy] + h * u;

        % Clip
        [q1, q2, ok] = clipSegmentToRect(p1, p2, Lx, Ly);
        if ~ok, continue; end

        % Enforce minimum post-clip length
        if norm(q2 - q1) >= Lmin
            accepted = accepted + 1;
            FractureData{accepted} = [q1; q2];
        end
    end

    if accepted < nfract
        warning('Generated %d of %d fractures. Increase MaxAttempts or reduce Lmin.', accepted, nfract);
        FractureData = FractureData(1:accepted);
    end
end

% ---------- helpers ----------
function [q1, q2, ok] = clipSegmentToRect(p1, p2, Lx, Ly)
% Liang–Barsky clipping to [0,Lx]x[0,Ly]
    d = p2 - p1;
    p = [-d(1),  d(1), -d(2),  d(2)];
    q = [ p1(1), Lx - p1(1),  p1(2), Ly - p1(2)];

    u1 = 0; u2 = 1; ok = true;
    for i = 1:4
        if p(i) == 0
            if q(i) < 0, ok = false; return; end
        else
            r = q(i)/p(i);
            if p(i) < 0
                if r > u2, ok = false; return; end
                if r > u1, u1 = r; end
            else
                if r < u1, ok = false; return; end
                if r < u2, u2 = r; end
            end
        end
    end

    q1 = p1 + u1 * d;
    q2 = p1 + u2 * d;

    % Clamp
    q1 = [min(max(q1(1),0.0),Lx), min(max(q1(2),0.0),Ly)];
    q2 = [min(max(q2(1),0.0),Lx), min(max(q2(2),0.0),Ly)];
end
