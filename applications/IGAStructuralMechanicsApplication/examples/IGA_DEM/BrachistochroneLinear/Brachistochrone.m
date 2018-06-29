r = 0.05172;

for s = 0:0.1:4
    s
x = r*(s - sin(s));
y = r*(1 - cos(s));
fprintf('%.15f,%.15f,%.15f\n',x,0.0,-y);
% (x + ', 0.0, ' + y + '\n');
end