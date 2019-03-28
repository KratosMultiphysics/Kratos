function [] = MainCode()

for i = [1:3]
    figure(i)
    hold off;
end

global M;
global ShapeN;
global SpacingR;

M = 1.0;
ShapeN = 3;
SpacingR = 4;

for i = 0:1000
    ThisExists = MakeThisFile([num2str(i), 'drained_triaxial.csv']);

    
    if (ThisExists == false )
        return
    end
     MakeThisFile([num2str(i),'undrained_triaxial.csv']);
%     MakeThisFile([num2str(i),'undrained_triaxial.csv']);
%     MakeThisFile([num2str(i),'oedometer.csv']);
%     MakeThisFile([num2str(i),'isotropic.csv']);
end
% MakeThisFile('drained_triaxial.csv')
% MakeThisFile('undrained_triaxial.csv')

% MakeThisFileOed('oedometer.csv')
% MakeThisFileOed('isotropic.csv')


function [ThisExists] = MakeThisFile(XFILE)
ThisExists = true;

if ( isfile(XFILE) == false)
    ThisExists = false;
    return
end
rawData = csvread(XFILE);

Stress = rawData(:,2:7);
Strain = rawData(:,8:13);

ps = rawData(:,14);
pt = rawData(:,15);
pc = rawData(:,16);
[p, J] = ComputeStressInvariants( Stress);
[eVol, eDev] = ComputeStrainInvariants( Strain);
% 
% Slope(1) = 0;
% for k = 2:length(eVol)
%     Slope(k) = log(p(k)/p(k-1));
%     Slope(k) = Slope(k)/(eVol(k)-eVol(k-1));
% end
% figure(121)
% semilogy(p, Slope)
% hold on


figure(1)
PlotYieldSurface( ps(1),  'k');
PlotYieldSurface( ps(end), 'g-.');

% Triaxial plane
plot(p, J, 'linewidth', 1.5);
xlabel('p')
ylabel('q')
axis equal

figure(2)
plot( eDev, J);
xlabel('eDev')
ylabel('sDev')
hold on


figure(3)
semilogx( p, -eVol);
xlabel('p')
ylabel('eVol')
hold on


function PlotYieldSurface( p0, SPEC)

global M;
global ShapeN;
global SpacingR;
n = ShapeN;
r = SpacingR;

p0 = -p0;


pp = linspace(0, p0, 1000);

qq = M*pp .* ( -log(pp/p0)/log(r)).^(1/n);
plot(pp, qq, SPEC, 'linewidth', 1.0)
hold on

plot([0, 0.6*p0], M*[0, 0.6*p0], 'r-.')

axis equal




function [p, J] = ComputeStressInvariants( Stress)

p = sum(Stress(:,1:3)')/3;

for i = 1:size(Stress,1)
    J(i) = 0;
    for e = 1:3
        J(i) = J(i) + (Stress(i,e)-p(i))^2;
    end
    for e = 4:6
        J(i) = J(i) + 2*(Stress(i,e))^2;
    end
    J(i) = sqrt(0.5*J(i))*sqrt(3);
end
p = - p;
    

function [eVol, eDev] = ComputeStrainInvariants( Strain)

eVol = sum(Strain(:,1:3)')/3;

for i = 1:size(Strain,1)
    eDev(i) = 0;
    for e = 1:3
        eDev(i) = eDev(i) + (Strain(i,e)-eVol(i))^2;
    end
    for e = 4:6
        eDev(i) = eDev(i) + 2*(Strain(i,e)/2.0)^2;
    end
    eDev(i) = sqrt(0.5*eDev(i))*sqrt(3);
end
eVol = - eVol*3;
