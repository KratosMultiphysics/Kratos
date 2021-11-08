function [] = MainCode()

for i = [1:3, 10]
    figure(i)
    hold off;
end

global M;

M = 1.4;

for i = 0:1000
    ThisExists = MakeThisFile([num2str(i), 'drained_triaxial.csv']);

    
    if (ThisExists == false )
        return
    end
    MakeThisFile([num2str(i),'drained_biaxial.csv']);
    MakeThisFile([num2str(i),'undrained_triaxial.csv']);
    MakeThisFile([num2str(i),'oedometer.csv']);
    MakeThisFile([num2str(i),'isotropic.csv']);
end
% MakeThisFile('drained_triaxial.csv')
% MakeThisFile('undrained_triaxial.csv')

% MakeThisFileOed('oedometer.csv')
% MakeThisFileOed('isotropic.csv')


function [] = MakeThisFileOed(XFILE)

rawData = csvread(XFILE);

Stress = rawData(:,2:7);
Strain = rawData(:,8:13);

ps = rawData(:,14);
pt = rawData(:,15);
pc = rawData(:,16);

[p, J] = ComputeStressInvariants( Stress);
[eVol, eDev] = ComputeStrainInvariants( Strain);

% Triaxial plane
figure(1)
plot(p, J);
xlabel('p')
ylabel('q')
axis equal

figure(2)
plot( eDev, J);
xlabel('eDev')
ylabel('sDev')
hold on


figure(3)
semilogx( p, eVol);
xlabel('p')
ylabel('eVol')
hold on



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
PlotYieldSurface( ps(1), pt(1), pc(1), 'k');
PlotYieldSurface( ps(end), pt(end), pc(end), 'g-.');

% Triaxial plane
plot(p, J)
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


figure(10)
plot(-Strain(:,2), -Stress(:,2))
hold on

function PlotYieldSurface( ps, pt, pc, SPEC)

global M;

pp = linspace(0, -pc, 400);
jj = M * sqrt(-pp.*(pp+pc));

plot( pp, M*pp, 'r-.')
hold on
plot(pp+pt, jj, SPEC, 'linewidth', 2.0)


pp = linspace(0, -ps, 400);
jj = M * sqrt(-pp.*(pp+ps));
plot(pp, jj, SPEC, 'linewidth', 1.0)

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
