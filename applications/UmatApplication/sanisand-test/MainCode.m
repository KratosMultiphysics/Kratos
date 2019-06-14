function [] = MainCode()

global NumberOfResults
NumberOfResults = 31;

global thisResult;
thisResult = NumberOfResults+1;

for i = [1:3, 2105, 2106, 2107]
    figure(i)
    hold off;
end


for i = 0:1000
    [ThisExists, ThisInfo]= MakeThisFile([num2str(i),'-NoImplex-drained_triaxial.csv']);
    [ThisExists, ThisInfo]= MakeThisFile([num2str(i),'-NoImplex-drained_triaxial_ext.csv']);
    if ( i == 0)
        Info = ThisInfo;
    elseif (ThisExists)
        Info(i+1)=ThisInfo;
    else
        break;
    end
end



Info = ComputeNorms(Info);


figure(2105)
PrintSlope(1e2, 1e-6)
ll = legend('Error$_\sigma$ Explicit', ...
    'Error$_\epsilon$ Explicit', ...
    'location', 'best');
set(ll, 'interpreter','latex')


figure(2106)
ll = legend('Error$_\sigma$ Explicit', ...
    'location', 'best');
set(ll, 'interpreter','latex');

figure(2107)
ll = legend('Explicit', ...
    'location', 'best');
set(ll, 'interpreter', 'latex');



figure(2105)
grid minor
pause(1)
set(gca, 'FontSize', 14)
% print('Triaxial-Element-1', '-dpdf')

figure(2106)
grid minor
set(gca, 'FontSize', 14)
pause(1)
% print('Triaxial-Element-2', '-dpdf')

figure(2107)
grid minor
set(gca, 'FontSize', 14)
pause(1)
% print('Triaxial-Element-3', '-dpdf')

figure(1)

yyy = ylim();
xxx = xlim();
axis equal
ylim([0, yyy(end)]);
xlim(xxx);

% print('Triaxial-PQ', '-dpdf')
figure(2)
ll = legend('Exact', '40', '10', '4', 'location', 'best');
set(ll, 'interpreter', 'latex')

% print('Triaxial-eDevJ', '-dpdf')


function [] = PrintSlope(x, y)
x = [x, x*10];
y = [y, y/10];
plot(x, y,'k', 'linewidth', 1.5)
plot(x(1)*[1,1], y, 'k-.', 'linewidth', 1.1)
plot(x, y(2)*[1,1], 'k-.', 'linewidth', 1.1)

text(x(1)/1.5, y(1)/4, '$1$', 'fontsize', 14, 'HorizontalAlignment', 'Center', 'interpreter', 'latex');
text(x(1)*3, y(1)/20, '$1$', 'fontsize', 14, 'HorizontalAlignment', 'Center', 'interpreter', 'latex');

function Info = ComputeNorms(Info)

n = length(Info);

denom1 = Info(end).stress * Info(end).stress';
denom1 = denom1 + Info(end).ps^2;
denom1 = denom1 + Info(end).pt^2;
denom1 = denom1 + Info(end).pc^2;
denom2 = Info(end).strain * Info(end).strain';

for i = 1:n-1
    Info(i).stressDifference =  Info(i).stress - Info(end).stress;
    Info(i).stressDifference(7) =Info(i).ps - Info(end).ps;
    Info(i).stressDifference(8) =Info(i).pc - Info(end).pc;
    Info(i).stressDifference(9) =Info(i).pt - Info(end).pt;
    Info(i).strainDifference = Info(i).strain- Info(end).strain;
    
    Info(i).Error1 = Info(i).stressDifference*Info(i).stressDifference';
    Info(i).Error1 = sqrt(Info(i).Error1/denom1);
    Info(i).Error2 = Info(i).strainDifference * Info(i).strainDifference';
    Info(i).Error2 = sqrt(Info(i).Error2/denom2);
    
end


return;

figure(2105)

loglog([Info(1:n-1).n], [Info(1:n-1).Error1], '*-.')
hold on
loglog([Info(1:n-1).n], [Info(1:n-1).Error2], '*-.')


xlabel('Number of steps', 'interpreter', 'latex')
ylabel('Error', 'interpreter', 'latex')




figure(2106)
loglog( [Info(1:n-1).ComputationalCost], [Info(1:n-1).Error1], '*-.')
hold on
xlabel('Computational time (s)','interpreter', 'latex')
ylabel('Error','interpreter', 'latex')

figure(2107)
loglog( [Info(1:n-1).n], [Info(1:n-1).ComputationalCost], '*-.')
hold on
xlabel('Number of steps','interpreter', 'latex')
ylabel('Computational time (s)','interpreter', 'latex')




function [ThisExists, Info] = MakeThisFile(XFILE)
global NumberOfResults;
global thisResult;
thisResult = thisResult - 1;

ThisExists = true;

if ( isfile(XFILE) == false)
    Info = 0;
    ThisExists = false;
    return
end
rawData = csvread(XFILE);

if ( rawData(end,1) < 1e-8)
    Info = 0;
    ThisExists = false;
    return;
end

Stress = rawData(:,2:7);
Strain = rawData(:,8:13);

ps = rawData(:,14);
pt = rawData(:,15);
pc = rawData(:,16);
[p, J] = ComputeStressInvariants( Stress);
[eVol, eDev] = ComputeStrainInvariants( Strain);



Info.n = rawData(end,1)+1;
Info.stress = Stress(end,:);
Info.ps = ps(end);
Info.pc = pc(end);
Info.pt = pt(end);
Info.strain = Strain(end,:);
Info.ComputationalCost = rawData(end,end);



wwiiddtthh = 1.1;
SPEC = '-';
if (thisResult == 3) % 1
    SPEC = 'bs-.';
elseif (thisResult == 9) % 10
    SPEC = 'g^-.';
elseif (thisResult == 11)
    SPEC = 'r-.*';
elseif (thisResult == NumberOfResults)
    SPEC = 'm';
    wwiiddtthh = 2;
else
%     return;
end





% Triaxial plane
figure(1)

plot(p, J, SPEC, 'linewidth', wwiiddtthh)
hold on
xlabel('$p$ (kPa)', 'interpreter', 'latex')
ylabel('$q$ (kPa)', 'interpreter', 'latex')

set(gca, 'FontSize', 14)

figure(2)
plot( eDev, J, SPEC, 'linewidth', wwiiddtthh)
hold on
xlabel('$\epsilon_d$', 'interpreter', 'latex')
ylabel('$q$ (kPa)', 'interpreter', 'latex')
set(gca, 'FontSize', 14)

figure(3)
semilogx( -Strain(:,2), -eVol, SPEC,'linewidth', wwiiddtthh)
xlabel('$\epsilon_z$', 'interpreter', 'latex')
ylabel('$\epsilon_v$ ', 'interpreter', 'latex')
hold on
set(gca, 'FontSize', 14)



if ( thisResult == NumberOfResults)
    figure(1)
    PlotYieldSurface( ps(1), pt(1), pc(1), 'k');
    PlotYieldSurface( ps(end), pt(end), pc(end), 'g');
end


function PlotYieldSurface( ps, pt, pc, SPEC)
return;
global M;

pp = linspace(0, -pc, 400);
jj = M * sqrt(-pp.*(pp+pc));

% plot( pp(1:230)+pt, M*pp(1:230), 'r--')
hold on
plot(pp+pt, jj, [SPEC '-'], 'linewidth', 1.0)


pp = linspace(0, -ps, 400);
jj = M * sqrt(-pp.*(pp+ps));
plot(pp, jj, [SPEC '-.'], 'linewidth', 1.0)

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




