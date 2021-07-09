clear
%Script Description%

% This script calculates emission profile of the dipole 
% embedded into low-index material in vicinity of high-index substrate

% Geometry of the system:
                       %z=0%       %z>0
%               ----------|----------|----------
%                         |          |
%                 glass   |   PMMA   |   Air
%               (Layer 1) | (Layer 2)| (Layer 3)
%                         |          |
%               ----------|----------|----------

addpath('TMM.Functions');

%% User Input Section%
lrEm=2; %Layer that contains a dipole
nSlices=301; %Number of locations of the dipole to average emission over
lam0=.563;%um -> operating wavelength

dnr=0.125/155; angM=90;%collection angle for sz, will be sind(ang) for data! So 90 deg. is MAX.


% setup geometry/material parameters
hPMMA=2;% thickness of PMMA layer, in microns 
epsPMMA=1.4912^2; 
epsGlass=1.5239^2; 

legs=[num2str(hPMMA),'\mum'];

epsStack=struct(...
    'epsXY',[epsGlass epsPMMA 1], ...
    'epsZZ',[epsGlass epsPMMA 1]);

dh=hPMMA/nSlices;
htPlot=(dh/2:dh:hPMMA);


% Sz calculations.
[szTotArr,nrArr]=szProfile(lam0,epsStack.epsXY,epsStack.epsZZ,hPMMA,lrEm,htPlot, ...
    dnr,angM);

save(['emission.h=',num2str(hPMMA),'um.mat'], 'szTotArr','nrArr','lam0','epsStack','hPMMA','lrEm','htPlot','dnr','angM','dh');

%% Post-Process the data.

figure(5)
clf
szPlt=sum(szTotArr,1)*dh/hPMMA;
plot(asind(nrArr),szPlt/max(szPlt),'linewidth',2);
xlabel('\theta,degree');
ylabel('S(\theta),arb. units');
xlim([0 90])
set(gca,'fontSize',18)
legend([num2str(hPMMA),'\mum'], 'Location','NorthWest')

box on;
grid off;
set(gca,'fontsize',18);


