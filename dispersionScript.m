% (c) 2021, University of Massachusetts - Lowell
% script can be used to calculate the dispersion of the guided modes 
% part of the ASE near total internal reflection project
% script reproduces part of Fig. 5a from APL (2021); see DOI:10.1063/5.0051901

clear 
addpath('TMM.Functions');

% setup the geometry/material parameters
hPMMA=2; % height of PMMA layer
epsPMMA=1.4912^2; % permittivity of doped PMMA layer (gain added later)
epsGlass=1.5239^2; % permittivity of high-index glass

lam0=.563; % operating wavelength

iPMMA=(-4:0.125:-2); % controls range of gain in the system

omg0=2*pi/lam0; %operating frequency/c
zj=[0 hPMMA hPMMA]; %last interface is virtual

% polarization-specific arrays of modal index 
% (modal index represents singular reflection)
% rlog arrays represent lg(1/R) and are used for monitoring purposes
nxTMArr=0*iPMMA; 
nxTEArr=0*iPMMA; 
rlogTM=0*iPMMA; 
rlogTE=0*iPMMA; 

for il=1:length(iPMMA)
    
    eps1j=[epsGlass,real(epsPMMA)-1i*10^iPMMA(il),1]; 
    eps2j=eps1j; 
    
    rTM=@(nx)(refInvFun(nx*omg0,zj,eps1j,eps2j,omg0,true)); 
    rTE=@(nx)(refInvFun(nx*omg0,zj,eps1j,eps2j,omg0,false)); 
    

    % use manual solving process where the 1/R is minimized across the
    % square mesh; mesh parameters are for 2.5um slab, 
    % adjust for other PMMA thicknesses
    nxRe=sind(75.5:0.00125:79.5)*sqrt(epsGlass); 
    nxIm=(-0.0025:0.0000625:0.005); 
    [n2Re,n2Im]=meshgrid(nxRe,nxIm); 
    nxC=n2Re+1i*n2Im; 
    rinvTM=0*n2Re; 
    rinvTE=0*n2Re; 
    parfor in=1:numel(n2Re) %calculate the reflectivity arrays
        rinvTM(in)=log10(abs(rTM(nxC(in)))^2); 
        rinvTE(in)=log10(abs(rTE(nxC(in)))^2); 
    end 
    [~,imin]=min(rinvTM(:)); % find minimum of 1/R for TM polarization
    nxTMArr(il)=nxC(imin); 
    rlogTM(il)=rinvTM(imin); 
    
    [~,imin]=min(rinvTE(:)); % find minimum of 1/R for TE polarization
    nxTEArr(il)=nxC(imin); 
    rlogTE(il)=rinvTE(imin); 
   
    % plot 1/R landscape for current gain level 
    figure(21)
    subplot(1,2,1)
    surf(nxRe,nxIm,(rinvTM),'EdgeColor','none')
    set(gca,'clim',[-2 0])
    colorbar
    view(2)
    title(['TM:',num2str(iPMMA(il))])

    subplot(1,2,2)
    surf(nxRe,nxIm,(rinvTE),'EdgeColor','none')
    set(gca,'clim',[-2 0])
    colorbar
    view(2)
    title(['TE:',num2str(iPMMA(il))])
    drawnow 
end 
%% save
save(['modes.h=',num2str(hPMMA),'um.mat'],'epsGlass','epsPMMA','hPMMA','iPMMA','lam0','omg0','nxTEArr','nxTMArr','rlogTE','rlogTM','nxC')

%% plotting
figure(2)
set(gca,'xscale','log')
set(gca,'FontSize',18)
xlabel('-\epsilon"')
yyaxis right
plot(10.^iPMMA,asind(real(nxTMArr)/sqrt(epsGlass)),...
    10.^iPMMA,asind(real(nxTEArr)/sqrt(epsGlass)), 'LineWidth',2)
ylabel('\theta')
ylim([77 77.2])
yyaxis left
plot(10.^iPMMA,2*omg0*imag(nxTMArr),...
    10.^iPMMA,2*omg0*imag(nxTEArr), 'LineWidth',2)
ylabel('1/L,\mum^{-1}')
ylim([0 0.05])
legend({[num2str(hPMMA),'\mum, TM'],[num2str(hPMMA),'\mum, TE']}, 'Location','SouthWest')

figure(3)
plot(10.^iPMMA,rlogTM,10.^iPMMA,rlogTE, 'LineWidth',2)
set(gca,'xscale','log')
set(gca,'FontSize',18)
xlabel('-\epsilon"')
ylabel('lg(1/R)')

