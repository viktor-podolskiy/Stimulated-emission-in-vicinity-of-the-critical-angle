% (c) 2014-2021, University of Massachusetts - Lowell
% non-commercial use only
% see enclosed license

% function calculates layer-specific plane-wave amplitudes and Poynting
% fluxes based on pre-calculated Transfer and Scattering matrices
function [propData]=localTMMcPLMIN(propData,c0TM,c0TE)

Tarr=propData.Tarr; 
Rarr=propData.Rarr; 
Parr=propData.Parr; 
TarrTE=propData.TarrTE; 
RarrTE=propData.RarrTE; 
ParrTE=propData.ParrTE; 

zi=propData.zi; 
nl=length(zi)-1;

cPl=zeros(1,nl+1); 
cMin=zeros(1,nl+1); 
cPlTE=zeros(1,nl+1); 
cMinTE=zeros(1,nl+1); 

%TM waves
phi0=exp(1i*propData.kzArr(1)*zi(1)); 

cPl(1)=c0TM;

for il=1:nl
    Rc=Rarr(il); 
    Tc=Tarr(il); 
    cMin(il)=phi0*Rc*phi0*cPl(il); 
    cPl(il+1)=Tc*phi0*cPl(il); 
    phi0=Parr(il); 
end

%TE waves
phi0=diag(exp(1i*propData.kzArrTE(1)*zi(1))); 

cPlTE(1)=c0TE;

for il=1:nl
    Rc=RarrTE(il); 
    Tc=TarrTE(il); 
    cMinTE(il)=phi0*Rc*phi0*cPlTE(il); 
    cPlTE(il+1)=Tc*phi0*cPlTE(il); 
    phi0=ParrTE(il); 
end

propData.cPlArr=cPl;
propData.cMinArr=cMin; 
propData.cPlArrTE=cPlTE;
propData.cMinArrTE=cMinTE; 
        
end 
