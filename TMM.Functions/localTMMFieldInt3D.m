% (c) 2014-2021, University of Massachusetts - Lowell
% non-commercial use only
% see enclosed license

function [fld,fldLFC,fldTM,fldTE] = localTMMFieldInt3D(propData,fType)

% calculates the total field at the origin based on pre-calculated
% amplitudes of the kr-dependent "plane-wave" spectrum; integration of the
% field over the polar angle in cylindrical coordinates is done
% analytically

% fType choices: 
% 1=Ez
% 2=Dz
% 3=Ex
% 4=Dx
% 5=Hy

cPlArr=propData.cPlArr; 
cMinArr=propData.cMinArr; 
kzArr=propData.kzArr; 
kr=propData.kr; 
zi=propData.zi; 

cPlArrTE=propData.cPlArrTE; 
cMinArrTE=propData.cMinArrTE; 
kzArrTE=propData.kzArrTE; 

zArr=[0]; 
nArr=[1]; 
nl=length(zi)-1; 
zi(nl+1)=zi(nl); 

% find the layer corresponding to a given z-position
for il=1:nl
    nArr(0>zi(il))=il+1; 
end 
for il=nl:-1:1
    zi(il+1)=zi(il);
end 
zi(1)=0;
z0Arr=zi(nArr);

fldTM=0; fldTE=0; fldLFC=0; %implements local field correction for 

% TM-waves
cPli=cPlArr(:); 
cMini=cMinArr(:); 
kzi=kzArr(:); 
ezi=propData.ezzNL(:); 
switch fType
    case 1
        prefP=-propData.exyNL(nArr)*kr./kzi(nArr)./ezi(nArr); 
        prefM=prefP; 
    case 2
        prefP=-propData.exyNL(nArr)*kr./kzi(nArr); 
        prefM=prefP; 
    case 3
        prefP=1;
        prefM=-prefP;
    case 4
        prefP=propData.exyNL(nArr);
        prefM=-prefP;
    case 5
        prefP=propData.exyNL(nArr)*(propData.omg)./kzi(nArr); 
        prefM=prefP; 
end 
fldTM=fldTM+...
    (prefP.*cPli(nArr).*exp(1i*kzi(nArr).*(zArr-z0Arr))+...
    prefM.*cMini(nArr).*exp(-1i*kzi(nArr).*(zArr-z0Arr)));

%TE-waves
cPli=cPlArrTE(:); 
cMini=cMinArrTE(:); 
kzi=kzArrTE(:); 
switch fType
    case 3
        prefP=-1;
        prefM=-prefP;
    case 4
        prefP=-propData.exyNL(nArr);
        prefM=-prefP;
    case 5
        prefP=-kzi(nArr)./propData.omg; 
        prefM=prefP; 
end 
fldTE=(prefP.*cPli(nArr).*exp(1i*kzi(nArr).*(zArr-z0Arr))+...
    prefM.*cMini(nArr).*exp(-1i*kzi(nArr).*(zArr-z0Arr)));

fld=fldTM+fldTE; 
end

