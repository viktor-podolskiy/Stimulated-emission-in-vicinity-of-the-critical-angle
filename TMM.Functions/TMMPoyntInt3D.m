% (c) 2014-2019, University of Massachusetts - Lowell
% non-commercial use only
% see enclosed license

function [sz,szTM,szTE] = TMMPoyntInt3D(propData,fType)

% calculates the total Poynting flux in the outer-most layer

% fType choices: 
% 1=z-oriented dipole
% 3=x-oriented dipole

cPlArr=propData.cPlArr; 
cMinArr=propData.cMinArr; 
kzArr=propData.kzArr; 
kr=propData.kr; 
zi=propData.zi; 
omg=propData.omg; 

cPlArrTE=propData.cPlArrTE; 
cMinArrTE=propData.cMinArrTE; 
kzArrTE=propData.kzArrTE; 

% zArr=[0]; 
nl=length(zi)-1; 
zi(nl+1)=zi(nl); 
nArr=[nl+1]; 

szTM=0; szTE=0; 

% TM-waves
cPli=cPlArr(:); 
kzi=kzArr(:); 
exy=propData.exyNL; 
switch fType
    case 1
        szTM=szTM+2*pi*abs(cPli(nArr).^2)*omg/kr*real(exy(nArr)/kzi(nArr)); 
    case 3
        szTM=szTM+4*pi*abs(cPli(nArr).^2)*omg/kr*real(exy(nArr)/kzi(nArr));
end 

%TE-waves
cPli=cPlArrTE(:); 
cMini=cMinArrTE(:); 
kzi=kzArrTE(:); 
if fType==3
    szTE=szTE+4*pi*abs(cPli(nArr).^2)*real(kzi(nArr))/omg/kr; 
end 

sz=szTM+szTE; 
end

