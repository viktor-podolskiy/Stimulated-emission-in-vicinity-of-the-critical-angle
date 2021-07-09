% (c) 2014-2021, University of Massachusetts - Lowell
% non-commercial use only
% see enclosed license

% calculates transfer matrices for the layered stack for plane wave
% excitation

function [propData,R0,R0te]=localTMMfun(kr, ...
    zi, exyArr, ezzArr,omg)

% cPl, cMin = array of amplitudes of modes propagating in +z and in -z
% directions; 
% propData = structure describing kz-s, and relative amplitudes of the
% modes in each sublayer

% omg=angular frequency/c
% kx = wavevector component along the interface
% zi,exyArr,ezzArr = geometry and permittivity distribution across the
% structure

nl=length(zi)-1;
zi(nl+1)=zi(nl); 

kzArr=zeros(2,nl+1); 
kzArrTE=zeros(1,nl+1); 
sArr=zeros(2,nl+1); 
sArrTE=zeros(1,nl+1); 
exyNL=zeros(1,nl+1); 
ezzNL=zeros(2,nl+1); 

cPl=zeros(2,nl+1); 
cMin=zeros(2,nl+1); 
cPlTE=zeros(1,nl+1); 
cMinTE=zeros(1,nl+1); 


% start the TMM calculations
Ttot1=1; % assume starting from the local medium
Rtot1=0; 

Ttot1TE=1; % assume starting from the local medium
Rtot1TE=0; 

Tarr=cell(nl,1); 
Rarr=cell(nl,1); 
Parr=cell(nl,1); 

TarrTE=zeros(nl,1); 
RarrTE=zeros(nl,1); 
ParrTE=zeros(nl,1); 


%TM fields
for il=nl:-1:1
    [dj1,ej1,kz1T,ez1T,exyT, s1T]=edfun(exyArr(il+1),ezzArr(il+1)); 
    [dj,ej,kz1B,ez1B,exyB,s1B]=edfun(exyArr(il),ezzArr(il)); 

    kzArr(il+1)=kz1T; 
    sArr(il+1)=s1T;
    ezzNL(il+1)=ez1T; 
    exyNL(il+1)=exyT; 
    phi1=exp(1i*kz1T*(zi(il+1)-zi(il))); 
    
    % note: ex-=-ex+, dz-=dz+
    Dj1=dj1*phi1*Rtot1*phi1+dj1; 
    Ej1=-ej1*phi1*Rtot1*phi1+ej1; 

    % implementing the rectangular matrices 
    XM=[-dj, Dj1; ...
        ej, Ej1]; 
    YM=[dj; ej]; 

    mm=[1 1]; 
    RT=mat2cell(XM\YM,mm);
    
    %interface-specific RT
    Rarr{il}=RT{1}; 
    Tarr{il}=RT{2}; 
    Parr{il}=phi1; 
    
    %overall RT
    Rtot1=RT{1}; 
    Ttot1=Ttot1*phi1*RT{2}; 
end
phi0=exp(1i*kz1B*zi(1)); 
kzArr(1)=kz1B; 
sArr(1)=s1B; 
exyNL(1)=exyB; 
ezzNL(1)=ez1B;  

Rall=zeros(length(Rarr),1); 
Pall=Rall; 
Tall=Rall; 
for il=1:length(Rarr)
    Rc=Rarr{il}; 
    Rall(il)=Rc; 
    Pc=Parr{il}; 
    Pall(il)=Pc; 
    Tc=Tarr{il}; 
    Tall(il)=Tc; 
end 

% R0=0;  
R0=Rall(1); 
R0=phi0*R0*phi0; 


%TE waves
for il=nl:-1:1
    [hj1,ej1,kzT,~, sT]=ehfunTE(exyArr(il+1)); 
    [hj,ej,kzB,~,sB]=ehfunTE(exyArr(il)); 

    kzArrTE(il+1)=kzT; 
    sArrTE(il+1)=sT;
    phi1=exp(1i*kzT*(zi(il+1)-zi(il))); 
    
    % note: ey-=-ey+, hx-=hx+
    Hj1=hj1*phi1*Rtot1TE*phi1+hj1; 
    Ej1=-ej1*phi1*Rtot1TE*phi1+ej1; 

    XM=[-hj(1,1), Hj1(1,1); ...
        ej(1,1), Ej1(1,1)]; 
    YM=[hj(1,1); ej(1,1)]; 
    phi1=phi1(1,1); 
    mm=[1 1]; 
    RT=mat2cell(XM\YM,mm);
    
    %interface-specific RT
    RarrTE(il)=RT{1}; 
    TarrTE(il)=RT{2}; 
    ParrTE(il)=phi1; 
    
    %overall RT
    Rtot1TE=RT{1}; 
    Ttot1TE=Ttot1TE*phi1*RT{2}; 
end
phi0=diag(exp(1i*kzB*zi(1))); 
kzArrTE(1)=kzB; 
sArrTE(1)=sB; 

RallTE=RarrTE; 
PallTE=ParrTE; 
TallTE=TarrTE; 

R0te=RallTE(1)*phi0^2; 

% combine the two
propData=struct('omg',omg,'kr',kr,'kzArr',kzArr,'kzArrTE',kzArrTE,...
    'ezzNL',ezzNL,'exyNL',exyNL, ...
    's',sArr, 'sTE',sArrTE,...
    'cPlArr',cPl,'cMinArr',cMin,'cPlArrTE',cPlTE,'cMinArrTE',cMinTE,...
    'zi',zi,...
    'Rarr',Rall,'Tarr',Tall,'Parr',Pall,'R0',R0,...
    'RarrTE',RallTE,'TarrTE',TallTE,'ParrTE',PallTE,'R0TE',R0te...
    ); 


% dispersion and field-dependencies of TM and TE-polarized modes
    function [dM,eM,k1,e1,exyT,s1]=edfun(exy,ezz)
        k1=mySqrt(exy*(omg^2-kr^2/ezz),-pi/2); 
%         if imag(k1)<0
%             k1=-k1;
%         end 
        dM(1,1)=-kr/k1*exy; 
        eM(1,1)=1; 

        e1=ezz; 
        exyT=exy; 
        s1=sqrt(abs((kr/e1/omg)^2+(k1/exyT/omg)^2)); 
    end

    function [hM,eM,k1,exyT,s1]=ehfunTE(exy)
        k1=mySqrt(exy*omg^2-kr^2,-pi/2); 
%         if imag(k1)<0
%             k1=-k1;
%         end 
        hM(1,1)=-k1/omg; 
        eM(1,1)=1; 

        exyT=exy; 
        s1=sqrt(abs(exyT)^2); 
    end 

end

