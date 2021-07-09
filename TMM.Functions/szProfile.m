function [szTotArr,nrArr] = szProfile(lam0,exy,ezz,hi,lnum,h0Arr,dnr,angM)

omg=2*pi/lam0; 
nrArr=(dnr/2:dnr:sind(angM));

sfxTMArr=zeros(length(h0Arr),length(nrArr)); 
sfxTEArr=0*sfxTMArr; 
sfzArr=0*sfxTMArr;

    
zi=zeros(1,length(hi)+1);
for il=1:length(hi)
    zi(il+1)=zi(il)+hi(end-il+1);
end
zi=zi-zi(end-lnum+2);

exyArr=fliplr(exy); %emission geometry is flipped with respect to incident geometry
ezzArr=fliplr(ezz); 

krArr=nrArr*sqrt(exy(1))*omg; 
dkr=dnr*sqrt(exy(1))*omg; 

for ih=1:length(h0Arr)

    h0=h0Arr(ih);
    zi=zi+h0;

    % -- emission script --

    % ------------------
    % TMM setup 

    % setting up the field
    sDipole0xTE=zeros(length(nrArr),1); 
    sDipole0xTM=0*sDipole0xTE; 
    sDipole0z=0*sDipole0xTE; 

    % split the layer stack into two, top stack above the dipole, and bottom
    % stack below the dipole. 

    % break the structure into two
    iB=find(zi<0);%negative values in structure. 
    iT=find(zi>=0);%positive values in structure.

    ziB=-fliplr([zi(1),zi(iB)]);
    exyB=fliplr([exyArr(iB),exyArr(iB(end)+1)]); 
    ezzB=fliplr([ezzArr(iB),ezzArr(iB(end)+1)]); 

    ziT=[zi(iT),zi(iT(end))]; 
    exyT=[exyArr(iT),exyArr(iT(end)+1)]; 
    ezzT=[ezzArr(iT),ezzArr(iT(end)+1)]; 

    parfor ikr=1:length(krArr)
        kr1=krArr(ikr); 

        % TMM solution
        [~,rB,rBte]=localTMMfun(kr1,ziB,exyB,ezzB,omg);
        [propDataTop,rT,rTte]=localTMMfun(kr1,ziT,exyT,ezzT,omg);

        for fType=[1,3]

            if fType==3
                % spectral amplitudes of emission 

                % electric dipole Px
                a0T=-1i/2*...
                    kr1*propDataTop.kzArr(1)/propDataTop.exyNL(1); 
                a0B=-a0T; 
                a0Tte=1i*omg^2/2*kr1/propDataTop.kzArrTE(1); 
                a0Bte=-a0Tte; 
            else 

                % electric dipole Pz...
                a0T=1i*kr1^2/propDataTop.ezzNL(1); 
                a0B=a0T; 
                a0Tte=0; 
                a0Bte=a0Tte;
            end 

            % renormalize spectral amplitudes to take into account multiple
            % reflections between top and bottom sub-stacks
            aT=(1-rB*rT)\(a0T+rB*a0B); 

            aTte=(a0Tte+rBte*a0Bte)/(1-rBte*rTte); 

            propDataTop=localTMMcPLMIN(propDataTop,aT,aTte); 

            % field calculation
            [sz,szTM,szTE]=TMMPoyntInt3D(propDataTop,fType); 
            if fType==3
                sDipole0xTM(ikr)=sDipole0xTM(ikr)+szTM*dkr; 
                sDipole0xTE(ikr)=sDipole0xTE(ikr)+szTE*dkr; 
            else 
                sDipole0z(ikr)=sDipole0z(ikr)+sz*dkr; 
            end 
        end 

    end 

    %% ---

  
    %store kr-dependent Poynting flux 
    sfxTMArr(ih,:)=(sDipole0xTM); 
    sfxTEArr(ih,:)=(sDipole0xTE); 
    sfzArr(ih,:)=(sDipole0z); 
    
    zi=zi-h0; 
end 

szTotArr=(2*sfxTMArr+2*sfxTEArr+sfzArr)/3/omg^4; 

end

