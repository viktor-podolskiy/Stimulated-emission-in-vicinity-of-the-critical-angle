% (c) 2021, University of Massachusetts - Lowell

function [rInv] = refInvFun(kxArr,zj,epsXY,epsZZ,omg,TMwaves)
%RINVFUN calculates 1/R for a given polarization and layer setup 

[~,RTM,RTE]=localTMMfun(kxArr,zj,epsXY,epsZZ,omg); 
if TMwaves
    rInv=1/RTM; 
else 
    rInv=1/RTE; 
end 

end

