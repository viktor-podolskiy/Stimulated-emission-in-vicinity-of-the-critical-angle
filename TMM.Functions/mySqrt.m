% (c) 2014-2021, University of Massachusetts - Lowell
function [output] = mySqrt(input,arg0)
%MYSQRT implements square root assuming arguments between arg0...arg0*2*pi

magnitude=abs(input); 
argument=angle(input); %angle returns [-pi...pi]

argument(argument<arg0)=argument(argument<arg0)+2*pi; 
argument(argument>=arg0+2*pi)=argument(argument>=arg0+2*pi)-2*pi; 

output=sqrt(magnitude).*exp(1i*argument/2); 

end

