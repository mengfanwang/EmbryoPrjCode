function h = gaussianKernel3D(sigma, r)
%This function creates a pre-defined 3-D Gaussian kernel.
%
%   Description :
%  
%   h = gaussianKernel3D(r,sigma) returns a rotationally symmetric Gaussian
%   kernel h of size 2*r+1 with standard deviation sigma (positive). The 
%   default size of r is 1 and the default sigma is 0.5.
%
%   Author : Mohammad Mustafa
%   By courtesy of The University of Nottingham and Mirada Medical Limited,
%   Oxford, UK
%
%   Published under a Creative Commons Attribution-Non-Commercial-Share Alike
%   3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%   
%   June 2012
    if nargin==1
        r=ceil(sigma*sqrt(-2*log(sqrt(2*pi)*sigma*(1e-2))));
    end
        
    [x,y,z] = meshgrid(-r:r,-r:r,-r:r);
    arg   = -(x.*x + y.*y + z.*z)/(2*sigma*sigma);
    h = exp(arg);
    
    if sum(h(:)) ~= 0,
      h  = h/sum(h(:));
    end;
end
