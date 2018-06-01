function ngen = rep_fun(gen,fit,config)
%**************************************************************************
% FUNCTION: ngen = rep_fun(gen,fit,config)
% INFO:     Applies replication operator
% INPUT:    gen = current population (cell array)
%           fit = current generation fitness (struct)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   ngen = next generation (cell array)
% AUTHOR:   A. Hensley, 06-Dec-2012
% HISTORY:
%**************************************************************************
% Rev 1.0       06-Dec-2012       Hensley       Initial Release
%**************************************************************************
%   MIT License
%   
%   Permission is hereby granted, free of charge, to any person obtaining a 
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
% 
%   The above copyright notice and this permission notice shall be included 
%   in all copies or substantial portions of the Software.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Prune
unfit = fit_chk(fit,config);
fitF = fit;
genF = gen;
fitF(unfit) = [];
genF(unfit) = [];

if isempty(fitF)
   disp('Population Extinct')
   return
end

%Replicate
fitpdf = fitF/sum(fitF);
fitcdf = [0 cumsum(fitpdf)];
fitsup = 0:length(fitpdf);
new = rand(1,length(fit));
idx = ceil(interp1(fitcdf,fitsup,new));
ngen = genF(idx);

