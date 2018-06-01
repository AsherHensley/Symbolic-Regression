function ngen = recombG_fun(surv,config)
%**************************************************************************
% FUNCTION: new = recombG_fun(surv,config)
% INFO:     Applies gene recombination operator to survivors
% INPUT:    surv = survivor chromosomes (cell array)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   ngen = updated generation (cell array)
% AUTHOR:   A. Hensley, 02-Dec-2012
% HISTORY:
%**************************************************************************
% Rev 1.0       02-Dec-2012       Hensley       Initial Release
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

%Setup
N = length(surv);
nrecomb = 2*round(config.recombG*N);
[~,idx] = sort(rand(length(surv),1));
idx = idx(1:nrecomb);
new = cell(1,nrecomb);

%Begin
for kk = 1:2:nrecomb
    
    %Select Parents
    parentA = surv{idx(kk)};
    parentB = surv{idx(kk+1)};
    
    %Recombination
    temp = config.gene_start(rand_idx(config.genes,1));
    xpt = [temp temp+config.genesize];
    childA = [parentA(1:xpt(1)-1) parentB(xpt(1):xpt(2)-1) parentA(xpt(2):end)];
    childB = [parentB(1:xpt(1)-1) parentA(xpt(1):xpt(2)-1) parentB(xpt(2):end)];
    
    %Update Output Variable
    new{kk} = childA;
    new{kk+1} = childB;
    
end

%Update
ngen = surv;
ngen(idx) = new;

