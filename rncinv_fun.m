function ngen = rncinv_fun(surv,config)
%**************************************************************************
% FUNCTION: ngen = rncinv_fun(surv,config)
% INFO:     Applies inversion to constants domain
% INPUT:    surv = survivor chromosomes (cell array)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   new = transposed chromosomes (cell array)
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
ninv = round(config.rncinv_rate*N);
[~,idx] = sort(rand(length(surv),1));
idx = idx(1:ninv);
new = cell(1,ninv);
setsize = length(config.rncinv_set);

%Begin
for kk = 1:ninv
    
    %Select Chromosome
    cur = surv{idx(kk)};
    
    %Select Gene
    gene_idx = config.gene_start(rand_idx(config.genes));
    a = gene_idx+config.headsize+config.tailsize;
    b = a+config.constsize-1;
    const_mask = a:b;
    const = cur(const_mask);
    
    %Select Insertion Point (Not allowed to be 1)
    pnt = rand_idx(config.constsize,1);

    %Get Sequence
    seq_len = config.rncinv_set(rand_idx(setsize));
    invmask = pnt+(1:seq_len);
    invmask(invmask>config.constsize) = [];
    seq = const(invmask);
    
    %Apply
    const_inv = const;
    const_inv(invmask) = fliplr(const_inv(invmask));
    
    %Update Output Variable
    cur(const_mask) = const_inv;
    new{kk} = cur;

end

%Update
ngen = surv;
ngen(idx) = new;
