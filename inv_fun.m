function ngen = inv_fun(surv,config)
%**************************************************************************
% FUNCTION: new = inv_fun(surv,config)
% INFO:     Applies inversion operator to survivors
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
ninv = round(config.inv_rate*N);
[~,idx] = sort(rand(length(surv),1));
idx = idx(1:ninv);
new = cell(1,ninv);
setsize = length(config.inv_set);

%Begin
for kk = 1:ninv
    
    %Select Chromosome
    cur = surv{idx(kk)};
    
    %Select Gene
    gene_idx = rand_idx(config.genes);
    a = config.gene_start(gene_idx);
    head_mask = a:a+config.headsize;
    head = cur(head_mask);
    
    %Select Inversion Sequence
    seq_str = 1;
    seq_end = config.inv_set(rand_idx(setsize));
    mask = seq_str:seq_end;
    mask(mask>length(head)) = [];
    
    %Apply
    head(mask) = fliplr(head(mask));
    
    %Update Chromosome
    cur(head_mask) = head;
    new{kk} = cur;

end

%Update
ngen = surv;
ngen(idx) = new;
