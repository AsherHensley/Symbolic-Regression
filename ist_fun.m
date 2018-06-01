function ngen = ist_fun(surv,config)
%**************************************************************************
% FUNCTION: new = ist_fun(surv,config)
% INFO:     Applies insertion sequence transpose operator to survivors
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
ntransp = round(config.IST_rate*N);
[~,idx] = sort(rand(length(surv),1));
idx = idx(1:ntransp);
new = cell(1,ntransp);
setsize = length(config.ISE_set);

%Begin
for kk = 1:ntransp
    
    %Select Chromosome
    cur = surv{idx(kk)};
    
    %Select Gene
    gene_idx = rand_idx(config.genes);
    a = config.gene_start(gene_idx);
    gene_mask = a:a+config.headsize+config.tailsize-1;
    gene = cur(gene_mask);
    
    %Select Transposon
    mask = (1:rand_idx(setsize,1))-1;
    head_tail = config.headsize+config.tailsize;
    pnt = rand_idx(head_tail-length(mask)+1);
    is = gene(pnt+mask);
    
    %Select Insertion Point (Not allowed to be 1)
    ins = rand_idx(config.headsize-1)+1;
    
    %Apply
    head = gene(1:config.headsize);
    headT = [head(1:ins-1) is head(ins:end)];
    headT = headT(1:config.headsize);
    
    %Update Chromosome
    gene(1:config.headsize) = headT;
    cur(gene_mask) = gene;
    new{kk} = cur;

end

%Update
ngen = surv;
ngen(idx) = new;

