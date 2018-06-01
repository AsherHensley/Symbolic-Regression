function new = genet_fun(surv,config)
%**************************************************************************
% FUNCTION: new = genet_fun(surv,config)
% INFO:     Applies gene transpose operator to survivors
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

%Config Chk
if config.genes==1
    disp('WARNING: Not enough genes to perform gene transpose operator')
    new = [];
    return
end

%Setup
N = length(surv);
ntransp = round(config.geneT_rate*N);
idx = rand_idx(N,ntransp);
new = cell(1,ntransp);

%Begin
for kk = 1:ntransp
    
    %Select Chromosome
    cur = surv{idx(kk)};
    
    %Select Gene
    gene_num = rand_idx(config.genes-1,1)+1;
    gene_idx = config.gene_start(gene_num);
    mask = gene_idx:gene_idx+config.genesize-1;
    cur_gene = cur(mask);
    
    %Update Output Variable
    cur(mask) = [];
    new{kk} = [cur_gene cur];

end

