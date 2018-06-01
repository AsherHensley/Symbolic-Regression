function ngen = rist_fun(surv,config)
%**************************************************************************
% FUNCTION: new = rist_fun(surv,config)
% INFO:     Applies root insertion sequence transpose operator to survivors
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
ntransp = round(config.RIST_rate*N);
[~,idx] = sort(rand(length(surv),1));
idx = idx(1:ntransp);
new = cell(1,ntransp);
setsize = length(config.RIST_set);

%Begin
for kk = 1:ntransp
    
    %Select Chromosome
    cur = surv{idx(kk)};
    
    %Select Gene
    gene_idx = config.gene_start(rand_idx(config.genes));
    head_mask = gene_idx:gene_idx+config.headsize-1;
    head = cur(head_mask);
    
    %Select Insertion Point (Not allowed to be 1)
    pnt = rand_idx(config.headsize-1,1)+1;
    
    %Find Next Function
    temp = head;
    temp(1:pnt-1)='#';
    pnt = find(ismember(temp,cell2mat(config.F)),1);
    if isempty(pnt)
        new{kk} = cur;
        continue
    end

    %Get Transposon
    mask = pnt+((1:rand_idx(setsize,1))-1);
    mask(mask>config.headsize) = [];
    ris = head(mask);
    
    %Apply
    headT = [ris head];
    headT = headT(1:config.headsize);
    
    %Update Output Variable
    cur(head_mask) = headT;
    new{kk} = cur;

end

%Update
ngen = surv;
ngen(idx) = new;

