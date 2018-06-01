function ngen = nextgen(gen,fit,config)
%**************************************************************************
% FUNCTION: ngen = nextgen(gen,fit,config)
% INFO:     Updates generation
% INPUT:    gen = current population (cell array)
%           fit = current generation fitness (struct)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   ngen = next generation (cell array)
% AUTHOR:   A. Hensley, 28-Nov-2012
% HISTORY:
%**************************************************************************
% Rev 1.0       28-Nov-2012       Hensley       Initial Release
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

%Elitism
mostfit = find(fit==max(fit),1);
safe = gen(mostfit);
leastfit = find(fit==min(fit),1);
gen(leastfit) = [];
fit(leastfit) = [];

%Replication
ngen = rep_fun(gen,fit,config);
clear gen

%Mutation
if config.switch.mutate
    ngen = mutate_fun(ngen,config);
end

%Recombination (1-pt)
if config.switch.recomb1
    ngen = recomb1_fun(ngen,config);
end

%Recombination (2-pt)
if config.switch.recomb2
    ngen = recomb2_fun(ngen,config);
end

%Recombination (Gene)
if config.switch.recombG
    ngen = recombG_fun(ngen,config);
end

%Insertion Sequence Transpose
if config.switch.IST
    ngen = ist_fun(ngen,config);
end

%Root Insertion Sequence Transpose
if config.switch.RIST
    ngen = rist_fun(ngen,config);
end

%Gene Transpose (NOT IMPLEMENTED)
% if config.switch.geneT
%     new = genet_fun(surv,config);
%     ngen = [ngen new]; 
% end

%Sequence Inversion
if config.switch.inv
    ngen = inv_fun(ngen,config);
end

%Constant Transpose
if config.switch.rnc_trnsp && config.switch.rnc
    ngen = rnctrnsp_fun(ngen,config); 
end

%Constant Inversion
if config.switch.rnc_inv && config.switch.rnc
    ngen = rncinv_fun(ngen,config);
end

%Add Most Fit
ngen = [ngen safe];


