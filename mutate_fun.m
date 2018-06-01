function ngen = mutate_fun(surv,config)
%**************************************************************************
% FUNCTION: new = mutate_fun(surv,config)
% INFO:     Applies mutation operator to survivors
% INPUT:    surv = survivor chromosomes (cell array)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   ngen = updated generation (cell array)
% AUTHOR:   A. Hensley, 01-Dec-2012
% HISTORY:
%**************************************************************************
% Rev 1.0       01-Dec-2012       Hensley       Initial Release
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
nmut = round(config.mutate*N);
idx = rand_idx(N,nmut);
new = cell(1,nmut);

%Begin
for kk = 1:nmut
    
    %Determine Mutation Points
    cur = surv{idx(kk)};
    pts = rand_idx(length(cur),config.mutate_pts);
    type = config.chrm_map(pts);

    %Apply
    for jj = 1:length(pts)
       
        switch type(jj)
            
            case 'H'
                
                ii = rand_idx(length(config.TF));
                cur(pts(jj)) = config.TF{ii};
                
            case 'T'
                
                ii = rand_idx(length(config.T));
                cur(pts(jj)) = config.T{ii};
                
            case 'C'
                
                ii = rand_idx(length(config.const_set));
                cur(pts(jj)) = config.const_set(ii);
                
            otherwise
                
                error('Bad Chrm Map')
                
        end
        
    end
    
    %Update Output Variable
    new{kk} = cur;
    
end

%Update Population
ngen = surv;
ngen(idx) = new;


