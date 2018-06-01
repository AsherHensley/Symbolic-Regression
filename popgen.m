function x = popgen(n,config)
%**************************************************************************
% FUNCTION: x = popgen(n,config)
% INFO:     Generates n candidtate solutions based config preferences
% INPUT:    n = number of candidate solutions to generate (scalar)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   x = candiate solution chromosomes (cell array)
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

%Init
x = cell(1,n);
for kk = 1:n
    for jj = 1:config.genes
        
        %Generate Head
        h = config.headsize;
        nh = length(config.TF);
        cur_head = [config.TF{rand_idx(nh,h)}];
        
        %Generate Tail
        t = config.tailsize;
        nt = length(config.T);
        cur_tail = [config.T{rand_idx(nt,t)}];
        
        %Generate Constants
        if config.switch.rnc
            nc = length(config.const_set);
            cur_const = config.const_set(rand_idx(nc,t));
        else
            cur_const = [];
        end
        
        %Assemble Chromosome
        x{kk} = [x{kk} cur_head cur_tail cur_const];
    end
    
    
end
