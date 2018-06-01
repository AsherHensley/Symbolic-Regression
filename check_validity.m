function v = check_validity(x,config)
%**************************************************************************
% FUNCTION: v = check_validity(x,config)
% INFO:     Checks chromosome validity
% INPUT:    x = test chromosomes (cell array)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   v = validity flag
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
v = true;
N = length(x);

%Test
for ii = 1:N
    
    %Head
    curG = x{ii};
    Hmap = config.chrm_map=='H';
    curHead = curG(Hmap);
    if any(~ismember(num2cell(curHead),config.TF))
        v = false;
        break
    end
    
    %Tail
    Tmap = config.chrm_map=='T';
    curTail= curG(Tmap);
    if any(~ismember(num2cell(curTail),config.T))
        v = false;
        break
    end
    
    %Constants
    Cmap = config.chrm_map=='C';
    curConst= curG(Cmap);
    if any(~ismember(curConst,config.const_set))
        v = false;
        break
    end
    
end


