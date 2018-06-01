function [f s parsimony] = popexp(x,config)
%**************************************************************************
% FUNCTION: f = popexp(x,config)
% INFO:     Converts chromosomes to functional expressions
% INPUT:    x = chromosomes (cell array)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   f = function handles (cell array)
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

%Setup
N = length(x);
f = cell(1,N);
s = cell(1,N);
parsimony = zeros(1,N);
for kk = 1:N
    
    %Segment/Reshape
    curG = reshape(x{kk}(1:config.chrmsize)',[],config.genes)';
    expG = cell(1,config.genes);
    
    for jj = 1:config.genes
        
        %Current Gene
        z = curG(jj,:);
        
        %Set Counter
        ii = 1;
        
        %Begin Gene Expression
        while ii<=config.headsize+config.tailsize
            
            %Lookup T/F
            cursym = symconv(z(ii),config);
               
            if ii==1
                
                %Init
                ft = ['@' cursym];
                ii = ii+1;
                
            else
                
                %Fill Empty Slots
                idxA = find(ft==config.empty,1);
                if ~isempty(idxA)
                    ft = [ft(1:idxA-1) cursym ft(idxA+1:end)];
                    ii = ii+1;
                else
                    
                    %Update Place Holders
                    idxB = strfind(ft,config.empty_temp);
                    if ~isempty(idxB)
                        ft(idxB) = config.empty;
                    else
                       break 
                    end
                    
                end
            end
            
        end
        parsimony(kk) = parsimony(kk)+ii-1;
        
        %Chk for Constants
        const_loc = strfind(ft,config.empty_const);
        if ~isempty(const_loc)
            c_idx = curG(jj,config.headsize+config.tailsize+1:end);
            c_val = config.const(:,jj);
            for uu = 1:length(const_loc)
                const_loc = const_loc(1);
                c = c_val(str2double(c_idx(uu))+1);
                c = ['(' num2str(c) ')'];
                ft = [ft(1:const_loc-1) c ft(const_loc+1:end)];
                const_loc = strfind(ft,'?');
            end
        end
        
        %Save Current Gene
        expG{jj} = ft(2:end);
        
    end
    
    %Assemble Final Expression
    expL = [];
    for ii = 1:config.genes-1
        expL = [expL expG{ii} config.linkop];
    end
    expL = [expL expG{end}];
    f{kk} = eval(['@(' config.arg ')' vectorize(expL)]);
    s{kk} = expL;
    
    %Error Chk
    if nargin(f{kk})~=config.nvars
       error('Bad Expression Generated') 
    end
end


