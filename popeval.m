function [fit varargout] = popeval(fun,X,y,config)
%**************************************************************************
% FUNCTION: f = popeval(fun,X,y,config)
% INFO:     Evaluates fitness of current population
% INPUT:    fun = current population (cell array)
%           X = predictors (matrix)
%           y = response (vector)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   f = fitness results(vector)
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

%Error Chk
if nargout>1 && length(fun)>1
    error('Error in popeval: fun variable has too many elements')
end

%Assign Inputs
for jj = 1:length(config.input_args)
    eval([config.input_args{jj} '=X(:,jj);'])
end

%Setup Eval
N = length(fun);
fit = zeros(1,N);
implicitSearch = false;
if all(y==0)
    implicitSearch = true;
end

%Compute Fitness
for kk = 1:N
    eval(['cur = fun{kk}(' config.arg ');']);
    if imag(sum(cur))~=0
        fit(kk) = 0;
        continue
    end
    if implicitSearch
        fit(kk) = fitfun(X,fun{kk},config);
    else
        fit(kk) = fitfun(y,cur,config);
    end
end

%Handle Varargout
if nargout>1
    varargout{1} = cur;
end

