function result = gepfun(X,y,config)
%**************************************************************************
% FUNCTION: result = gepfun(X,y,config)
% INFO:     Run GEP algorithm
% INPUT:    X = design matrix with colums corresponding to variables and 
%               rows corresponding to observations.
%           y = response vector
%           config = configuration structure
% OUTPUT:   result = configuration struct for GEP algorithm 
% AUTHOR:   A. Hensley, 11-Jan-2013
% HISTORY:
%**************************************************************************
% Rev 1.1.      13-Jun-2018       Hensley       Updated for multivariate 
%                                               data.
% Rev 1.0       11-Jan-2013       Hensley       Initial Release
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

%Storage
result.fit = zeros(1,config.trials);
result.fun = cell(1,config.trials);
result.it = zeros(1,config.trials);
result.expr = cell(1,config.trials);
result.mdl = cell(1,config.trials);
result.history = cell(1,config.trials);
result.lastgen = cell(1,config.trials);
result.lastfit = cell(1,config.trials);

for kk = 1:config.trials
    
    %Init Population
    for jj = 1:100
        gen = popgen(config.popsize,config);
        [fun,symb,parsimony] = popexp(gen,config);
        fit = popeval(fun,X,y,config);
        unfit = fit_chk(fit,config);
        %disp(fit)
        if sum(~unfit)>=config.founder_pop
           break 
        end
        disp(['Founder Pop ' num2str(jj) ' Failed'])
    end
    if jj==100
        disp('Unable to Init Founder Population')
        return
    end
    
    %Setup
    history.fit = zeros(1,config.maxgen);
    history.fun = cell(1,config.maxgen);
    history.std = zeros(1,config.maxgen);
    converge = false;
    hw = waitbar(0,'Initializing...');
    figure,
    
    for it = 1:config.maxgen
        
        %Update
        if it>1
            [fun,symb,parsimony] = popexp(gen,config);
            fit = popeval(fun,X,y,config);
        end
        
        %Most Fit
        history.fit(it) = max(fit);
        history.std(it) = std(fit);
        mostfit = find(fit==max(fit),1);
        mdl = fun(mostfit);
        [~,mdl_eval] = popeval(mdl,X,y,config);
        
        %Update Figure
        subplot(2,1,1)
        plot(X,y,'.-'),hold on,plot(X,mdl_eval,'r--'),grid on,hold off
        legend('Data','Model')
        title('Model')
        subplot(2,1,2)
        plot(X,y-mdl_eval,'b.-'),grid on
        title('Error')
        
        %Current Solution
        syms a b c d e f g h
        eval(['f = ' symb{mostfit} ';']);
        history.fun{kk} = vpa(expand(f),3);
        disp(history.fun{kk})
        disp(max(fit))
        
        %Update Waitbar
        waitbar(it/config.maxgen,hw,...
            ['Generation ' num2str(it) ' (Fitness: ' num2str(max(fit)) ')']);
        
        %Chk for Convergence
        if max(fit)>=config.convcrit
            converge = true;
            disp('Algorithm Converged')
            break
        end
        
        %Next Generation
        if it<config.maxgen
            gen = nextgen(gen,fit,config);
        end
    end
    delete(hw)
    %close
    
    %Check for Convergence
    if ~converge
        disp('Algorithm Could Not Converge')
    end
    
    %Save Results/Simplfy Expression
    syms a a b c d e f g h
    eval(['f = ' symb{mostfit} ';']);
    result.expr{kk} = symb{mostfit};
    result.fit(kk) = max(fit);
    result.fun{kk} = vpa(expand(f),3);
    result.it(kk) = it;
    result.mdl{kk} = mdl;
    result.history{kk} = history;
    result.lastgen{kk} = gen;
    result.lastfit{kk} = fit;
    disp('Final Answer:')
    disp(result.expr{kk})
    disp(result.fun{kk})
    
end

