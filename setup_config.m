function config = setup_config(parm)
%**************************************************************************
% FUNCTION: config = setup_config(parm)
% INFO:     Sets up config structure for GEP algorithm
% INPUT:    parm = struct of parameters with the following fields:
%               .nvars = number of input variables (mandatory)
%               .genes = number of genes (optional, default = 1)
%               .headsize = chromosome headsize (optional, default = 10)
%               .rnc = random numerical constant flag (optional, default =
%                   false).
%               .popsize = population size (optional, default = 100);
%               .maxgen = maximum number of generations (optional, default
%                   = 100).
%               .convcrit = convergence criteria between 0 and 1 which is
%                   percent of maximum fitness (1000) required to declare
%                   convergence (optional, default = 0.999)
%               .selthr = minimum finess value for survival (optional,
%                   default = 1e-2).
%               .library = function library cell array (optional, default =
%                   {'+','-','*','/'}). Supported library functions are as 
%                   follows:
%                       '+' = addition
%                       '-' = subtraction
%                       '*' = multiplication
%                       '/' = division
%                       'E' = exp()
%                       'L' = log()
%                       'S' = sin()
%                       'C' = cos()
%                       'T' = tan()
%                       'P' = ()^() 
%                       'Q' = sqrt()
%               .founders = minimim number of survivors for generation 0.
%               .trials = number of times to run GEP algorithm.
%               .linkop = linking operation(s)
%               .mutate = mutation rate
%               .rncfun = numerical constants function handle
%               .fitfun = fitness function(default = mean-square)
%                         {'numhits','mean-square','r-square','complex-mse'}
%
% OUTPUT:   config = configuration struct for GEP algorithm 
% AUTHOR:   A. Hensley, 11-Jan-2013
% HISTORY:
%**************************************************************************
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

%Set Defaults
nvars = nan;
genes = 1;
headsize = 10;
rnc = false;
popsize = 100;
maxgen = 100;
convcrit = 0.9999;
selthr = 1e-2;
library = {'+','-','*','/'};
trials = 1;
founders = 1;
linkop = '+';
mutate = 0.1;
rncfun = @(m,n)2*rand(m,n);
fitfun = 'mean-square';

%Extract Info From Parameter Struct
fields = fieldnames(parm);
for kk = 1:length(fields)
    
    cur = parm.(fields{kk});
    
    switch fields{kk}
        
        case 'nvars'
            nvars = cur;
            
        case 'genes'
            genes = cur;
            
        case 'headsize'
            if headsize<4
                error('Minimum Headsize is 4')
            end
            headsize = cur;
            
        case 'rnc'
            rnc = cur;
            
        case 'popsize'
            popsize = cur;
            
        case 'maxgen'
            maxgen = cur;
            
        case 'convcrit'
            convcrit = cur;
            
        case 'selthr'
            if selthr==0
                error('selthr must be > 0')
            end
            selthr = cur;
            
        case 'library'
            library = cur;
            
        case 'founders'
            founders = cur;
            
        case 'trials'
            trials = cur;
            
        case 'linkop'
            linkop = cur;
            
        case 'mutate'
            mutate = cur;
            
        case 'rncfun'
            rncfun = cur;
            
        case 'fitfun'
            fitfun = cur;
            
        otherwise
            error('Unrecognized Field >> setup_config.m')
            
    end
end

%Error Chk
if isnan(nvars)
    error('Number of variables is undefined >> setup_config.m')
end
if nvars>26
    errordlg('Unable to handle more than 26 preditor variables')
end

%Switches
config.switch.mutate = true;
config.switch.recomb1 = true;
config.switch.recomb2 = true;
config.switch.recombG = true;
config.switch.IST = true;
config.switch.RIST = true;
config.switch.geneT = false;
config.switch.inv = true;
config.switch.rnc = rnc;
config.switch.rnc_trnsp = true;
config.switch.rnc_inv = true;

%Setup Variables & Operators
config.nvars = nvars;
input_args = num2cell(char(97:97+nvars-1));
config.input_args = input_args;
if config.switch.rnc
    config.T = [input_args '?'];
else
    config.T = input_args;
end
config.F = library;
[config.Fmap,maxarg] = set_library(library,'#'); %{'(#+#)','(#-#)','(#*#)','(#/#)'};
config.Ftemp = set_library(library,'$'); %{'($+$)','($-$)','($*$)','($/$)'};
config.TF = [config.F config.T];
config.TFmap = [config.Fmap config.T];
config.TFtemp = [config.Ftemp config.T];
config.maxarg = maxarg;
config.arg = [cell2mat(input_args') repmat(',',size(input_args',1),1)]';
config.arg = config.arg(:)';
config.arg = config.arg(1:end-1);
config.empty = '#';
config.empty_temp = '$';

%Gene Preferences
config.genes = genes;
config.headsize = headsize;
config.tailsize = config.headsize*(config.maxarg-1)+1;
if config.switch.rnc
    config.constsize = config.tailsize;
else
    config.constsize = 0;
end
config.genesize = config.headsize+config.tailsize+config.constsize;
config.chrmsize = config.genesize*config.genes;
config.totalsize = config.chrmsize;
config.linkop = linkop;

%Make Chrm Map
H = repmat('H',1,config.headsize);
T = repmat('T',1,config.tailsize);
C = repmat('C',1,config.constsize);
config.chrm_map = repmat([H T C],1,config.genes);
config.gene_start = 1:config.genesize:config.chrmsize;

%Evolution Preferences
config.trials = trials;
config.popsize = popsize;
config.maxgen = maxgen;
config.mutate = mutate;
config.mutate_pts = 2;
config.recomb1 = 0.3;
config.recomb2 = 0.3;
config.recombG = 0.3;
config.IST_rate = 0.1;
config.ISE_set = 1:3;
config.RIST_rate = 0.1;
config.RIST_set = 1:3;
config.geneT_rate = 0.1;
config.inv_rate = 0.1;
config.inv_set = 2:4;
config.rnctrnsp_rate = 0.1;
config.rnctrnsp_set = 1:3;
config.rncinv_rate = 0.1;
config.rncinv_set = 2:4;
config.founder_pop = founders;

%Convergence/Selection Preferences
config.fitfun = fitfun; %{'numhits','mean-square','r-square','complex-mse'};
config.fitfunmax = 1000;
config.hit_prec = 0.01;
config.selthr = selthr;
config.convcrit = convcrit*config.fitfunmax;

%Constant Setup
config.rncfun = rncfun;
cnum = 10;
config.const_set = num2str(0:cnum-1);
config.const_set(isspace(config.const_set)) = [];
config.empty_const = '?';
if config.switch.rnc
    config.const = config.rncfun(cnum,config.genes);
    %onfig.const = config.crng(1) + config.crng(2)*randn(cnum,config.genes);
else
    config.const = [];
end


%Support Function
function [Fout,varargout] = set_library(Fin,chr)

N = length(Fin);
Fout = cell(1,N);
maxarg = 0;
for kk = 1:N
    cur = Fin{kk};
    switch cur
        
        case '+'
            Fout{kk} = ['(' chr '+' chr ')'];
            maxarg = max(maxarg,2);
            
        case '-'
            Fout{kk} = ['(' chr '-' chr ')'];
            maxarg = max(maxarg,2);
            
        case '*'
            Fout{kk} = ['(' chr '*' chr ')'];
            maxarg = max(maxarg,2);
            
        case '/'
            Fout{kk} = ['(' chr '/' chr ')'];
            maxarg = max(maxarg,2);
            
        case 'E'
            Fout{kk} = ['exp(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'L'
            Fout{kk} = ['log(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'S'
            Fout{kk} = ['sin(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'C'
            Fout{kk} = ['cos(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'T'
            Fout{kk} = ['tan(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'P'
            Fout{kk} = ['(2)^(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'Q'
            Fout{kk} = ['sqrt(abs(' chr '))'];
            maxarg = max(maxarg,1);
            
        case 'Y'
            Fout{kk} = ['real(sqrt(' chr '))'];
            maxarg = max(maxarg,1);
            
        case 'Z'
            Fout{kk} = ['sqrt(' chr ')'];
            maxarg = max(maxarg,1);
            
        case 'A'
            Fout{kk} = ['(' chr ')^2'];
            maxarg = max(maxarg,0);
            
        case 'B'
            Fout{kk} = 'const(2)';
            maxarg = max(maxarg,0);
            
        case 'D'
            Fout{kk} = 'const(3)';
            maxarg = max(maxarg,0);
            
        otherwise
            error('Unrecognized Function')
    end
end

if nargout>1
    varargout{1} = maxarg;
end





