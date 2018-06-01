function f = fitfun(T,P,config)
%**************************************************************************
% FUNCTION: f = fitfun(T,P,type)
% INFO:     Executes chosen fitness function
% INPUT:    T = Target values (vector)
%           X = Predicted Values (matrix)
%           config = configuration parameters for GEP algorithm (struct)
% OUTPUT:   f = fitness score (vector)
% AUTHOR:   A. Hensley, 06-Dec-2012
% HISTORY:
%**************************************************************************
% Rev 1.0       06-Dec-2012       Hensley       Initial Release
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
K = config.fitfunmax;

%Implicit
if isa(P,'function_handle')
    
    %Only 2 Vars For Now
    if size(T,2)>2
        error('Implicit proc not configured for more than 2 vars')
    end
    
    dT = diff(T);
    x = T(1:end-1,1);
    dx = dT(:,1);
    y = T(1:end-1,2);
    dy = dT(:,2);
    dP_da = (P(x+dx,y)-P(x,y))./dx;
    dP_db = -(P(x,y+dy)-P(x,y))./dy;

    e1 = dT(:,2)./dT(:,1)-dP_da./dP_db;
    e2 = dT(:,1)./dT(:,2)-dP_db./dP_da;
    
    if all(isnan(e1))||all(isnan(e2))
        f = 0;
        return
    end
    
    f = K/(1+nansum(e1.^2)+nansum(e2.^2));
    
    if imag(f)~=0
        f = 0;
    end
    
    return
end

%Error Chk
if any(isnan(T)) || any(isnan(P))
    f = 0;
    return
end
if any(isinf(T)) || any(isinf(P))
    f = 0;
    return
end

%Explicit
switch config.fitfun
    
    case 'relative'
        
        E = mean(sqrt(abs((P-T)./T)));
        
        if isreal(E)
            f = K./(1+E);
        else
            f = 0;
        end
    
    case 'numhits'
        
        E = abs(T-P)<=config.hit_prec;
        f = sum(E/length(T))*K;
    
    case 'mean-square'
        
        E = mean((P-T).^2);
        
        if isreal(E)
            f = K./(1+E);
        else
            f = 0;
        end
        
    case 'r-square'
        
        R = corrcoef(T,P);
        f = R(2,1)^2*K;
        
    case 'complex-mse'
        
        Ei = mean((real(P)-real(T)).^2);
        Er = mean((imag(P)-imag(T)).^2);
        f = K./(1+Er+Ei);
        
    case 'complex-mse-mag'
        
        E = mean(abs(P-T).^2);
        f = K./(1+E);
        
    case 'hamming-mse'
        
        H = hamming(length(T));
        H = H/sum(H);
        E = H'*(abs(P-T).^2);
        f = K./(1+E);
        
    case 'kernel'
        
        N = length(T);
        
        if length(P)==1
            P = P*ones(N,1);
        end

        dT = diff(T(:));
        dP = diff(P(:));

        Kt = dT*(1./dT');
        Kp = dP*(1./dP');
        
        E = sqrt(nanmean(nansum((Kt-Kp).^2)));
        f = K/(1+E);
        
        if ~isreal(f)
            f = 0;
        end
        
    otherwise
        errordlg('Fitness Method Not Implemented')
end

r3 = @(z)round(z*1000)/1000;
if length(unique(r3(P)))==1
    f = 0;
end

