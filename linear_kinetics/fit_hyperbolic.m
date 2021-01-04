function [fitobj, f, dfdx, fp] = fit_hyperbolic(Xin, Yin)
% FIT_HYPERBOLIC computes the parameters of the hyperbolic function that
% fits the input data

% [FITOBJ F DFDX FP] = FIT_HYPERBOLIC(XIN,YIN) takes the independenet
% variable vector XIN and the dependent variable vector YIN and computes
% the parameters of a hyperbolic function y = f(x) = a - b/(x+c)
% algebraically (rather than by least-squares minimzation). This is done
% because Yin are in fact hyperbolic function of Xin in this model. NOTE:
% Numerical difficulties arise when the hyperbolic function has very low
% curvature; in this case, a linear function is fitted (by linear
% regression) instead.
% 
% OUTPUTS:
%
% FITOBJ -- cfit object describing the fit
%
% F -- function handle describing the fitted hyperbolic function
% 
% DFDX -- function handle describing the derivative
%
% FP -- fixed points (at most two). Fixed points are ordered by absolute
% value from smallest to largest



% Use only values that are not nan in the fit
TF = ~isnan(Yin);

X = Xin(TF);
Y = Yin(TF);

% Cannot fit if there are too few points
if nnz(TF) < 3
    fitobj = nan;
    f = nan;
    dfdx = nan;
    fp = nan;
    return;
end



% Y data can have a huge range which can cause issues with fitting. Rescale
% Y data to the absolute value of the maximum. Will rescale the fitted
% parameters back

scale_factor = max(abs(Y));
Y = Y / scale_factor;

x1 = X(1:end-2);
x2 = X(2:end-1);
x3 = X(3:end);

y1 = Y(1:end-2);
y2 = Y(2:end-1);
y3 = Y(3:end);


% The auxiliary alpha parameter determines whether we fit a hyperbolic or
% linear function

alpha = ((y3-y1)./(x3-x1)) ./ ((y2-y1)./(x2-x1));

if min(abs(alpha - 1)) < 1e-3 || ( any(alpha - 1 < 0) && any(alpha - 1 > 0) )
    coeffs = regress(Y, [ones(size(X)), X]);
    a = coeffs(1) * scale_factor;
    b = coeffs(2) * scale_factor;
        
    % Fixed point calculation
    fp = a/ (1-b);

    ft = fittype('a + b*x');
    fitobj = cfit(ft, a, b);
    
    f = @linear_f;
    dfdx = @linear_dfdx;

else
    c = (alpha .* x3 - x2) ./ (1-alpha);
    b = - (y3 - y1) ./ ( 1./(x3+c) - 1./(x1+c) );
    a = y1 + b ./ (x1 + c);
    
    a = mean(a);
    b = mean(b);
    c = mean(c);
        
    a = a * scale_factor;
    b = b * scale_factor;

    % Fixed point calculation
    fp = nan(1,2);
    
    if (a - c)^2 - 4 * (b - a * c) >= 0
        fp(1) = 0.5 * (a - c + sqrt( (a - c)^2 - 4 * (b - a * c)  ) );
        fp(2) = 0.5 * (a - c - sqrt( (a - c)^2 - 4 * (b - a * c)  ) );
        
        if abs(fp(2)) < abs(fp(1))
            fp = flip(fp);
        end
    end
        
    ft = fittype('a - b/(x + c)');
    fitobj = cfit(ft, a, b, c);
    
    f = @hyperbolic_f;
    dfdx = @hyperbolic_dfdx;
   
end



    function r = hyperbolic_f(x)
        r = a - b ./ ( x + c );
    end

    function r = hyperbolic_dfdx(x)
        r = b ./ ( x + c ).^2;
    end

    function r = linear_f(x)
        r = a + b * x;
    end

    function r = linear_dfdx(x)
        r = b;
    end


end

