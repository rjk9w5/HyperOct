%% Implementation of newton method

function [x_current] = NewtonMethod(x1,func,tol)

	if(nargin < 3)
		tol = 1e-10;
	end

    x_current = x1;
    x_old = x1*.01;

	%dfunc = @(r,dr) (feval(func,(r+dr)) - feval(func,r))/dr;
    dfunc = @(xo, xn) (feval(func,xn) - feval(func, xo)) / (xn - xo);

	%h = .001;
	it = 0;
	delta = 1;
	while(delta > tol)
		it = it+1;
		%x = x1 - func(x1)/dfunc(x1,h);
        x_new = x_current - func(x_current)/dfunc(x_current,x_old);
		%delta = abs(x1-x)/x;
        
        delta = abs(x_new - x_current) / x_new;

		if(it > 1000000)
			input('1,000,000 iterations. Continue?')
            error("NewtonMethod: Could not converge!")
            break
			tol = 2*tol;
		end
        
        x_old = x_current;
        x_current = x_new;

		%x1 = x; 
	end

end
