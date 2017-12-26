function r = find_r(func,theta,r_guess,err_tol,center)

% "func" is the function handle for a function on whose zero you
% want to know at r(theta).  "center" is the origin for your polar
% coordinate system.

r = r_guess;        %% initial guess for r
uvec = [cos(theta), sin(theta)];    %% unit vector pointing in the direction of the given theta

x = r*uvec + center;    %% (x,y) coordinates at (r,theta)
curr_err = func(x);     %% error is the size of the function you're trying to find the zero of at the current point
while abs(curr_err) > err_tol
    rp = r+err_tol;
    rm = r-err_tol;
    xp = rp*uvec;
    xm = rm*uvec;
    
    dfr = (func(center+xp) - func(center+xm))/(2*err_tol);  %% computes the directional derivative of func
    
    r = r - curr_err/dfr;       %% update r
    
    x = r*uvec + center;        %% update x
    
    curr_err = func(x);         %% update error
    x
curr_err
    if(imag(curr_err)~=0)
        break
    end
end 
