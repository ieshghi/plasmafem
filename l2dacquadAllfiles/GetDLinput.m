function [src_loc, src_val, targ_loc, targ_norm] = GetDLinput(p,boxes,pbarc,n,psi,psi_t,psi_s,psi_ts,darc)

% p - grid points (in s,t)
% boxes - boxes associated with p
% pbarc - points on boundary, equispaced in arc-length (in x,y)
% n - inward unit normal on pbarc
% psi, psi_t, psi_s, psi_ts - PDE solution and its derivatives

global xc yc;
global deg xlg wlg;
global RHS;
global fbfunc;

N = deg;

Nbdry = length(pbarc(:,1));
numbox = length(boxes(:,1));


src_loc = zeros(N^2*numbox,2);
src_val = zeros(N^2*numbox,1);

targ_loc = pbarc;
targ_norm = -1.*n;

int_pts = zeros(N^2,2);
int_wts = zeros(N^2,1);
int_pts_xy = zeros(N^2,2);
Jacob = zeros(N^1,1);

for boxind = 1:numbox
        box = boxes(boxind,:);
        boxpts = p(box,:);
        
        % map standard Gauss-Legendre data onto this box
        low_s = boxpts(3,2);
        big_s = boxpts(1,2);
        low_t = boxpts(2,1);
        big_t = boxpts(1,1);
        [s,ws] = lgmap(xlg,wlg,low_s,big_s);
        [t,wt] = lgmap(xlg,wlg,low_t,big_t);

        fb = fbfunc(t);
        
        % Loop over integration points to vectorize them properly
        for k=1:N^2
            i = mod(k-1,N)+1;
            j = floor((k-1)/N)+1;
            int_pts(k,:) = [t(i), s(j)];
            int_wts(k) = wt(i)*ws(j);
            int_pts_xy(k,:) = [(xc + s(j)*fb(i)*cos(t(i))), (yc + s(j)*fb(i)*sin(t(i)))];
            Jacob(k) = s(j)*fb(i)^2;
        end

        %psi at integration points
        psi_int = interpU(boxpts,psi(box),psi_t(box),psi_s(box),psi_ts(box),int_pts);
        % RHS at integration points
        RHSint = RHS(int_pts_xy(:,1),int_pts_xy(:,2),psi_int);
        
        RHSint
        int_wts
        Jacob
        
        integrand = RHSint.*int_wts.*Jacob;
        
        outind_low = N^2*(boxind-1) + 1;
        outind_high = N^2*boxind;
        
        src_loc(outind_low:outind_high,:) = int_pts_xy;
        src_val(outind_low:outind_high,:) = integrand;
        
end
