function t = tok2(p,A,d1,d2,d3,d4)
	x = p(1);
	y = p(2);
	t = x^4*1./8+A*(1./2*x^2*log(x)-x^4*1./8)+d1+d2*x^2+d3*(x^4-4*x^2*y^2)+d4*y;
end

