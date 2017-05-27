function t = tok(p,c,d1,d2,d3)
	x = p(1);
	y = p(2);
	t = c/8*x^4 + d2*x^2 + d3*(x^4-4*x^2*y^2)+d1;
end

