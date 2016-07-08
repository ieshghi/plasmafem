function q = parswitch(eps,del,kap)
	a = [[1,(1+eps)^2,(1+eps)^4];[1,(1-eps)^2,(1-eps)^4];[1,(1-del*eps)^2,(1-del*eps)^4-4*(1-del*eps)^2*kap^2*eps^2]];
	b = -1/8*[(1+eps)^4;(1-eps)^4;(1-del*eps)^4];
	
	size(a)
	size(b)	


	q = linsolve(a,b);
end
