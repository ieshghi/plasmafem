function [] = makeplot(del,eps,kap)
	d = parswitch(del,eps,kap);
	c=1;
	N = 10000;

	theta = linspace(0,2*pi,10000);

	tuk = @(x) tok(x,c,d(1),d(2),d(3));

	q = zeros(1,N);

	for i=1:N
		q(i)=find_r(tuk,theta(i),0.7,1e-6,[1,0]);
	end
	
	x = 1 + q.*cos(theta);
	y = q.*sin(theta);
	
	save bdry.txt x y -ASCII
end
