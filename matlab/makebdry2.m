function [] = makebdry2(eps,del,kap,c,gamma)
	N = 10000;
	theta = linspace(0,2*pi,10000);
	q = zeros(1,N);

	d = parswitch2(eps,del,kap,c,gamma);
	tuk = @(x) tok2(x,c,d(1),d(2),d(3),d(4)); %tokamak
	for i=1:N
		q(i)=find_r(tuk,theta(i),0.7,1e-10,[1,0]); %tokamak
	end
	x = 1 + q.*cos(theta);
    y = q.*sin(theta);

    d1 = d(1);
    d2 = d(2);
    d3 = d(3);
    d4 = d(4);

	save bdry.txt x y -ASCII
	save ../infiles/params.txt d1 d2 d3 d4 c gamma -ASCII
end
