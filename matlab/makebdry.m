function [] = makebdry(eps,del,kap,c)
	N = 10000;
	theta = linspace(0,2*pi,10000);
	q = zeros(1,N);

	if (del == 0) && (eps==0) && (kap==0)
		tuk = @(x) func(x);
		for i=1:N
			q(i)=1.0;
			%q(i)=find_r(tuk,theta(i),0.7,1e-10,[0,0]); %circle
		end
		x = q.*cos(theta);
		y = q.*sin(theta);
	else
		d = parswitch(eps,del,kap,c);
		tuk = @(x) tok(x,c,d(1),d(2),d(3)); %tokamak
		for i=1:N
			q(i)=find_r(tuk,theta(i),0.7,1e-10,[1,0]); %tokamak
		end
		x = 1 + q.*cos(theta);
		y = q.*sin(theta);
	end
	
	
	save bdry.txt x y -ASCII
	save ../infiles/params.txt eps del kap c -ASCII
end
