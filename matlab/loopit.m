function [] = loopit(m,n)
  for i=m:n
	  size = 0.1/(sqrt(2)^(i-1));
	  mkmesh(size)
  end
end
