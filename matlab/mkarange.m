function [] = mkarange() 
  as = linspace(-0.8,0.0,10);
  for i=1:10
      i
	  size = 0.04;
      makebdry2(0.32,0.33,1.7,as(i),1.2)
	  mkmesh3(size,(-1)*i)
  end
end
