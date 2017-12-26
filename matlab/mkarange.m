function [] = mkarange() 
  as = linspace(0.0,1.0,20);
  for i=7:20
      i
	  size = 0.04;
      makebdry2(0.32,0.33,1.7,as(i),1.2)
	  mkmesh3(size,i)
  end
end
