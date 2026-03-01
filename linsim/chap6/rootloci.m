p=-1:0.01:1

  for j=1:length(p)
	  ev(:,j)=eig([0 1;p(j) 0]);
end
