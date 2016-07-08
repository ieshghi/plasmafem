function [] = fuck(size)
    my_bdry = @(x,y) (x^4)/8+1+x^2+(x^4-4*x^2*y^2);
    md_dist_func = @(p) dexpr(p,my_bdry);
    fd=md_dist_func;
    pfix=[0,0;0,10;10,0;10,10];
    [p,t]=distmesh2d(fd,@huniform,size,pfix);
    bedge = boundedges(p,t);
    b = unique(bedge);
    save ../infiles/b.txt b -ASCII;
    save ../infiles/p.txt p -ASCII;
    save ../infiles/t.txt t -ASCII;
end

