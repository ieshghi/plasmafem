function [] = oldfuck(size)
    fd=inline('drectangle(p,1,2,1,2)','p');
    pfix=[1,1;1,2;2,1;2,2];
    [p,t]=distmesh2d(fd,@huniform,size,[1,1;2,2],pfix);
    bedge = boundedges(p,t);
    b = unique(bedge);
    save ../infiles/b.txt b -ASCII;
    save ../infiles/p.txt p -ASCII;
    save ../infiles/t.txt t -ASCII;
end
