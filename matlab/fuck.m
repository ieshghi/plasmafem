function [] = fuck(size)
    a = importdata('bdry.txt');
    fd = @(p) dpoly(p,a');

    pfix=[0,-1;2,-1;2,1;0,1];
    [p,t]=distmesh2d(fd,@huniform,size,[0,-1;2,1],pfix);
    bedge = boundedges(p,t);
    b = unique(bedge);
    save ../infiles/b.txt b -ASCII;
    save ../infiles/p.txt p -ASCII;
    save ../infiles/t.txt t -ASCII;
end

