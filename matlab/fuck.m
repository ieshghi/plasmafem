function [] = fuck(size)
    a = importdata('bdry.txt');
    fd = @(p) dpoly(p,a');

    pfix=[0,-1;2,-1;2,1;0,1];
    [p,t]=distmesh2d(fd,@huniform,size,[0,-1;2,1],pfix);
    bedge = boundedges(p,t);
    b = unique(bedge);

    n=int2str(length(b));

    s1 = '../infiles/';
    mkdir(strcat(s1,n))
    sb = strcat(s1,n,'/b.txt');
    sp = strcat(s1,n,'/p.txt');
    st = strcat(s1,n,'/t.txt');

    save(sb,'b','-ascii');
    save(sp,'p','-ascii');
    save(st,'t','-ascii');
end

