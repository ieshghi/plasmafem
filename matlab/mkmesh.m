function [] = mkmesh(size)
    a = importdata('bdry.txt');
    params = importdata('../infiles/params.txt');
    fd = @(p) dpoly(p,a');

    pfix=[0,-1;2,-1;2,1;0,1];
    [p,t]=distmesh2d(fd,@huniform,size,[0,-1;2,1],pfix);
    bedge = boundedges(p,t);
    b = unique(bedge);
    nt = length(t(:,1))
    hs = zeros([1,nt]);

    for i=1:nt
            t1 = sqrt((p(t(i,1),1)-p(t(i,2),1))^2+(p(t(i,1),2)-p(t(i,2),2))^2);
            t2 = sqrt((p(t(i,3),1)-p(t(i,2),1))^2+(p(t(i,3),2)-p(t(i,2),2))^2);
            t3 = sqrt((p(t(i,1),1)-p(t(i,3),1))^2+(p(t(i,1),2)-p(t(i,3),2))^2);
	    hs(i) = max([t1,t2,t3]);
    end
    edge = max(hs)


    n=int2str(length(b));

    s1 = '../infiles/';
    mkdir(strcat(s1,n))
    sb = strcat(s1,n,'/b.txt');
    sp = strcat(s1,n,'/p.txt');
    st = strcat(s1,n,'/t.txt');
    sh = strcat(s1,n,'/h.txt');
    sr = strcat(s1,n,'/params.txt');

    save(sb,'b','-ascii');
    save(sp,'p','-ascii');
    save(st,'t','-ascii');
    save(sh,'size','-ascii');
    save(sr,'params','-ascii');
end

