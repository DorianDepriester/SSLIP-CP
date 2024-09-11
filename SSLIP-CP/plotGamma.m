function plotGamma(ebsd_oI, grains_oI, sS, step)
    cmap=flipud(lbmap(256,'BrownBlue'));
    sSLocalGrain=sS*grains_oI.meanOrientation;
    gammas=ebsd_oI.prop.gamma(:,:,step);
    minmax=3*std(gammas(:));
    N=size(ebsd_oI.prop.gamma,2);
    shape=[3,4];
    mtexFig = newMtexFigure('figSize','large','layout',shape);
    for i=1:N
        [a,b]=ind2sub(shape,i);
        nextAxis(mtexFig,a,b)
        j=sub2ind(flip(shape),b,a);
        j=i;
        plot(ebsd_oI, gammas(:,j),'micronbar','on')
        hold on
        quiver(grains_oI,sSLocalGrain(:,j).trace,'color','gray','linewidth',1)
    %    quiver(grains_oI,sSLocalGrain(:,j).b,'color','gray','linewidth',2)
        plot(grains_oI.boundary)    
        mtexColorMap(mtexFig.gca,cmap)
        mtexTitle(slipSytem2str(sS(j)))
       %nextAxis
    end
    setColorRange('equal')
    CLim(gcm,[-minmax,minmax])
    mtexColorbar
end