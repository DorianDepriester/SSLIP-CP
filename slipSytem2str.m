function s = slipSytem2str(sS)
    d = [sS(1).b.uvw sS(1).n.hkl];
    d(abs(d) < 1e-10) = 0;
    s='[$';
    for i=1:6
        if d(i)<0
            s = [s sprintf('\\bar{%i}',-round(d(i)))];
        else
            s = [s sprintf('%i',round(d(i)))];
        end
        if i==3
            s = [s ']('];
        end
    end
    s = [s ')$'];
end