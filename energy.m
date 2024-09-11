function [E, grad, Hess] = energy(g, H, s0)
    gabs=abs(g);
    E=gabs'*(H/2*gabs+s0);
    s=sign(g);
    
    if nargout>1
        grad=(gabs'*H+s0')*diag(s);
    
        if nargout>2
            Hess=H.*sign(g*g');
        end
    end
end
