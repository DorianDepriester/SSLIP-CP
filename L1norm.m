function [L1, grad, Hess] = L1norm(g)
    L1=sum(abs(g));
    if nargout>1
        grad=sign(g);
        if nargout>2
            Hess=zeros(size(g));
        end
    end
end

