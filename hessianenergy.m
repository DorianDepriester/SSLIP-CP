function Hout = hessianenergy(gamma, lambda, H, varargin)
%Mahalanobis(               gamma, C, E, SchmidMat, F, DeltaF, th)
    Hout=H.*sign(gamma*gamma');
    
    if nargin>3
        C=varargin{1};
        E=varargin{2};
        SchmidMat=varargin{3};
        F=varargin{4};
        DeltaF=varargin{5};
        A=zeros(4);
        B=A;
        fbar=zeros(4,1);
        f=DeltaF(1:2,1:2);
        f=f(:);
        Ft=F(1:2,1:2);
        Ft=Ft(:);
        ns=size(SchmidMat,3);
        for alpha=1:ns
            Sat=kron(eye(2), SchmidMat(1:2,1:2,alpha));
            A=A+gamma(alpha)*(C*Sat'+Sat*C);
            fbar=fbar-gamma(alpha)*Sat*E;
            f=f-gamma(alpha)*Sat*Ft;
            for beta=1:ns
                Sb=kron(eye(2), SchmidMat(1:2,1:2,beta));
                Sbt=kron(eye(2), Sb(1:2,1:2));
                B=B+gamma(alpha)*gamma(beta)*Sat*C*Sbt';
            end
        end
        cov=4*C+2*A+B;

        Cinv=C^(-1);
        Hess=zeros(12,12);
        S=cov^(-1);
        gradc2=zeros(12,1);
        for i=1:12
            Sit=kron(eye(2), SchmidMat(1:2,1:2,i));
            gradc2(i)=(f-fbar)'*(2*S*Sit*(E-Ft)-1/8*(C^-1*Sit+(C^-1*Sit)')*(f-fbar));
            Ki=Cinv*Sit+Sit'*Cinv;
            for j=i:12
                Sjt=kron(eye(2), SchmidMat(1:2,1:2,j));
                Kj=Cinv*Sjt+Sjt'*Cinv;
                T=Sjt'*Ki+Sit'*Kj;
                Hij=(E-Ft)'*( 2*Sjt'*S*Sit*(E-Ft) - 1/4*T*(f-fbar) );
                Hess(i,j)=Hij;
                Hess(j,i)=Hij;
            end
        end
        Hout = Hout + lambda.ineqnonlin*Hess;
    end
end

