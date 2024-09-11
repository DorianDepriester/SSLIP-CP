function [gammas, w, res, new_orientations, varargout] = computeSSLIPgeneral(ebsd, sS , method, varargin)
    ori=ebsd.orientations;
    sSLocal = ori * sS;
    [npt, N]=size(sSLocal);
    gammas=zeros(npt,N);
    res=zeros(npt,1);
    F=ebsd.prop.F;
    DeltaF=ebsd.prop.DeltaF.matrix;
    DeltaF11=DeltaF(1,1,:);
    DeltaF12=DeltaF(1,2,:);
    DeltaF21=DeltaF(2,1,:);
    DeltaF22=DeltaF(2,2,:);
    w=zeros(npt,1);
    CRSS=ebsd.prop.CRSS;
    spinMat = zeros(3,3,npt);
    updt_crss= nargout==5;
    if strcmpi(method,'energy') || updt_crss
        if nargin<7
            error('Please, provide the h0, q, sinf and a')
        else
            h0=varargin{1};
            q=varargin{2};
            sinf=varargin{3};
            a=varargin{4};
            h=LatentHaderningMatrix(q, sS.n);
        end
    end
    parfor i=1:npt
        SchmidMat = sSLocal(i,:).deformationTensor;
        CRSS_i=CRSS(i,:)';
        Mtens=EinsteinSum(SchmidMat,[1,-1,3],F(i),[-1,2]);
        M=Mtens.matrix;
        M11=squeeze(M(1,1,:))';
        M12=squeeze(M(1,2,:))';
        M21=squeeze(M(2,1,:))';
        M22=squeeze(M(2,2,:))';
        A=[M11; M12; M21; M22];
        beq=[DeltaF11(i); DeltaF12(i); DeltaF21(i); DeltaF22(i)];
        gamma0=pinv(A)*beq;
        if strcmpi(method,'energy')
            k1=powerLaw(h0,CRSS_i,h, sinf, a);
            fun= @(x) energy(x, k1, CRSS_i);
            hessianfcn=@(x,lambda) hessianenergy(x, lambda, k1);
        elseif strcmpi(method,'L1')
            fun= @(x) L1norm(x);
            hessianfcn=@(x, lambda) zeros(12);
        else
            error('Unknown method')
        end
        fminconoption = optimoptions('fmincon','Display','none','SpecifyObjectiveGradient',true,'HessianFcn',hessianfcn);        
        [x,f] = fmincon(fun,gamma0, [], [], A, beq, -ones(N,1), ones(N,1), [], fminconoption);
        d=length(gammas(i,:));
        w(i)=f;
        gammas(i,:)=x;
        error_deltaF=A*x-beq;
        res(i)=norm(error_deltaF);
        spin=SchmidMat.antiSym*x;
        spinMat(:,:,i) = spin.matrix;
        c=length(CRSS(i,:));
        
        %% RK4
        if updt_crss
            k1=powerLaw(h0,CRSS_i,h, sinf, a);
            k2=powerLaw(h0,CRSS_i+k1*abs(x)/2,h, sinf, a);
            k3=powerLaw(h0,CRSS_i+k2*abs(x)/2,h, sinf, a);
            k4=powerLaw(h0,CRSS_i+k3*abs(x),h, sinf, a);
            k=(k1+2*k2+2*k3+k4)/6;
            CRSS(i,:)=CRSS_i+k*abs(x);
        end
    end
    new_orientations=ori.* orientation(-spinTensor(spinMat));
    if updt_crss
        varargout={CRSS};
    end
end

