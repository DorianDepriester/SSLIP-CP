%% Parameters
method='L1'; % should be 'L1' or 'Energy'

% Crystal plasticity parameters (only used for the Energy method)
q=1.702;
h0=309.5;
CRSS=12.22;
sinf=121.8;
a=2.192;

%% Files locations
ebsd_file='EBSD_clean.ctf';
DIC_fold='DIC';

%% Load EBSD data
% Ignore the warning (MTEX-generated file)
CS = crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Copper', 'color', [0.53 0.81 0.98]);
ebsd = EBSD.load(ebsd_file,CS,'interface','ctf');

%% Compute grains
[grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'));

%% Slip systems
sS = slipSystem.fcc(CS);
sS = sS.symmetrise('antipodal');
N=size(sS,1);
sS = sS([6 5 4 11 10 12 7 8 9 3 1 2]); % reorder SSs

%% Initialize new fields
npt=length(ebsd);
nstep=3;
ebsd.prop.gamma=zeros(npt,length(sS),nstep);
ebsd.prop.w=zeros(npt,nstep);
ebsd.prop.res=zeros(npt,nstep);

%% Perform slip estimation for each step
for step=1:nstep
    %% Load gradient values from tabular data
    fname=sprintf('%s/Step_%i.txt', DIC_fold, step);
    S=load(fname);
    npt=size(S,1);
    F=zeros(3,3,npt);
    F(1,1,:)=S(:,1);
    F(2,1,:)=S(:,2);
    F(1,2,:)=S(:,3);
    F(2,2,:)=S(:,4);
    F(3,3,:)=1./(F(1,1,:).*F(2,2,:)-F(1,2,:).*F(2,1,:));
    if step==1
        DeltaF=F-repmat(eye(3),1,1,npt);
    else
        DeltaF=F-ebsd.prop.F.matrix;
    end
    ebsd.prop.F=tensor(F,'rank',2);
    ebsd.prop.DeltaF=tensor(DeltaF,'rank',2);
    if step==1
        ebsd.prop.CRSS=CRSS*ones(npt,N);
    end

    %% Do perform slip estimation
    [ebsd.prop.gamma(:,:,step), ebsd.prop.w(:,step), ebsd.prop.res(:,step), ebsd.orientations, ebsd.prop.CRSS]=computeSSLIPgeneral(ebsd, sS, method, h0, q, sinf, a);

    %% Plot results
    figure
    plotGamma(ebsd, grains, sS, step)
    pause(0.1) % Just wait to let the figure appear
end

%% A enlever
ebsd_66=ebsd(ebsd.grainId==66);
gammas=ebsd.prop.gamma;
gamma66=gammas(ebsd.grainId==66,:,:);
ebsd_66.prop.gamma=gamma66;
figure
plotGamma(ebsd_66, grains(66), sS, 3)