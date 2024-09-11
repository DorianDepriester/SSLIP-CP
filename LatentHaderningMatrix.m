function M = LatentHaderningMatrix(q, n)
% 	Q=q*ones(3,3);
% 	I=ones(3,3);
% 	M=[];
% 	for i=1:4
% 		row=[];
% 		for j=1:4
% 			if i==j
% 				row=[row I];
% 			else
% 				row=[row Q];
% 			end
% 		end
% 		M=[M ; row];
% 	end
% 	M=repmat(M,[2,2]);


theta_min=0.1*degree;
hkl=n.hkl;
T=hkl*hkl';	% Matrix of dot products
% x1=hkl(:,1);
% y1=hkl(:,2);
% z1=hkl(:,3);
% x2=x1';
% y2=y1';
% z2=z1';
% norm_cross_prod= (y1*z2-y2*z1).^2 + z1*x2-z2*



d=sqrt(sum(hkl.^2,2));
dmat=d*d';
theta = acos(T./dmat);	% Must be 0 or pi if perfectly parallel
para = theta < theta_min | pi-theta < theta_min;
M=ones(size(para));
M(~para)=q;
end

