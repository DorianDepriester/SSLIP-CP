function h = powerLaw(h0, CRSS, latent_matrix, satu_stress, a)
	n_ss=length(CRSS);
    k=zeros(size(CRSS));
    k(CRSS<satu_stress)=h0.*(1-CRSS(CRSS<satu_stress)./satu_stress).^a;
	h=repmat(k',[n_ss 1]).*latent_matrix;
end

