%%% Calcul du volume, de l'aire et de l'épaisseur dans chaque bin de densité
Vol = NaN(720,295,length(bin)-1);
Area = NaN(720,295,length(bin)-1);
Height = NaN(720,295,length(bin)-1);
for i=1:720
	i
	for j=1:295
		for l=1:length(bin_mid)	
			if(bin(l) <= Gamma_max(i,j) & bin(l+1) > Gamma_min(i,j))				
				Area(i,j,l) = dxdy(i,j,1);
				if(isnan(zsurf(i,j,l+1)))
					Height(i,j,l) = depth_bot(i,j) - zsurf(i,j,l);
					Vol(i,j,l) = dxdy(i,j,1)*(depth_bot(i,j) - zsurf(i,j,l));
				elseif(isnan(zsurf(i,j,l)))
					Height(i,j,l) = zsurf(i,j,l+1);
					Vol(i,j,l) = dxdy(i,j,1)*zsurf(i,j,l+1);
				else
					Height(i,j,l) = zsurf(i,j,l+1) - zsurf(i,j,l);
					Vol(i,j,l) = dxdy(i,j,1)*(zsurf(i,j,l+1) - zsurf(i,j,l));
				end
			end		
		end
	end
end

Vol_l = squeeze(nansum(Vol));  %%% Moyenne zonale du volume binné en densité
Vol_k = squeeze(nansum(dxdydz));  %%% Moyenne zonale du volume en z

%%% Calcul de la somme cumulée à partir du fond des deux volumes ci-dessus
Vol_rk = NaN(295,44);
for j=1:295
	for k=1:44
		Vol_rk(j,k) = nansum(Vol_k(j,45-k:44),2);
	end
end
Vol_rl = NaN(295,length(bin_mid));
for j=1:295
	for l=1:length(bin_mid)
		Vol_rl(j,l) = nansum(Vol_l(j,length(bin_mid)+1-l:length(bin_mid)),2);
	end
end
Vol_rl(Vol_rl==0)=NaN;
Vol_rk(Vol_rk==0)=NaN;

%%% Calcul de la prodondeur des surfaces de densité correspondant au remplissage ordonné en densité depuis le fond de l'océan
Iso_depth = NaN(295,length(bin_mid));
for j=1:295
	for l=1:length(bin_mid)
		K = find(Vol_rk(j,:) <= Vol_rl(j,l),1,'last');
		if(isempty(K))
			Iso_depth(j,l) = zax(45) - (zax(45) - zax(44))*Vol_rl(j,l)/Vol_rk(j,1);
		elseif(K==44)
			Iso_depth(j,l) = 0;
		else
			Iso_depth(j,l) = zax(45-K) - (zax(45-K) - zax(45-K-1))*(Vol_rl(j,l) - Vol_rk(j,K))/(Vol_rk(j,K+1) - Vol_rk(j,K));
		end
	end
end
Iso_depth = fliplr(Iso_depth);

%%% Exemple de champ 'nrho_lz' prêt à être plotté (moyenne zonale le long des isopycnes, remappée en z)
axis_lz = 0:5:6000;
nrho_lz = NaN(295,length(axis_lz));
for j=1:295
	for l=1:length(bin_mid)
		ind = axis_lz <= nanmax(nanmax(depth(:,j,:),[],3),[],1) & axis_lz >= Iso_depth(j,l);
		nrho_lz(j,ind) = bin_mid(l);
	end
end

