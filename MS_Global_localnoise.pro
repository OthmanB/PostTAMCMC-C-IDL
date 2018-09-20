function MS_Global_localnoise, stat_synthese_freq, stat_synthese_unsorted, plength

lmax=n_elements(stat_synthese_freq[0,*,0])-1
Nmax=n_elements(stat_synthese_freq[0,0,*])

local_noise=dblarr(lmax+1,Nmax)

for i=0,Nmax-1 do begin
	for j=0,lmax do begin
		freq0=Stat_synthese_freq[3,j,i]
		local_noise[j,i]=Harvey_noise(freq0,stat_synthese_unsorted[3, total(plength[0:7]): total(plength[0:8])-1])
	endfor
endfor
	return, local_noise
end

function Harvey_noise, x, noise_params

noise=0.
cpt=0
for i=0,2 do begin
	noise=noise+noise_params[cpt]/(1.+(noise_params[cpt+1]*(1d-3 * x))^noise_params[cpt+2]) ;modele veritable de Harvey
cpt=cpt+3
endfor
	noise=noise + noise_params[n_elements(noise_params)-1]
	return, noise
end