function [Rho, Rho_tilde] = Rho_Rho_Tilde_Evaluation_p681(paths_energy,Delays,chIT)
Q=size(paths_energy,1);
Rho=cell(Q,1);
Rho_tilde=cell(Q,1);

for qq=1:Q
    P_qq=size(paths_energy{qq,1},1);
    paths_energy_qq=paths_energy{qq,1}(:,chIT);
    Delays_qq=Delays{qq,1};
    Rho_qq=paths_energy_qq.*exp(1j*2*pi*rand(P_qq,1)); % Uniform Phase gains
    Rho_tilde_qq=Rho_qq.*exp(1j*2*pi*Delays_qq);
    
    Rho{qq,1}=Rho_qq;
    Rho_tilde{qq,1} = Rho_tilde_qq;
    
end
end

