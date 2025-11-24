function [paths_power] = Evaluate_shadowing(paths_power_NoShadowing,elSAT,shadow,scenario,band,propag)

Q=size(paths_power_NoShadowing,1);
paths_power=cell(Q,1);

for ut_IdX=1:Q
    P_qq=length(paths_power_NoShadowing{ut_IdX,1});
    paths_power{ut_IdX}=zeros(P_qq,1);
    for sat_IdX=1:P_qq
        paths_energyTmp=paths_power_NoShadowing{ut_IdX,1}(sat_IdX);
        elev = elSAT{ut_IdX,1}(sat_IdX);
        if ~isnan(paths_energyTmp)
            if ~isnan(elev)
                % codice da aggiungere per il calcolo della potenza di segnale
                chScen = shadowing(scenario,band);
                
                % leggi elevation dal file scenario in gradi
                eIdx = round(abs(elev),-1)/10;      % arrotondato alla decina
                if strcmp(propag,'los')
                    flag = true(size(elev));
                    sigma = chScen.sgL(eIdx(flag));
                    loss = randn(size(eIdx)) .* sigma;
                else
                    flag = 100*rand(size(elev)) <= chScen.lp(eIdx);
                    sigma(flag) = chScen.sgL(eIdx(flag));
                    sigma(~flag) = chScen.sgNL(eIdx(~flag));
                    loss = randn(size(eIdx)) .* sigma;
                    loss(~flag) = loss(~flag) + chScen.cl(eIdx(~flag));
                end
                adl = 10.^(0.1*loss);
                if shadow
                    paths_energy = paths_energyTmp./sqrt(adl); % pathEnergyTmp è la potenza che dovrebbe avere il segnale senza shadowing
                else
                    paths_energy = paths_energyTmp;
                end
                
                paths_power{ut_IdX,1}(sat_IdX)=paths_energy;
            end
        end
        
    end
end
end

