function [Antenna_pos] = Antenna_pos_generation(Nx,Ny,spacing)
% positions on the (x,y) plane with z_pos=0;

if Nx*Ny==1
    
    Antenna_pos=[0,0,0];

else


    Nx_pos_positive=[];
    for nx=1:Nx/2
        Nx_pos_positive=[Nx_pos_positive; (nx-1)*spacing+spacing/2];
    end
    Nx_pos=[Nx_pos_positive; -Nx_pos_positive];
    Ny_pos_positive=[];
    for ny=1:Ny/2
        Ny_pos_positive=[Ny_pos_positive; (ny-1)*spacing+spacing/2];
    end
    Ny_pos=[Ny_pos_positive; -Ny_pos_positive];
    Antenna_pos=[];
    for nx=1:length(Nx_pos)
        for ny=1:length(Ny_pos)
            Antenna_pos=[Antenna_pos; Nx_pos(nx), Ny_pos(ny), 0];
        end
    end
end

