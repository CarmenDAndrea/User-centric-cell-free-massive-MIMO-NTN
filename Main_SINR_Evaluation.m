clear all
close all
clc
%% SIMULATION PARAMETERS%%
M=32; %number of subcarriers
N=16; %number of OTFS symbols
f_star=20e9; %frequency at which the data are generated
fCarrier=20e9; % carrier frequency
lambda=3e8/fCarrier; %wavelength
n=3; % n of the subcarrier spacing, n=0: 15 kHz, n=1: 30 kHz, n=2: 60 kHz, n=3: 120 kHz
Delta_f=15e3*2^n; %subcarrier spacing
T=1/Delta_f; %symbol duration
B=Delta_f*M; % bandwidth
F = 5; %Downlink receiver Noise figure in dB
N0 = -174;  %PSD noise in dBm/Hz
noise_variance=B*10^(0.1*F)*10^(0.1*N0)*10^-3; % F*N0*B
% OTFS pulses
definePulses = 0;    % 1: Pulses defined. 0: Default rectangular pulses
rect = @(fx) 1.*(fx<1).*(fx>0);
pulseTx = rect;
pulseRx = rect;

% Array antenna design
Nx=16; %number of antennas on the x axis on the satellite
Ny=16; %number of antennas on the y axis on the satellite
NS=Nx*Ny; %number of antennas at the satellites
Nx_UT=4; %number of antennas on the x axis on the UT
Ny_UT=4; %number of antennas on the y axis on the UT
N_UT=Nx_UT*Ny_UT; %number of antennas at the UT
d=0.5*lambda; % antenna spacing
Gamma_dBi_SAT=23; %dBi
Gamma_dBi_UT=10; %dBi

[Antenna_pos_SAT] = Antenna_pos_generation(Nx,Ny,d);

[Antenna_pos_UT] = Antenna_pos_generation(Nx_UT,Ny_UT,d);

Pt_SAT_dBW=15; %total power budget on the satellites

N_UC=3; %values for the number of satellites serving each user

Num_Deployments=10; %number of deployments for the Starlink constellation

index_UT_interest=1; %UT for which the performance is evaluated

% Constraints for the association
N_max_UTs_per_satellite=10; %maximum number of UTs served from each satellite
min_Angle_SAT=2*pi/9; %minimum angle between users served from the same satellite

numChIt=100; %100 %number of channel iterations per snapshot

Beamformer_SAT='ZF'; % UPA or ZF Beamformer for the satellite

Combiner_UT='ZF'; % UPA or ZF Combiner for the UT

% 0: no constraints, 
% 1: only on satellite load, 
% 2: only on minimum angle, 
% 3: both satellite load and minimum angle
Association_constraint=2; 

name_base='Results_SINR_Starlink_';

filename=[name_base,'NRX_antennas_',num2str(N_UT),'_SAT_',Beamformer_SAT, '_UT_',Combiner_UT, '_Constraint',num2str(Association_constraint)],

data_path='./data_starlink/';


SINR_OTFS_per_user_symbol=cell(Num_Deployments,numChIt,N*M);
SINR_OFDM_per_user_symbol=cell(Num_Deployments,numChIt,N*M);


for ch=0:Num_Deployments-1
    %% Channel Generation (load from dati_starlink)
    load([data_path,'SlantRange_000',num2str(ch)]);
    load([data_path,'Latency_000',num2str(ch)]);
    load([data_path,'DoppShift_000',num2str(ch)]);
    load([data_path,'elevationSAT_000',num2str(ch)]); %angle at which UTs see the satellites
    load([data_path,'elevationUT_000',num2str(ch)]); %angle at which satellites see the UTs
    load([data_path,'azimuthSAT_000',num2str(ch)]); %angle at which UTs see the satellites
    load([data_path,'azimuthUT_000',num2str(ch)]); %angle at which satellites see the UTs

    Num_Total_satellites=size(srSATUT,1); %number of total satellites in the constellation

    % Generation of p681 LMS Channel
    [channLMS,totDopplerFreq] = p681LMS_Channel(fCarrier,fShift,numChIt,elSAT);

    % Remove NaN and reduce dimension of data
    [index_satellites,sr,channelLMS,del,fs,elevationSAT,azimuthSAT,elevationUT,azimuthUT] = Remove_NaN_from_Data_p681(srSATUT,channLMS,delay,totDopplerFreq,elSAT,azSAT,elUT,azUT);

    index_satellites_eff=index_satellites(2:end); %remove the first UT whihc is the ideal one

    % Data Elaboration
    [paths_energy,Delays,Doppler_frequencies,theta_angles_SAT_id, phi_angles_SAT_id,theta_angles_SAT_eff, phi_angles_SAT_eff,...
        theta_angles_UT_id, phi_angles_UT_id,theta_angles_UT_eff, phi_angles_UT_eff] = Data_Elaboration_p681(sr,channelLMS,del,fs,elevationSAT,azimuthSAT,elevationUT,azimuthUT,fCarrier,f_star);

    % Antenna Arrays evaluation

    [A_matrix_SAT,A_matrix_UT] = Array_response_Evaluation_SAT_UT(Antenna_pos_SAT, Antenna_pos_UT,theta_angles_SAT_eff,...
        phi_angles_SAT_eff,theta_angles_UT_eff,phi_angles_UT_eff,lambda, Gamma_dBi_SAT, Gamma_dBi_UT);

    sr_eff=sr(2:end);

    % Association rules
    if (Association_constraint==0)
        % No constraints
        [Alpha_association] = User_association_User_centric_distances(sr_eff,N_UC);
        [Satellites_Load] = Evaluate_Satellites_Load(Alpha_association,index_satellites_eff,Num_Total_satellites);
    elseif (Association_constraint==1)
        % Constraint on the Satellite Load
        [Alpha_association,Satellites_Load] = User_association_User_centric_distances_Satellites_Load(sr_eff,N_UC,N_max_UTs_per_satellite,Num_Total_satellites,index_satellites_eff);
    elseif (Association_constraint==2)
        % Constraint on minimum angle
        [Alpha_association, Satellites_Load] = User_association_User_centric_Min_Angle(sr_eff,N_UC,Num_Total_satellites,index_satellites_eff,theta_angles_UT_eff,phi_angles_UT_eff,min_Angle_SAT);
    elseif (Association_constraint==3)
        % Constraint on the satellite Load and minimum angle
        [Alpha_association, Satellites_Load] = User_association_User_centric_Satellites_Load_Min_Angle(sr_eff,N_UC,N_max_UTs_per_satellite,Num_Total_satellites,index_satellites_eff,theta_angles_UT_eff,phi_angles_UT_eff,min_Angle_SAT);
    end
    % Power allocation (Uniform)
    [Powers_DL] = Uniform_power_Allocation(Pt_SAT_dBW,Alpha_association,M,N,index_satellites_eff,Num_Total_satellites);
    % Beamformers and combiners (ideal positions)
    if strcmp(Beamformer_SAT,'UPA')
        [W] = Beamformers_UPA(index_UT_interest,theta_angles_UT_id, phi_angles_UT_id,theta_angles_UT_eff, phi_angles_UT_eff,...
            Antenna_pos_SAT,Gamma_dBi_SAT,lambda,Alpha_association);
    elseif strcmp(Beamformer_SAT,'ZF')
        [W] = Beamformers_ZF_selective(index_UT_interest,theta_angles_UT_id, phi_angles_UT_id,theta_angles_UT_eff, phi_angles_UT_eff,...
            Antenna_pos_SAT,Gamma_dBi_SAT,lambda,Alpha_association,index_satellites_eff);
    end

    if strcmp(Combiner_UT,'UPA')
        [V_OFDM] = Combiners_UPA(index_UT_interest, theta_angles_SAT_id, phi_angles_SAT_id,theta_angles_SAT_eff, phi_angles_SAT_eff,...
            Antenna_pos_UT,Gamma_dBi_UT,lambda,Alpha_association);
        [V_OTFS] = Combiners_UPA(index_UT_interest,theta_angles_SAT_id, phi_angles_SAT_id,theta_angles_SAT_eff, phi_angles_SAT_eff,...
            Antenna_pos_UT,Gamma_dBi_UT,lambda,Alpha_association);
    elseif strcmp(Combiner_UT,'ZF')
        [V_OFDM] = Combiner_ZF(index_UT_interest,theta_angles_SAT_id, phi_angles_SAT_id,theta_angles_SAT_eff, phi_angles_SAT_eff,...
            Antenna_pos_UT,Gamma_dBi_UT,lambda,Alpha_association);
        [V_OTFS] = Combiner_ZF(index_UT_interest, theta_angles_SAT_id, phi_angles_SAT_id,theta_angles_SAT_eff, phi_angles_SAT_eff,...
            Antenna_pos_UT,Gamma_dBi_UT,lambda,Alpha_association);
    end


    % Processings

    % OTFS e OFDM processings

    Psi_OTFS = Psi_matrix_OTFS_reduced_size(index_UT_interest,N,M,Delays, Doppler_frequencies,Delta_f,T,definePulses,pulseTx,pulseRx);

    Psi_OFDM = Psi_matrix_OFDM_reduced_size(index_UT_interest,N,M,Delays, Doppler_frequencies,Delta_f,T);

    % Delta_matrix_definitions

    P_index_UT_interest=size(Delays{index_UT_interest},1);

    Delta_OTFS=cell(P_index_UT_interest,1);

    Delta_OFDM=cell(P_index_UT_interest,1);
% No phase compensation
    for pp=1:P_index_UT_interest
        if Alpha_association{index_UT_interest}(pp)==1
            Delta_OTFS{pp,1}=eye(N*M);
            Delta_OFDM{pp,1}=eye(N*M);
        else
            Delta_OTFS{pp}=zeros(N*M);
            Delta_OFDM{pp}=zeros(N*M);
        end
    end


    for chIT=1:numChIt

        [Rho, Rho_tilde] = Rho_Rho_Tilde_Evaluation_p681(paths_energy,Delays,chIT);

        % LMMSE Post Coding

        [D_matrices_OTFS,N_matrices_OTFS] = Practical_LMMSE_Detection_MultiAntenna_UTs(index_UT_interest,Alpha_association, Delta_OTFS,Psi_OTFS, A_matrix_SAT, A_matrix_UT, W, V_OTFS, Rho_tilde,Powers_DL, noise_variance);

        [D_matrices_OFDM,N_matrices_OFDM] = Practical_LMMSE_Detection_MultiAntenna_UTs(index_UT_interest,Alpha_association, Delta_OFDM, Psi_OFDM, A_matrix_SAT, A_matrix_UT, W, V_OFDM, Rho,Powers_DL, noise_variance);
        
        % SINR Evaluation per user and symbol

        SINR_OTFS_mat=SINR_Evaluation_reduced_size_Psi(index_UT_interest,Alpha_association,Delta_OTFS,Psi_OTFS, A_matrix_SAT, A_matrix_UT, W, V_OTFS, Rho_tilde, D_matrices_OTFS,N_matrices_OTFS,Powers_DL, noise_variance,index_satellites_eff);

        SINR_OTFS_per_user_symbol(ch+1,chIT,:)=SINR_OTFS_mat(:);

        SINR_OFDM_mat=SINR_Evaluation_reduced_size_Psi(index_UT_interest,Alpha_association,Delta_OFDM,Psi_OFDM, A_matrix_SAT, A_matrix_UT, W, V_OFDM, Rho, D_matrices_OFDM,N_matrices_OFDM, Powers_DL, noise_variance,index_satellites_eff);

        SINR_OFDM_per_user_symbol(ch+1,chIT,:)= SINR_OFDM_mat(:);

        chIT,
    end


    ch,

end

save(filename)
