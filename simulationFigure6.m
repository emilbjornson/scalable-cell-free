%This Matlab script generates Figure 6(a) and Figure 6(b) in the paper:
%
%Emil Bjornson, Luca Sanguinetti, “Scalable Cell-Free Massive MIMO
%Systems,” IEEE Transactions on Communications, vol. 68, no. 7, pp.
%4247-4261, July 2020.
%
%Download article: http://arxiv.org/pdf/1908.03119
%
%This is version 1.01 (Last edited: 2020-12-06)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


close all;
clear;

%Select which figure to generate:
%selectSimulationSetup = 1 gives Figure 6(a)
%selectSimulationSetup = 2 gives Figure 6(b)
selectSimulationSetup = 2;


%% Define simulation setup

%Number of Monte-Carlo setups
nbrOfSetups = 25;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

if selectSimulationSetup == 1
    
    %Number of APs per setup
    L = 400;
    
    %Number of antennas per AP
    N = 1;
    
elseif selectSimulationSetup == 2
    
    %Number of APs per setup
    L = 100;
    
    %Number of antennas per AP
    N = 4;
    
end

%Number of UEs in the network
K = 100;

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 10;


%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 1000;

%Prepare to save simulation results
SE_MR_tot = zeros(K,nbrOfSetups);
SE_LP_MMSE_tot = zeros(K,nbrOfSetups);
SE_P_MMSE_tot = zeros(K,nbrOfSetups);
SE_MR_perfect_tot = zeros(K,nbrOfSetups);
SE_LP_MMSE_perfect_tot = zeros(K,nbrOfSetups);
SE_P_MMSE_perfect_tot = zeros(K,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    [gainOverNoisedB,R,pilotIndex,D] = generateSetup(L,K,N,tau_p,1);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    %Compute the equal power allocation for centralized precoding
    rho_central = (rho_tot/tau_p)*ones(K,1);
    
    %Compute the power allocation in Eq. (43) for distributed precoding
    rho_dist = zeros(L,K);
    gainOverNoise = db2pow(gainOverNoisedB);
    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in Eq. (43)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            
        end
        
    end
    
    
    %Compute SE using the capacity bound in Proposition 3
    [SE_MR,SE_LP_MMSE,SE_MR_perfect,SE_LP_MMSE_perfect,SE_P_MMSE,SE_P_MMSE_perfect] = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,rho_dist,R,pilotIndex,rho_central);
    
    %Save the SE values
    SE_MR_tot(:,n) = SE_MR;
    SE_LP_MMSE_tot(:,n) = SE_LP_MMSE;
    SE_P_MMSE_tot(:,n) = SE_P_MMSE;
    SE_MR_perfect_tot(:,n) = SE_MR_perfect;
    SE_LP_MMSE_perfect_tot(:,n) = SE_LP_MMSE_perfect;
    SE_P_MMSE_perfect_tot(:,n) = SE_P_MMSE_perfect;
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end


%% Plot simulation results
figure;
hold on; box on;
plot(sort(SE_P_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_LP_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_MR_tot(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'P-MMSE (Scalable)','LP-MMSE (Scalable)','MR (Scalable)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 10]);
