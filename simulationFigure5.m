%This Matlab script generates Figure 5(a) and Figure 5(b) in the paper:
%
%Emil Bjornson, Luca Sanguinetti, “Scalable Cell-Free Massive MIMO
%Systems,” IEEE Transactions on Communications, to appear. 
%
%Download article: http://arxiv.org/pdf/1908.03119
%
%This is version 1.0 (Last edited: 2020-04-08)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


close all;
clear;

%Select which figure to generate:
%selectSimulationSetup = 1 gives Figure 5(a)
%selectSimulationSetup = 2 gives Figure 5(b)
selectSimulationSetup = 1;


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

%Prepare to save simulation results
SE_scalable_MR_tot = zeros(K,nbrOfSetups);
SE_scalable_LP_MMSE_tot = zeros(K,nbrOfSetups);
SE_scalable_P_MMSE_tot = zeros(K,nbrOfSetups);
SE_scalable_MMSE_tot = zeros(K,nbrOfSetups);
SE_all_MR_tot = zeros(K,nbrOfSetups);
SE_all_LP_MMSE_tot = zeros(K,nbrOfSetups);
SE_all_P_MMSE_tot = zeros(K,nbrOfSetups);
SE_all_MMSE_tot = zeros(K,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    [gainOverNoisedB,R,pilotIndex,D] = generateSetup(L,K,N,tau_p,1);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    
    %% Proposed Scalable Cell-Free Massive MIMO
    
    %Compute SE using Propositions 1 and 2
    [SE_MR,SE_LP_MMSE,SE_P_MMSE,SE_MMSE] = functionComputeSE_uplink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);
    
    %Save SE values
    SE_scalable_MR_tot(:,n) = SE_MR;
    SE_scalable_LP_MMSE_tot(:,n) = SE_LP_MMSE;
    SE_scalable_P_MMSE_tot(:,n) = SE_P_MMSE;
    SE_scalable_MMSE_tot(:,n) = SE_MMSE;
    
    
    %% Original Cell-Free Massive MIMO
    
    %Define the case when all APs serve all UEs
    D_all = ones(L,K);
    
    %Compute SE using Propositions 1 and 2
    [SE_MR,SE_LP_MMSE,SE_P_MMSE,SE_MMSE] = functionComputeSE_uplink(Hhat,H,D_all,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);
    
    %Save SE values
    SE_all_MR_tot(:,n) = SE_MR;
    SE_all_LP_MMSE_tot(:,n) = SE_LP_MMSE;
    SE_all_P_MMSE_tot(:,n) = SE_P_MMSE;
    SE_all_MMSE_tot(:,n) = SE_MMSE;
    
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end


%% Plot simulation results
figure;
hold on; box on;

plot(sort(SE_all_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_scalable_P_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_all_LP_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_scalable_LP_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_all_MR_tot(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',4);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'MMSE (All)','P-MMSE (Scalable)','L-MMSE (All)','LP-MMSE (Scalable)','MR (All)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 10]);
