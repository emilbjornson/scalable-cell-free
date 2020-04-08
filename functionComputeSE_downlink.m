function [SE_MR,SE_LP_MMSE,SE_MR_perfect,SE_LP_MMSE_perfect,SE_P_MMSE,SE_P_MMSE_perfect] = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,rho_dist,R,pilotIndex,rho_central)
%Compute DL SE for different precoding schemes using the hardening bound in
%Proposition 3.
%
%This function was developed as a part of the paper:
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
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = Matrix with dimension L x K where (l,k) is one if AP l
%                    serves UE k and zero otherwise
%B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the estimate between
%                    AP l and UE k in setup n, normalized by noise
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k in setup n,
%                    normalized by noise
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per AP
%K                 = Number of UEs
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%rho_dist          = Matrix with dimension L x K where (l,k) is the
%                    downlink transmit power that AP l allocates to UE k
%                    when using distributed precoding
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k
%                    in setup n, normalized by noise
%pilotIndex        = Vector containing the pilot assigned to each UE
%rho_central       = Vector containing the total power allocated to each UE
%                    when using centralized precoding
%
%OUTPUT:
%SE_MR              = Vector with the SE achieved with MR precoding
%SE_LP_MMSE         = Same as SE_MR but with LP-MMSE precoding
%SE_MR_perfect      = Same as SE_MR but with perfect CSI at the UEs
%SE_LP_MMSE_perfect = Same as SE_LP_MMSE but with perfect CSI at the UEs
%SE_P_MMSE          = Same as SE_MR but with P-MMSE precoding
%SE_P_MMSE_perfect  = Same as SE_P_MMSE but with perfect CSI at the UEs



%Store the N x N identity matrix
eyeN = eye(N);

%Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_p/tau_c);


%Prepare to store simulation results
signal_MR = zeros(K,1);
interf_MR = zeros(K,1);
cont_MR  = zeros(K,K);
scaling_MR = zeros(L,K);
interUserGains_MR = zeros(K,K,nbrOfRealizations);

if nargout > 1
    signal_LP_MMSE = zeros(K,1);
    interf_LP_MMSE = zeros(K,1);
    scaling_LP_MMSE = zeros(L,K);
    interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);
end


if nargout > 4
    signal_P_MMSE = zeros(K,1);
    interf_P_MMSE = zeros(K,1);
    scaling_P_MMSE = zeros(K,1);
    interUserGains_P_MMSE = zeros(K,K,nbrOfRealizations);
end



%% Compute scaling factors for precoding

%Computation for MR precoding
for l = 1:L
    
    %Extract which UEs are served by AP l
    servedUEs = find(D(l,:)==1);
    
    for ind = 1:length(servedUEs)
        
        %Compute scaling factor using the spatial correlation matrix of the
        %channel estimate
        scaling_MR(l,servedUEs(ind)) = trace(B(:,:,l,servedUEs(ind)));
        
    end
    
end


%Computation for LP-MMSE and P-MMSE precoding
if nargout > 1
    
    for n = 1:nbrOfRealizations
        
        %LP-MMSE precoding
        for l = 1:L
            
            %Extract which UEs are served by the AP
            servedUEs = find(D(l,:)==1);
            
            %Compute sum of estimation error covariance matrices of the UEs
            %served by AP l
            Cserved = sum(C(:,:,l,servedUEs),4);
            
            %Compute LP-MMSE precoding
            V_MR = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
            V_LP_MMSE = p*(p*(V_MR(:,servedUEs)*V_MR(:,servedUEs)'+Cserved)+eyeN)\V_MR(:,servedUEs);
            
            %Compute scaling factor by Monte Carlo methods
            scaling_LP_MMSE(l,servedUEs) = scaling_LP_MMSE(l,servedUEs) + sum(abs(V_LP_MMSE).^2,1)/nbrOfRealizations;
            
        end
        
        %P-MMSE precoding
        if nargout > 4
            
            %Go through all UEs
            for k = 1:K
                
                %Determine the set of serving APs
                servingAPs = find(D(:,k)==1);
                La = length(servingAPs);
                
                %Determine which UEs that are served by partially the same set
                %of APs as UE k
                servedUEs = sum(D(servingAPs,:),1)>=1;
                
                
                %Extract channel realizations and estimation error correlation
                %matrices for the APs that involved in the service of UE k
                Hhatallj_active = zeros(N*La,K);
                C_tot_blk = zeros(N*La,N*La);
                C_tot_blk_partial = zeros(N*La,N*La);
                
                for l = 1:La
                    Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
                    C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
                    C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
                end
                
                %Compute P-MMSE precoding
                V_P_MMSE = p*(p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k);
                
                %Compute scaling factor by Monte Carlo methods
                scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;
                
            end
            
        end
        
    end
    
end


%Compute all terms in Corollary 3 for MR precoding
for l = 1:L
    
    %Extract which UEs are served by the AP
    servedUEs = find(D(l,:)==1);
    
    %Go through all UEs served by the AP
    for ind = 1:length(servedUEs)
        
        %Extract UE index
        k = servedUEs(ind);
        
        %Desired signal term
        signal_MR(k) = signal_MR(k) + sqrt(rho_dist(l,k)*real(trace(B(:,:,l,k))));
        
        for i = 1:K
            
            %Non-coherent interference from UE k to UE i
            interf_MR(i) = interf_MR(i) + rho_dist(l,k)*real(trace(B(:,:,l,k)*R(:,:,l,i)))/real(trace(B(:,:,l,k)));
            
            if pilotIndex(k) == pilotIndex(i)
                
                %Coherent interference from UE k to UE i
                cont_MR(i,k) = cont_MR(i,k) + sqrt(rho_dist(l,k))*real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i)))/sqrt(real(trace(B(:,:,l,k))));
                
            end
            
        end
        
    end
    
end


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    %Matrix to store Monte-Carlo results for this realization
    if nargout > 1
        interf_LP_MMSE_n = zeros(K,K);
    end
    
    if nargout > 4
        interf_P_MMSE_n = zeros(K,K);
    end
    
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract channel estimates from all UEs to AP l
        Hhatallj = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract which UEs are served by the AP
        servedUEs = find(D(l,:)==1);
        
        %Compute MR precoding
        V_MR = Hhatallj(:,servedUEs);
        
        if nargout > 1 %Compute LP-MMSE precoding
            
            %Compute sum of estimation error covariance matrices of the UEs
            %served by AP l
            Cserved = sum(C(:,:,l,servedUEs),4);
            
            V_LP_MMSE = p*(p*(V_MR*V_MR'+Cserved)+eyeN)\V_MR;
            
        end
        
        
        %Go through all UEs served by the AP
        for ind = 1:length(servedUEs)
            
            %Extract UE index
            k = servedUEs(ind);
            
            %%Normalize MR precoding
            w = V_MR(:,ind)*sqrt(rho_dist(l,k)/scaling_MR(l,k));
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_MR(:,k,n) = interUserGains_MR(:,k,n) + Hallj'*w;
            
            
            %%LP-MMSE precoding
            if nargout > 1
                
                %Normalize LP-MMSE precoding
                w = V_LP_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_LP_MMSE(l,k));
                
                %Compute realizations of the terms inside the expectations
                %of the signal and interference terms
                signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
                interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + Hallj'*w;
                
                %Compute gain of the signal from UE that arrives at other UEs
                interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*w;
                
            end
            
        end
        
    end
    
    
    %Consider the centralized P-MMSE precoding scheme
    if nargout > 4
        
                %Go through all UEs
        for k = 1:K
            
                        %Determine the set of serving APs
            servingAPs = find(D(:,k)==1);
                        La = length(servingAPs);
                        
            %Determine which UEs that are served by partially the same set
            %of APs as UE k
            servedUEs = sum(D(servingAPs,:),1)>=1;
            
            %Extract channel realizations and estimation error correlation
            %matrices for the APs that involved in the service of UE k
            Hhatallj_active = zeros(N*La,K);
            C_tot_blk = zeros(N*La,N*La);
            C_tot_blk_partial = zeros(N*La,N*La);
            
            for l = 1:La
                Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
                C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
                C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
            end
            
            %Compute P-MMSE precoding
            w = p*(p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k);
            
            %Apply power allocation
            w = w*sqrt(rho_central(k)/scaling_P_MMSE(k));
            
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms
            signal_P_MMSE(k) = signal_P_MMSE(k) + (Hhatallj_active(:,k)'*w)/nbrOfRealizations;
            interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + Hhatallj_active'*w;
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hhatallj_active'*w;
            
        end
        
    end
    
    
    %Compute interference power in one realization
    if nargout > 1
        interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,2)/nbrOfRealizations;
    end
    
    if nargout > 1
        interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;
    end
    
end


%% Compute the SEs with different precoding schemes

%Compute SE with MR using the closed-form expressions in Corollary 3
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR).^2) ./ (interf_MR + sum(abs(cont_MR).^2,2) - abs(signal_MR).^2 + 1)));

%Compute SE with LP-MMSE
if nargout > 1
    SE_LP_MMSE = prelogFactor*real(log2(1+(abs(signal_LP_MMSE).^2) ./ (interf_LP_MMSE - abs(signal_LP_MMSE).^2 + 1)));
end

%Compute SE with P-MMSE
if nargout > 1
    SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - abs(signal_P_MMSE).^2 + 1)));
end


%Prepare to compute SEs with perfect CSI at the UEs
if nargout > 2
    SE_MR_perfect = zeros(K,1);
end

if nargout > 3
    SE_LP_MMSE_perfect = zeros(K,1);
end

if nargout > 5
    SE_P_MMSE_perfect = zeros(K,1);
end

for k = 1:K
    
    %Compute SE with MR precoding, assuming perfect CSI at the UE
    if nargout > 2
        
        SE_MR_perfect(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_MR(k,k,:)).^2 ./ ( sum(abs(interUserGains_MR(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    end
    
    %Compute SE with LP-MMSE precoding, assuming perfect CSI at the UE
    if nargout > 3
        
        SE_LP_MMSE_perfect(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_LP_MMSE(k,k,:)).^2 ./ ( sum(abs(interUserGains_LP_MMSE(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    end
    
    %Compute SE with P-MMSE precoding, assuming perfect CSI at the UE
    if nargout > 5
        
        SE_P_MMSE_perfect(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_P_MMSE(k,k,:)).^2 ./ ( sum(abs(interUserGains_P_MMSE(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
    end
    
end
