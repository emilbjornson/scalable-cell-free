function [SE_MR,SE_LP_MMSE,SE_P_MMSE,SE_MMSE] = functionComputeSE_uplink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex)
%Compute UL SE for different receive combining schemes using the capacity
%bound in Proposition 1 for the centralized schemes and the capacity bound
%in Propostion 2 for the distributed schemes.
%
%This function was developed as a part of the paper:
%
%Emil Bjornson, Luca Sanguinetti, “Scalable Cell-Free Massive MIMO
%Systems,” IEEE Transactions on Communications, vol. 68, no. 7, pp.
%4247-4261, July 2020.
%
%Download article: http://arxiv.org/pdf/1908.03119
%
%This is version 1.1 (Last edited: 2020-12-06)
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
%N                 = Number of antennas per AP
%K                 = Number of UEs per cell
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k
%                    in setup n, normalized by noise
%pilotIndex        = Vector containing the pilot assigned to each UE
%
%OUTPUT:
%SE_MR              = Vector with the SE achieved with MR combining
%SE_LP_MMSE         = Same as SE_MR but with LP-MMSE combining
%SE_P_MMSE          = Same as SE_MR but with P-MMSE combining
%SE_MMSE            = Same as SE_MR but with MMSE combining


%Store the N x N identity matrix
eyeN = eye(N);

%Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_p/tau_c);



%Prepare to store simulation results
signal_MR = zeros(K,1);
interf_MR = zeros(K,1);
cont_MR  = zeros(K,K);
scaling_MR  = zeros(K,L);


if nargout > 1
    signal_LP_MMSE = zeros(K,1);
    interf_LP_MMSE = zeros(K,1);
    scaling_LP_MMSE = zeros(L,K,nbrOfRealizations);
    interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);
end

if nargout > 2
    SE_P_MMSE = zeros(K,1);
end

if nargout > 3
    SE_MMSE = zeros(K,1);
end



%% Compute scaling factors for the combining vectors

%Compute all terms in Corollary 2 for MR combining

for l = 1:L
    %Extract which UEs are served by the AP
    servedUEs = find(D(l,:)==1);
    
    %Go through all UEs served by the AP
    for ind = 1:length(servedUEs)
        
        %Extract UE index
        k = servedUEs(ind);
        
        
        %Noise scaling
        scaling_MR(k,l) = trace(B(:,:,l,k));

        %Desired signal term
        signal_MR(k) = signal_MR(k) + sqrt(p)*real(trace(B(:,:,l,k)));
        
        for i = 1:K
            
            %Non-coherent interference from UE i to UE k
            interf_MR(k) = interf_MR(k) + p*real(trace(B(:,:,l,k)*R(:,:,l,i))); %/real(trace(B(:,:,l,k)));
            
            if pilotIndex(k) == pilotIndex(i)
                
                %Coherent interference from UE i to UE k
                cont_MR(i,k) = cont_MR(i,k) + sqrt(p)*real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i))); %/sqrt(real(trace(B(:,:,l,k))));
                
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

    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract channel estimates from all UEs to AP l
        Hhatallj = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute sum of estimation error covariance matrices of the UEs
        %served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);
        
        %Compute MR combining
        V_MR = Hhatallj(:,servedUEs);
        
        if nargout > 1 %Compute LP-MMSE combining
            V_LP_MMSE = p*(p*(V_MR*V_MR'+Cserved)+eyeN)\V_MR;
        end
        
        
        %Go through all UEs served by the AP
        for ind = 1:length(servedUEs)
            
            %Extract UE index
            k = servedUEs(ind);
            
            
            %%LP-MMSE combining
            if nargout > 1
                
                %Compute LP-MMSE combining
                v = V_LP_MMSE(:,ind);
                
                %Compute realizations of the terms inside the expectations
                %of the signal and interference terms
                signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*v)/nbrOfRealizations;
                interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + Hallj'*v;
                scaling_LP_MMSE(l,k,n) = norm(v).^2;
                
                %Compute gain of the signal from UE that arrives at other UEs
                interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*v;
                
            end
            
        end
        
    end
    
    
    %Compute interference power in one realization
    if nargout > 1
        interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,1)'/nbrOfRealizations;
    end
    
    
    %Consider the centralized schemes
    if nargout > 2
        
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
            
            %Compute P-MMSE combining
            v = p*(p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)'+C_tot_blk_partial)+eye(La*N))\Hhatallj_active(:,k);

            %Compute numerator and denominator of instantaneous SINR
            numerator = p*abs(v'*Hhatallj_active(:,k))^2;
            denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_P_MMSE(k) = SE_P_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            
            
            %Compute MMSE combining
            v = p*(p*(Hhatallj_active*Hhatallj_active'+C_tot_blk)+eye(La*N))\Hhatallj_active(:,k);

            %Compute numerator and denominator of instantaneous SINR
            numerator = p*abs(v'*Hhatallj_active(:,k))^2;
            denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MMSE(k) = SE_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            
        end
        
    end
        
end


%% Compute the SEs with different combining

%Compute SE with MR using the closed-form expression in Corollary 2
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR).^2) ./ (interf_MR + sum(abs(cont_MR).^2,1)' - abs(signal_MR).^2 + sum(scaling_MR,2))));

%Compute SE with LP-MMSE using Monte Carlo method
if nargout > 1
    SE_LP_MMSE = prelogFactor*real(log2(1+(abs(signal_LP_MMSE).^2) ./ (interf_LP_MMSE - abs(signal_LP_MMSE).^2 + sum(mean(scaling_LP_MMSE,3),1)'/p )));
end
