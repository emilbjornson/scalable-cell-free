function [gainOverNoisedB,R,pilotIndex,D,APpositions,UEpositions,distances] = generateSetup(L,K,N,tau_p,nbrOfSetups)
%This function generates realizations of the simulation setup described in
%Section VI.
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
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%
%OUTPUT:
%gainOverNoisedB = Matrix with dimension L x K x nbrOfSetups where
%                  element (l,k,n) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k,n) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilot assigned to the users
%D               = Matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is zero AP l serves UE k in setup n and zero otherwise
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distances       = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs


%Size of the coverage area (as a square with wrap-around)
squareLength = 2000; %meter


%% Define propagation model

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss exponent
alpha = 3.76;

%Standard deviation of shadow fading
sigma_sf = 10;

%Average channel gain in dB at a reference distance of 1 meter. Note that
%-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Angular standard deviation around the nominal angle (measured in degrees)
ASDdeg = 20;

%Set threshold for when a non-master AP decides to serve a UE
threshold = -40; %dB


%Prepare to save results
gainOverNoisedB = zeros(L,K,nbrOfSetups);
R = zeros(N,N,L,K,nbrOfSetups);
distances = zeros(L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
D = zeros(L,K,nbrOfSetups);
masterAPs = zeros(K,1);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Random AP locations with uniform distribution
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    %Random UE locations with uniform distribution
    UEpositions = (rand(K,1) + 1i*rand(K,1)) * squareLength;
    
    
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    
    
    
    for k = 1:K
        
        %Compute distances assuming that the APs are 10 m above the UEs
        [distancetoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
        distances(:,k,n) = sqrt(10^2+distancetoUE.^2);
        
        %Compute the channel gain divided by noise power
        gainOverNoisedB(:,k,n) = constantTerm - alpha*10*log10(distances(:,k,n)) + sigma_sf*randn(size(distances(:,k,n))) - noiseVariancedBm;

        %Determine the master AP for UE k by looking for AP with best
        %channel condition
        [~,master] = max(gainOverNoisedB(:,k,n));
        D(master,k,n) = 1;
        masterAPs(k) = master;
        
        %Assign orthogonal pilots to the first tau_p UEs
        if k <= tau_p
            
            pilotIndex(k,n) = k;
            
        else %Assign pilot for remaining UEs
            
            %Compute received power from to the master AP from each pilot
            pilotinterference = zeros(tau_p,1);
            
            for t = 1:tau_p
                
                pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1,n)==t,n)));
                
            end
            
            %Find the pilot with the least receiver power
            [~,bestpilot] = min(pilotinterference);
            pilotIndex(k,n) = bestpilot;
            
        end
        
        

        %Go through all APs
        for l = 1:L
            
            %Compute nominal angle between UE k and AP l
            angletoUE = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l)));
            
            %Generate normalized spatial correlation matrix using the local
            %scattering model
            R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);
            
        end
        
    end
    
    
    %Each AP serves the UE with the strongest channel condition on each of
    %the pilots where the AP isn't the master AP, but only if its channel
    %is not too weak compared to the master AP
    for l = 1:L
        
        for t = 1:tau_p
            
            pilotUEs = find(t==pilotIndex(:,n));
            
            if sum(D(l,pilotUEs,n)) == 0 %If the AP is not a master AP
                
                %Find the UE with pilot t with the best channel
                [gainValue,UEindex] = max(gainOverNoisedB(l,pilotUEs,n));
                
                %Serve this UE if the channel is at most "threshold" weaker
                %than the master AP's channel
                %[gainValue gainOverNoisedB(masterAPs(pilotUEs(UEindex)),pilotUEs(UEindex),n)]   
                if gainValue - gainOverNoisedB(masterAPs(pilotUEs(UEindex)),pilotUEs(UEindex),n) >= threshold
                    
                    D(l,pilotUEs(UEindex),n) = 1;
                    
                end
                
            end
            
        end
        
    end
    
end
