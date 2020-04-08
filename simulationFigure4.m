%This Matlab script generates Figure 4 in the paper:
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


%Number of channel realizations
nbrOfRealizations = 10000000;

%Generate CN(0,1) realizations
h = (randn(nbrOfRealizations,1) + 1i*randn(nbrOfRealizations,1))/sqrt(2);


%% Compute precoding vectors 
y = abs(h).^2;
normalization = mean(y./(y+1).^2);
w_MR = h/1; %MR
w_LP_MMSE = (h./(y+1))/sqrt(normalization); %LP-MMSE

%Compute effective channel gains (the imaginary part is zero in theory)
gain_MR = real(conj(h).*w_MR);
gain_LP_MMSE = real(conj(h).*w_LP_MMSE);


%% Plot the simulation results
figure; box on;
h1 = histogram(gain_MR);
hold on
h2 = histogram(gain_LP_MMSE);
h1.Normalization = 'pdf';
h1.BinWidth = 0.05;
h2.Normalization = 'pdf';
h2.BinWidth = 0.05;
xlabel('Channel gain','Interpreter','Latex');
ylabel('PDF','Interpreter','Latex');
xlim([0,5])
plot([0.25 0.6],[0.8 0.86],'k-');
text(0.65,0.86,'MR','Interpreter','Latex');
plot([1.3 1.65],[0.62 0.68],'k-');
text(1.7,0.68,'LP-MMSE','Interpreter','Latex');
