%This Matlab script can be used to generate Figures 4-6 and 8-14 in
%the article:
%
%Emil Bjornson, Erik G. Larsson, Merouane Debbah, "Massive MIMO for Maximal
%Spectral Efficiency: How Many Users and Pilots Should Be Allocated?,"
%to appear in IEEE Transactions on Wireless Communications.
%
%Download article: http://arxiv.org/pdf/1412.7102
%
%This is version 1.0 (Last edited: 2015-10-06)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.


%Initialization
close all;
clear all;


%%Select simulation case:
%
%simulationCase = 1: Figures 4-6 and 14
%simulationCase = 2: Figures 8-11
%simulationCase = 3: Figures 12
%simulationCase = 4: Figures 13
simulationCase = 1;


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

%Pathloss exponent
kappa = 3.7;

%Number of directions to look for interfering cells (for hexagonal cells)
directions = 6;

%Percentage of the radius inside the cell where no UEs are allowed
forbiddenRegion = 0.14;

%Parameters for the Monte Carlo simulations
monteCarloUEs = 1000000; %Number of random UE locations per cell

%Compute various combinations of the mu-parameters Eq. (18), using
%Monte Carlo simulations
[muValues1Mean,muValues2Mean,reuseMu1Mean,reuseMu1Mean2,reuseMu1MeanNext,reuseMu1Mean2Next,reuseMu2Mean,reuseMuMeanVariance,muValues1Worst,muValues2Worst,reuseMu1Worst,reuseMu1Worst2,reuseMu1WorstNext,reuseMu1Worst2Next,reuseMu2Worst,reuseMuWorstVariance,muValues1Best,muValues2Best,reuseMu1Best,reuseMu1Best2,reuseMu1BestNext,reuseMu1Best2Next,reuseMu2Best,reuseMuBestVariance,reuseFactor] = computeEnvironment(kappa,forbiddenRegion,monteCarloUEs);


if simulationCase == 1 %Simulations from Section IV.A and Section V (Figures 4-6 and 14)
    
    %Select range of BS antennas
    nbrOfMvalues = 1000; %Number of different cases
    Mvalues = round(logspace(1,5,nbrOfMvalues)); %Spread out antenna numbers equally in log-scale
    
    %Coherence block length
    S = 400 * ones(1,2);
    
    %Inverse SNR value
    sigma2rho = 1/10^(5/10) * ones(1,2); %5 dB
    
    %EVM value
    epsilon2 = [0 0.1^2];
    
elseif simulationCase == 2 %Some simulations from Section IV.B (Figures 8-11)
    
    %Select range of BS antennas
    Mvalues = 1:1000; %All antenna numbers from 1 to 1000
    
    %Coherence block length
    S = 400;
    
    %Inverse SNR value
    sigma2rho = 1/10^(5/10); %5 dB
    
    %EVM value
    epsilon2 = 0;
    
elseif simulationCase == 3 %One simulation from Section IV.B (Figure 12)
    
    %Select range of BS antennas
    Mvalues = [100 500];
    
    %Range of SNR values in dB
    SNRvaluesdB = -10:0.1:20;
    
    %Inverse SNR values in linear scale
    sigma2rho = 1./10.^(SNRvaluesdB/10);
    
    %Coherence block length
    S = 400 * ones(size(sigma2rho));
    
    %EVM value
    epsilon2 = zeros(size(sigma2rho));
    
elseif simulationCase == 4 %One simulation from Section IV.B (Figure 13)
    
    %Select range of BS antennas
    Mvalues = [100 500];
    
    %Coherence block length
    S = 10:10:2000;
    
    %Inverse SNR values in linear scale
    sigma2rho = 1/10^(5/10)*ones(size(S)); %SNR is 5 dB
    
    %EVM value
    epsilon2 = zeros(size(sigma2rho));
    
end


%Define the range of UEs to consider
Kvalues = 1:max(S);


%Compute the sum of all mu values in Eq. (18)
mu1all_mean = 1+directions*(sum(muValues1Mean(:))-1);
mu1all_worst = 1+directions*(sum(muValues1Worst(:))-1);
mu1all_best = 1+directions*(sum(muValues1Best(:))-1);


%Extract only reuse factors smaller or equal to 7
reuseIndices = find(reuseFactor>0 & reuseFactor<=directions+1);
for j = 1:length(reuseIndices);
    if sum(reuseFactor(reuseIndices(j))==reuseFactor(reuseIndices(1:j-1)))>0
        reuseIndices(j)=1;
    end
end
reuseIndices = reuseIndices(reuseIndices>1);



%%Compute spectral efficiencies according to Theorem 1 and 2.

%Placeholders for storing spectral efficiencies
SE_MR_mean = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));
SE_ZF_mean = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));
SE_PZF_mean = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));

SE_MR_worst = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));
SE_ZF_worst = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));
SE_PZF_worst = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));

SE_MR_best = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));
SE_ZF_best = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));
SE_PZF_best = zeros(length(Mvalues),max(S),length(reuseIndices),length(S));

SE_asymptotic_mean = zeros(max(S),length(reuseIndices),length(S));
SE_asymptotic_worst = zeros(max(S),length(reuseIndices),length(S));
SE_asymptotic_best = zeros(max(S),length(reuseIndices),length(S));


%Go through the different reuse factors
for j = 1:length(reuseIndices);
    
    %Extract the reuse factor
    currentReuseFactor = reuseFactor(reuseIndices(j));
    
    %Extract sum of mu-values for current reuse factor for mean interference
    mu1reuse_mean = directions*reuseMu1Mean(reuseIndices(j));
    mu2reuse_mean = directions*reuseMu2Mean(reuseIndices(j));
    variance_mean = directions*reuseMuMeanVariance(reuseIndices(j));
    
    %Extract sum of mu-values for current reuse factor for worst interference
    mu1reuse_worst = directions*reuseMu1Worst(reuseIndices(j));
    mu2reuse_worst = directions*reuseMu2Worst(reuseIndices(j));
    variance_worst = directions*reuseMuWorstVariance(reuseIndices(j));
    
    %Extract sum of mu-values for current reuse factor for best interference
    mu1reuse_best = directions*reuseMu1Best(reuseIndices(j));
    mu2reuse_best = directions*reuseMu2Best(reuseIndices(j));
    variance_best = directions*reuseMuBestVariance(reuseIndices(j));
    
    %Number of neighbors that use each of the other sets of pilots
    neighborsPerOtherPilot = directions/(currentReuseFactor-1);
    
    
    %Go through different number of antennas
    for n = 1:length(Mvalues)
        
        %Go through different other cases (varying SNR, coherence block, etc.)
        for m = 1:length(S)
            
            %Go through different number of UEs (limited by coherence block length)
            for K = 1:S(m)
                
                %Compute length of pilot signal
                B = currentReuseFactor*K;
                
                if B < S(m)
                    
                    
                    %Compute asymptotic limits according to Corollary 2
                    if n == 1
                        SE_asymptotic_mean(K,j,m) = K*(1-B/S(m))*log2(1+(1-epsilon2(m))/(mu2reuse_mean+epsilon2(m)));
                        SE_asymptotic_worst(K,j,m) = K*(1-B/S(m))*log2(1+(1-epsilon2(m))/(mu2reuse_worst+epsilon2(m)));
                        SE_asymptotic_best(K,j,m) = K*(1-B/S(m))*log2(1+(1-epsilon2(m))/(mu2reuse_best+epsilon2(m)));
                    end
                    
                    
                    %Maximum ratio (MR) combining/precoding
                    %
                    %Achievable spectral efficiency using the formula in
                    %Theorem 1, for mean, worst, and best case interference
                    SINR_MR_mean = B*(1-epsilon2(m))/(epsilon2(m)*B + (mu1all_mean*K + sigma2rho(m))*(B*(mu1reuse_mean+1)+sigma2rho(m))/((1-epsilon2(m))*Mvalues(n)) + mu2reuse_mean*B + B*variance_mean*(1/((1-epsilon2(m))*Mvalues(n))));
                    SE_MR_mean(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_MR_mean);
                    
                    SINR_MR_worst = B*(1-epsilon2(m))/(epsilon2(m)*B + (mu1all_worst*K + sigma2rho(m))*(B*(mu1reuse_worst+1)+sigma2rho(m))/((1-epsilon2(m))*Mvalues(n)) + mu2reuse_worst*B + B*variance_worst*(1/((1-epsilon2(m))*Mvalues(n))));
                    SE_MR_worst(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_MR_worst);
                    
                    SINR_MR_best = B*(1-epsilon2(m))/(epsilon2(m)*B + (mu1all_best*K + sigma2rho(m))*(B*(mu1reuse_best+1)+sigma2rho(m))/((1-epsilon2(m))*Mvalues(n)) + mu2reuse_best*B + B*variance_best*(1/((1-epsilon2(m))*Mvalues(n))));
                    SE_MR_best(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_MR_best);
                    
                    
                    %Zero-forcing (ZF) combining/precoding
                    %
                    %Achievable spectral efficiency using the formula in
                    %Theorem 1, for mean, worst, and best case interference
                    if Mvalues(n)-K>0
                        
                        %Compute one of the terms in Theorem 1
                        term2_ZF_mean = (directions*reuseMu1Mean2(reuseIndices(j))+1^2)/(B*(mu1reuse_mean+1)+sigma2rho(m));
                        term2_ZF_worst = (directions*reuseMu1Worst2(reuseIndices(j))+1^2)/(B*(mu1reuse_worst+1)+sigma2rho(m));
                        term2_ZF_best = (directions*reuseMu1Best2(reuseIndices(j))+1^2)/(B*(mu1reuse_best+1)+sigma2rho(m));
                        
                        SINR_ZF_mean = B*(1-epsilon2(m))/(epsilon2(m)*B + mu2reuse_mean*B + B*variance_mean/(Mvalues(n)-K)/(1-epsilon2(m)) +  (K*(mu1all_mean - (1-epsilon2(m))*B*term2_ZF_mean) + sigma2rho(m) )*(B*(mu1reuse_mean+1)+sigma2rho(m))/(Mvalues(n)-K)/(1-epsilon2(m)) );
                        SE_ZF_mean(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_ZF_mean);
                        
                        SINR_ZF_worst = B*(1-epsilon2(m))/(epsilon2(m)*B + mu2reuse_worst*B + B*variance_worst/(Mvalues(n)-K)/(1-epsilon2(m)) +  (K*(mu1all_worst - (1-epsilon2(m))*B*term2_ZF_worst) + sigma2rho(m) )*(B*(mu1reuse_worst+1)+sigma2rho(m))/(Mvalues(n)-K)/(1-epsilon2(m)) );
                        SE_ZF_worst(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_ZF_worst);
                        
                        SINR_ZF_best = B*(1-epsilon2(m))/(epsilon2(m)*B + mu2reuse_best*B + B*variance_best/(Mvalues(n)-K)/(1-epsilon2(m)) +  (K*(mu1all_best - (1-epsilon2(m))*B*term2_ZF_best) + sigma2rho(m) )*(B*(mu1reuse_best+1)+sigma2rho(m))/(Mvalues(n)-K)/(1-epsilon2(m)) );
                        SE_ZF_best(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_ZF_best);
                        
                    end
                    
                    
                    %Full-pilot zero-forcing (P-ZF) combining/precoding
                    %
                    %Achievable spectral efficiency using the formula in
                    %Theorem 2, for mean, worst, and best case interference
                    if Mvalues(n)-B>0
                        
                        %Compute one of the terms in Theorem 2
                        term2_PZF_mean = (directions*reuseMu1Mean2(reuseIndices(j))+1^2)/(B*(mu1reuse_mean+1)+sigma2rho(m)) + directions*reuseMu1Mean2Next(reuseIndices(j)) /(B*(neighborsPerOtherPilot*reuseMu1MeanNext(reuseIndices(j)))+sigma2rho(m));
                        term2_PZF_worst = (directions*reuseMu1Worst2(reuseIndices(j))+1^2)/(B*(mu1reuse_worst+1)+sigma2rho(m)) + directions*reuseMu1Worst2Next(reuseIndices(j)) /(B*(neighborsPerOtherPilot*reuseMu1WorstNext(reuseIndices(j)))+sigma2rho(m));
                        term2_PZF_best = (directions*reuseMu1Best2(reuseIndices(j))+1^2)/(B*(mu1reuse_best+1)+sigma2rho(m)) + directions*reuseMu1Best2Next(reuseIndices(j)) /(B*(neighborsPerOtherPilot*reuseMu1BestNext(reuseIndices(j)))+sigma2rho(m));
                        
                        SINR_PZF_mean = B*(1-epsilon2(m))/(epsilon2(m)*B + mu2reuse_mean*B + B*variance_mean/(Mvalues(n)-B)/(1-epsilon2(m)) +  (K*(mu1all_mean - (1-epsilon2(m))*B*term2_PZF_mean) + sigma2rho(m) )*(B*(mu1reuse_mean+1)+sigma2rho(m))/(Mvalues(n)-B)/(1-epsilon2(m)) );
                        SE_PZF_mean(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_PZF_mean);
                        
                        SINR_PZF_worst = B*(1-epsilon2(m))/(epsilon2(m)*B + mu2reuse_worst*B + B*variance_worst/(Mvalues(n)-B)/(1-epsilon2(m)) +  (K*(mu1all_worst - (1-epsilon2(m))*B*term2_PZF_worst) + sigma2rho(m) )*(B*(mu1reuse_worst+1)+sigma2rho(m))/(Mvalues(n)-B)/(1-epsilon2(m)) );
                        SE_PZF_worst(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_PZF_worst);
                        
                        SINR_PZF_best = B*(1-epsilon2(m))/(epsilon2(m)*B + mu2reuse_best*B + B*variance_best/(Mvalues(n)-B)/(1-epsilon2(m)) +  (K*(mu1all_best - (1-epsilon2(m))*B*term2_PZF_best) + sigma2rho(m) )*(B*(mu1reuse_best+1)+sigma2rho(m))/(Mvalues(n)-B)/(1-epsilon2(m)) );
                        SE_PZF_best(n,K,j,m) = K*(1-B/S(m))*log2(1+SINR_PZF_best);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end



%%Compute optimal number of UEs, K, for different system parameters

%Placeholders for storing simulation results for optimal number of UEs
optimalK_MR_mean = zeros(length(Mvalues),3,length(S));
optimalK_ZF_mean = zeros(length(Mvalues),3,length(S));
optimalK_PZF_mean = zeros(length(Mvalues),3,length(S));

optimalK_MR_worst = zeros(length(Mvalues),3,length(S));
optimalK_ZF_worst = zeros(length(Mvalues),3,length(S));
optimalK_PZF_worst = zeros(length(Mvalues),3,length(S));

optimalK_MR_best = zeros(length(Mvalues),3,length(S));
optimalK_ZF_best = zeros(length(Mvalues),3,length(S));
optimalK_PZF_best = zeros(length(Mvalues),3,length(S));


%Go through different number of antennas
for n = 1:length(Mvalues)
    
    %Go through different reuse factors
    for j = 1:length(reuseIndices)
        
        currentReuseFactor = reuseFactor(reuseIndices(j));
        
        %Go through different other cases (varying SNR, coherence block, etc.)
        for m = 1:length(S)
            
            %Optimize for average interference case
            [maxValue,maxIndex] = max(SE_MR_mean(n,:,j,m));
            if maxValue > optimalK_MR_mean(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_MR_mean(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            [maxValue,maxIndex] = max(SE_ZF_mean(n,:,j,m));
            if maxValue > optimalK_ZF_mean(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_ZF_mean(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            [maxValue,maxIndex] = max(SE_PZF_mean(n,:,j,m));
            if maxValue > optimalK_PZF_mean(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_PZF_mean(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            
            %Optimize for worst interference case
            [maxValue,maxIndex] = max(SE_MR_worst(n,:,j,m));
            if maxValue > optimalK_MR_worst(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_MR_worst(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            [maxValue,maxIndex] = max(SE_ZF_worst(n,:,j,m));
            if maxValue > optimalK_ZF_worst(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_ZF_worst(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            [maxValue,maxIndex] = max(SE_PZF_worst(n,:,j,m));
            if maxValue > optimalK_PZF_worst(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_PZF_worst(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            
            %Optimize for best interference case
            [maxValue,maxIndex] = max(SE_MR_best(n,:,j,m));
            if maxValue > optimalK_MR_best(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_MR_best(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            [maxValue,maxIndex] = max(SE_ZF_best(n,:,j,m));
            if maxValue > optimalK_ZF_best(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_ZF_best(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
            [maxValue,maxIndex] = max(SE_PZF_best(n,:,j,m));
            if maxValue > optimalK_PZF_best(n,2,m)
                %Store optimal number of UEs along with the optimized SE
                %and the corresponding reuse factor
                optimalK_PZF_best(n,:,m) = [maxIndex maxValue currentReuseFactor];
            end
            
        end
        
    end
    
end



%%Plot simulation results

if simulationCase == 1 %Simulations from Section IV.A
    
    
    %Plot Figure 4(a) and (b)
    figure(4);
    
    subplot(2,1,1);
    hold on; box on;
    
    plot(Mvalues,max(max(SE_asymptotic_mean(:,:,1)))*ones(size(Mvalues)),'k:','LineWidth',1);
    plot(Mvalues,optimalK_PZF_mean(:,2,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_mean(:,2,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_mean(:,2,1),'b-.','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('Asymptotic limit','P-ZF','ZF','MR','Location','NorthWest');
    set(gca,'Xscale','log');
    axis([10 1e5 0 400]);
    
    
    subplot(2,1,2);
    hold on; box on;
    
    plot(Mvalues,optimalK_PZF_mean(:,1,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_mean(:,1,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_mean(:,1,1),'b-.','LineWidth',1);
    plot(Mvalues,Mvalues,'k:','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Optimal Number of Scheduled UEs (K)');
    legend('P-ZF','ZF','MR','Location','NorthWest');
    set(gca,'Xscale','log');
    axis([10 1e5 0 200]);
    
    
    
    %Plot Figure 5(a) and (b)
    figure(5);
    
    subplot(2,1,1);
    hold on; box on;
    
    plot(Mvalues,max(max(SE_asymptotic_best(:,:,1)))*ones(size(Mvalues)),'k:','LineWidth',1);
    plot(Mvalues,optimalK_PZF_best(:,2,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_best(:,2,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_best(:,2,1),'b-.','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('Asymptotic limit','P-ZF','ZF','MR','Location','NorthWest');
    set(gca,'Xscale','log');
    axis([10 1e5 0 2600]);
    
    
    subplot(2,1,2);
    hold on; box on;
    
    plot(Mvalues,optimalK_PZF_best(:,1,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_best(:,1,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_best(:,1,1),'b-.','LineWidth',1);
    plot(Mvalues,Mvalues,'k:','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Optimal Number of Scheduled UEs (K)');
    legend('P-ZF','ZF','MR','Location','NorthWest');
    set(gca,'Xscale','log');
    axis([10 1e5 0 200]);
    
    
    
    %Plot Figure 6(a) and (b)
    figure(6);
    
    subplot(2,1,1);
    hold on; box on;
    
    plot(Mvalues,max(max(SE_asymptotic_worst(:,:,1)))*ones(size(Mvalues)),'k:','LineWidth',1);
    plot(Mvalues,optimalK_PZF_worst(:,2,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_worst(:,2,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_worst(:,2,1),'b-.','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('Asymptotic limit','P-ZF','ZF','MR','Location','NorthWest');
    set(gca,'Xscale','log');
    
    
    subplot(2,1,2);
    hold on; box on;
    
    plot(Mvalues,optimalK_PZF_worst(:,1,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_worst(:,1,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_worst(:,1,1),'b-.','LineWidth',1);
    plot(Mvalues,Mvalues,'k:','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Optimal Number of Scheduled UEs (K)');
    legend('P-ZF','ZF','MR','Location','SouthEast');
    set(gca,'Xscale','log');
    axis([10 1e5 0 60]);
    
    
    
    %Plot Figure 14
    figure(14);
    hold on; box on;
    
    for m = 1:length(epsilon2)
        plot(Mvalues,max(max(SE_asymptotic_mean(:,:,m)))*ones(size(Mvalues)),'k:','LineWidth',1);
        plot(Mvalues,optimalK_PZF_mean(:,2,m),'r-','LineWidth',1);
        plot(Mvalues,optimalK_ZF_mean(:,2,m),'k--','LineWidth',1);
        plot(Mvalues,optimalK_MR_mean(:,2,m),'b-.','LineWidth',1);
    end
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('Asymptotic limit','P-ZF','ZF','MR','Location','SouthEast');
    set(gca,'Xscale','log');
    axis([10 1e5 0 400]);
    
    
elseif simulationCase == 2
    
    
    %Plot Figure 8
    figure(8); hold on; box on;
    
    for j = [1 3];
        
        optAreaRates = max(SE_PZF_mean(:,:,j),[],2);
        plot(Mvalues,optAreaRates,'r','LineWidth',1);
        text(200,optAreaRates(200)+5,num2str(reuseFactor(reuseIndices(j))));
        
        optAreaRates = max(SE_ZF_mean(:,:,j),[],2);
        plot(Mvalues,optAreaRates,'k--','LineWidth',1);
        text(200,optAreaRates(200)+5,num2str(reuseFactor(reuseIndices(j))));
        
        
        optAreaRates = max(SE_MR_mean(:,:,j),[],2);
        plot(Mvalues,optAreaRates,'b-.','LineWidth',1);
        text(400,optAreaRates(400)+3,num2str(reuseFactor(reuseIndices(j))));
        
    end
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('P-ZF','ZF','MR','Location','NorthWest');
    
    
    %Plot Figure 9
    figure(9); hold on; box on;
    
    plot(Mvalues,optimalK_PZF_mean(:,2)./optimalK_PZF_mean(:,1),'r-','LineWidth',1);
    plot(Mvalues,optimalK_ZF_mean(:,2)./optimalK_ZF_mean(:,1),'k--','LineWidth',1);
    plot(Mvalues,optimalK_MR_mean(:,2)./optimalK_MR_mean(:,1),'b-.','LineWidth',1);
    
    xlabel('Number of BS Antennas (M)');
    ylabel('Per-UE Spectral Efficiency [bit/s/Hz/user]');
    legend('P-ZF','ZF','MR','Location','NorthEast');
    axis([0 1000 0 3]);
    
    
    %Plot Figure 10
    figure(10); hold on; box on;
    
    plot(Mvalues,Mvalues./optimalK_PZF_mean(:,1)','r-','LineWidth',1);
    plot(Mvalues,Mvalues./optimalK_ZF_mean(:,1)','k--','LineWidth',1);
    plot(Mvalues,Mvalues./optimalK_MR_mean(:,1)','b-.','LineWidth',1);
    
    plot([0 1000],[10 10],'k:');
    
    xlabel('Number of BS Antennas (M)');
    ylabel('BS Antennas per UE (M/K)');
    legend('P-ZF','ZF','MR','Location','NorthEast');
    axis([0 1000 0 14]);
    
    
    %Plot Figure 11
    figure(11); hold on; box on;
    
    mValues = [100 500];
    
    for m = 1:length(mValues)
        
        plot(Kvalues,max(SE_PZF_mean(mValues(m),:,:),[],3),'r-','LineWidth',1);
        plot(Kvalues,max(SE_ZF_mean(mValues(m),:,:),[],3),'k--','LineWidth',1);
        plot(Kvalues,max(SE_MR_mean(mValues(m),:,:),[],3),'b-.','LineWidth',1);
        
        [SEmax,kmax] = max(max(SE_PZF_mean(mValues(m),:,:),[],3));
        plot(kmax,SEmax,'r*','LineWidth',1);
        
        [SEmax,kmax] = max(max(SE_ZF_mean(mValues(m),:,:),[],3));
        plot(kmax,SEmax,'k*','LineWidth',1);
        
        [SEmax,kmax] = max(max(SE_MR_mean(mValues(m),:,:),[],3));
        plot(kmax,SEmax,'b*','LineWidth',1);
        
    end
    
    xlabel('Number of UEs (K)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('P-ZF','ZF','MR','Location','Best');
    axis([0 S(1)/2 0 140]);
    
    
elseif simulationCase == 3
    
    
    %Plot Figure 12
    figure(12); hold on; box on;
    
    for m = 1:length(Mvalues)
        plot(SNRvaluesdB,reshape(max(max(SE_PZF_mean(m,:,:,:),[],3),[],2),[length(SNRvaluesdB) 1]),'r-','LineWidth',1);
        plot(SNRvaluesdB,reshape(max(max(SE_ZF_mean(m,:,:,:),[],3),[],2),[length(SNRvaluesdB) 1]),'k--','LineWidth',1);
        plot(SNRvaluesdB,reshape(max(max(SE_MR_mean(m,:,:,:),[],3),[],2),[length(SNRvaluesdB) 1]),'b-.','LineWidth',1);
    end
    
    xlabel('Signal-to-Noise Ratio (SNR), \rho/\sigma^2 ');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('P-ZF','ZF','MR','Location','SouthEast');
    axis([min(SNRvaluesdB) max(SNRvaluesdB) 0 150]);
    
    
elseif simulationCase == 4
    
    
    %Plot Figure 13
    figure(13); hold on; box on;
    
    for m = 1:length(Mvalues)
        plot(S,reshape(max(max(SE_PZF_mean(m,:,:,:),[],3),[],2),[length(S) 1]),'r-','LineWidth',1);
        plot(S,reshape(max(max(SE_ZF_mean(m,:,:,:),[],3),[],2),[length(S) 1]),'k--','LineWidth',1);
        plot(S,reshape(max(max(SE_MR_mean(m,:,:,:),[],3),[],2),[length(S) 1]),'b-.','LineWidth',1);
    end
    
    xlabel('Coherence Block Length (S)');
    ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');
    legend('P-ZF','ZF','MR','Location','NorthWest');
    axis([0 max(S) 0 300]);
    
    
end
