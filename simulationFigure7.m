%This Matlab script can be used to generate Figure 7 in the article:
%
%Emil Bjornson, Erik G. Larsson, Merouane Debbah, "Massive MIMO for Maximal
%Spectral Efficiency: How Many Users and Pilots Should Be Allocated?,"
%vol. 15, no. 2, pp. 1293-1308, February 2016.
%
%Download article: http://arxiv.org/pdf/1412.7102
%
%This is version 1.2 (Last edited: 2016-08-22)
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


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

%Maximal number of Monte-Carlo realizations with random user locations
monteCarloRealizations = 5000;

%Select number of UEs
K = 10;

%Select range of BS antennas
Mvalues = 10:30:1000;

%Pathloss exponent
kappa = 3.7;

%SNR value (5 dB)
SNR = 10^(5/10);

%Coherence block length
S = 400;

%Number of tiers of hexagonals that are simulated, around the desired cell
tiers = 7;

%Percentage of the radius inside the cell where no UEs are allowed
forbiddenRegion = 0.14;

%Define intersite distance in a normalized scale
intersiteDistance = 2;
intersiteDistanceHalf = intersiteDistance/2;

dmax = intersiteDistanceHalf; %Normalized cell radius
dmin = dmax * forbiddenRegion; %Normalized shortest distance from a BS



%%Begin Monte-Carlo simulations

%Placeholders for storing spectral efficiencies from Monte-Carlo
%simulations
averageSEsReuse1_MR = zeros(monteCarloRealizations,length(Mvalues));
averageSEsReuse3_MR = zeros(monteCarloRealizations,length(Mvalues));
averageSEsReuse4_MR = zeros(monteCarloRealizations,length(Mvalues));

averageSEsReuse1_ZF = zeros(monteCarloRealizations,length(Mvalues));
averageSEsReuse3_ZF = zeros(monteCarloRealizations,length(Mvalues));
averageSEsReuse4_ZF = zeros(monteCarloRealizations,length(Mvalues));

averageSEsReuse1_PZF = zeros(monteCarloRealizations,length(Mvalues));
averageSEsReuse3_PZF = zeros(monteCarloRealizations,length(Mvalues));
averageSEsReuse4_PZF = zeros(monteCarloRealizations,length(Mvalues));




%Go through all Monte-Carlo realizations
for n = 1:monteCarloRealizations
    
    %Display simulation progress
    disp(['Realization ' num2str(n) ' out of ' num2str(monteCarloRealizations)]);
    
    
    %Vector with BS positions
    BSpositions = zeros(4,1);
    
    %Vectors that store which set of pilots that each BS uses (when there
    %is non-universal pilot reuse)
    reusePattern3 = zeros(4,1); %Pilot reuse 3
    reusePattern4 = zeros(4,1); %Pilot reuse 4
    
    
    %BS index of the next BS to be deployed
    itr = 1;
    
    %Deploy hexagonal cells in "tiers" number of tiers
    for alpha1 = 0:tiers
        for alpha2 = 0:tiers-alpha1
            
            if (alpha1 == 0) || (alpha2>=1)
                
                %Compute a BS location according to Eq. (30)
                BSloc = sqrt(3)*alpha2*intersiteDistanceHalf*1i + sqrt(3)*alpha1*intersiteDistanceHalf*exp(1i*pi*(30/180));
                
                
                %Special: The first BS is placed in the origin (this is
                %where the performance is computed)
                if (alpha1 == 0) && (alpha2 == 0)
                    BSpositions(itr) = BSloc;
                    reusePattern3(itr) = 1;
                    reusePattern4(itr) = 1;
                    itr = itr+1;
                    
                else
                    
                    %Compute the current BS location
                    basis = [3*intersiteDistanceHalf/2 0; sqrt(3)*intersiteDistanceHalf/2 sqrt(3)*intersiteDistanceHalf];
                    rotation = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)];
                    
                    %Deploy BSs in all six directions from the center cell
                    for m = 1:6
                        
                        %Compute the reuse pattern identity for reuse 3
                        if mod(2*alpha1+alpha2,3)==0
                            reusePattern3(itr) = 1;
                        elseif mod(mod(2*alpha1+alpha2,3)+m,2)==1
                            reusePattern3(itr) = 2;
                        elseif mod(mod(2*alpha1+alpha2,3)+m,2)==0
                            reusePattern3(itr) = 3;
                        end
                        
                        %Compute the reuse pattern identity for reuse 4
                        if mod(alpha1,2)==0 && mod(alpha2,2)==0
                            reusePattern4(itr) = 1;
                        elseif m == 1
                            reusePattern4(itr) = 1+mod(2*alpha1+mod(alpha2,2),4);
                        else
                            dims = round(basis\(rotation^(m-1)*basis*[alpha1; alpha2]));
                            reusePattern4(itr) = 1+mod(2*dims(1)+mod(dims(2),2),4);
                        end
                        
                        %Deploy a BS
                        BSpositions(itr) = BSloc;
                        itr = itr+1;
                        
                        %Rotate the BS location to consider the next
                        %direction of the six directions from the center
                        %cell
                        BSloc = BSloc .* exp(1i*2*pi*60/360);
                        
                    end
                end
                
            end
            
        end
    end
    
    %Compute the final number of BSs
    nbrBSs = itr-1;
    
    %Prepare to put out UEs in the cells
    UEpositions = zeros(K,nbrBSs);
    
    %Initiate matrices where first and second order interference are computed
    interference1reuse1 = zeros(K,1);
    interference2reuse1 = zeros(K,1);
    interference1reuse3 = zeros(K,3);
    interference2reuse3 = zeros(K,3);
    interference1reuse4 = zeros(K,4);
    interference2reuse4 = zeros(K,4);
    
    
    %Go through all the cells
    for j = 1:nbrBSs
        
        %Generate UE locations randomly with uniform distribution inside the cells
        nbrToGenerate = K; %Number of UE locations left to generate
        notFinished = true(K,1); %Indices of the UE locations that are left to generate
        
        
        %Iterate the generation of UE locations until all of them are inside a
        %hexagonal cell
        while nbrToGenerate>0
            
            %Generate new UE locations uniformly at random in a circle of radius dmax
            UEpositions(notFinished,j) = sqrt( rand(nbrToGenerate,1)*(dmax^2-dmin^2)+ dmin^2 ) .* exp(1i*2*pi*rand(nbrToGenerate,1));
            
            %Check which UEs that are inside a hexagonal and declare as finished
            finished = checkHexagonal(UEpositions(:,j)',dmax);
            
            %Update which UEs that are left to generate
            notFinished = (finished==false);
            
            %Update how many UEs that are left to generate
            nbrToGenerate = sum(notFinished);
            
        end
        
        %Finalize UE locations by translating them around the serving BS
        UEpositions(:,j) = UEpositions(:,j) + BSpositions(j);
        
        %Compute the distance from the users in cell j to BS j
        distancesSquaredBSj = abs(UEpositions(:,j) - BSpositions(j));
        
        
        %Focus on interference caused to the first BS (in the center)
        l = 1;
        
        %Compute the distance from the users in cell j to BS l
        distancesSquaredBSl = abs(UEpositions(:,j) - BSpositions(l));
        
        
        %Compute inteference terms of the types that show up in Eq. (7)
        interference1reuse1(:,l) = interference1reuse1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse1(:,l) = interference2reuse1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
        interference1reuse3(:,reusePattern3(j)) = interference1reuse3(:,reusePattern3(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse3(:,reusePattern3(j)) = interference2reuse3(:,reusePattern3(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
        
        interference1reuse4(:,reusePattern4(j)) = interference1reuse4(:,reusePattern4(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse4(:,reusePattern4(j)) = interference2reuse4(:,reusePattern4(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
    end
    
    
    %Go through all different number of BS antennas
    for mind = 1:length(Mvalues)
        
        m = Mvalues(mind);
        
        %Compute different pilot lengths (depending on pilot reuse factor)
        B1 = K;
        B3 = 3*K;
        B4 = 4*K;
        
        %Prepare to store the SINRs
        SINRreuse1_MR = zeros(K,1);
        SINRreuse3_MR = zeros(K,1);
        SINRreuse4_MR = zeros(K,1);
        
        SINRreuse1_ZF = zeros(K,1);
        SINRreuse3_ZF = zeros(K,1);
        SINRreuse4_ZF = zeros(K,1);
        
        SINRreuse1_PZF = zeros(K,1);
        SINRreuse3_PZF = zeros(K,1);
        SINRreuse4_PZF = zeros(K,1);
        
        
        %Compute uplink performance in center cell
        j = 1;
        
        %Compute the SINRs in accordance to Lemma 2 for MR combining
        SINRreuse1_MR(:,j) = m *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) + 1/SNR)*(interference1reuse1(1:K,j)+1/(B1*SNR)) + m*(interference2reuse1(1:K,j)-1));
        SINRreuse3_MR(:,j) = m *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) + 1/SNR)*(interference1reuse3(1:K,j)+1/(B3*SNR)) + m*(interference2reuse3(1:K,j)-1));
        SINRreuse4_MR(:,j) = m *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) + 1/SNR)*(interference1reuse4(1:K,j)+1/(B4*SNR)) + m*(interference2reuse4(1:K,j)-1));
        
        %Compute the SINRs in accordance to Lemma 2 for ZF combining
        if (m-K)>0 %Check if ZF exists
            SINRreuse1_ZF(:,j) = (m-K) *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) - sum(interference2reuse1(1:K,j)./(interference1reuse1(1:K,j)+1/(B1*SNR))) + 1/SNR)*( interference1reuse1(1:K,j)  +1/(B1*SNR)) + (m-K)*(interference2reuse1(1:K,j)-1));
            SINRreuse3_ZF(:,j) = (m-K) *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) - sum(interference2reuse3(1:K,j)./(interference1reuse3(1:K,j)+1/(B3*SNR))) + 1/SNR)*( interference1reuse3(1:K,j) +1/(B3*SNR)) + (m-K)*(interference2reuse3(1:K,j)-1));
            SINRreuse4_ZF(:,j) = (m-K) *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) - sum(interference2reuse4(1:K,j)./(interference1reuse4(1:K,j)+1/(B4*SNR))) + 1/SNR)*( interference1reuse4(1:K,j) +1/(B4*SNR)) + (m-K)*(interference2reuse4(1:K,j)-1));
        end
        
        %Compute the SINRs in accordance to Lemma 2 for PZF combining
        if (m-B1)>0 %Check if PZF exists
            SINRreuse1_PZF(:,j) = (m-B1) *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) - sum(interference2reuse1(1:K,j)./(interference1reuse1(1:K,j)+1/(B1*SNR))) + 1/SNR)*( interference1reuse1(1:K,j)  +1/(B1*SNR)) + (m-B1)*(interference2reuse1(1:K,j)-1));
        end
        
        if (m-B3)>0 %Check if PZF exists
            SINRreuse3_PZF(:,j) = (m-B3) *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) - sum(sum(interference2reuse3(1:K,:)./(interference1reuse3(1:K,:)+1/(B3*SNR)))) + 1/SNR)*( interference1reuse3(1:K,j) +1/(B3*SNR)) + (m-B3)*(interference2reuse3(1:K,j)-1));
        end
        
        if (m-B4)>0 %Check if PZF exists
            SINRreuse4_PZF(:,j) = (m-B4) *ones(K,1) ./ ( (sum(interference1reuse1(1:K,j)) - sum(sum(interference2reuse4(1:K,:)./(interference1reuse4(1:K,:)+1/(B4*SNR)))) + 1/SNR)*( interference1reuse4(1:K,j)  +1/(B4*SNR)) + (m-B4)*(interference2reuse4(1:K,j)-1));
        end
        
        
        %Compute the average SEs over the network by Monte-Carlo simulation
        averageSEsReuse1_MR(n,mind) = mean(log2(1+SINRreuse1_MR(:)));
        averageSEsReuse3_MR(n,mind) = mean(log2(1+SINRreuse3_MR(:)));
        averageSEsReuse4_MR(n,mind) = mean(log2(1+SINRreuse4_MR(:)));
        
        averageSEsReuse1_ZF(n,mind) = mean(log2(1+SINRreuse1_ZF(:)));
        averageSEsReuse3_ZF(n,mind) = mean(log2(1+SINRreuse3_ZF(:)));
        averageSEsReuse4_ZF(n,mind) = mean(log2(1+SINRreuse4_ZF(:)));
        
        averageSEsReuse1_PZF(n,mind) = mean(log2(1+SINRreuse1_PZF(:)));
        averageSEsReuse3_PZF(n,mind) = mean(log2(1+SINRreuse3_PZF(:)));
        averageSEsReuse4_PZF(n,mind) = mean(log2(1+SINRreuse4_PZF(:)));
        
    end
    
end



%%Begin analytical computation using Theorem 1


%Parameters for the Monte Carlo computations of mu-parameters
monteCarloUEs = 1000000; %Number of random UE locations per cell

%Compute various combinations of the mu-parameters Eq. (18), using
%Monte Carlo simulations
[muValues1Mean,muValues2Mean,reuseMu1Mean,reuseMu1Mean2,reuseMu1MeanNext,reuseMu1Mean2Next,reuseMu2Mean,reuseMuMeanVariance,muValues1Worst,muValues2Worst,reuseMu1Worst,reuseMu1Worst2,reuseMu1WorstNext,reuseMu1Worst2Next,reuseMu2Worst,reuseMuWorstVariance,muValues1Best,muValues2Best,reuseMu1Best,reuseMu1Best2,reuseMu1BestNext,reuseMu1Best2Next,reuseMu2Best,reuseMuBestVariance,reuseFactor] = computeEnvironment(kappa,forbiddenRegion,monteCarloUEs);

%Number of directions to look for interfering cells (for hexagonal cells)
directions = 6;

%Compute the sum of all mu values in Eq. (18)
mu1all_mean = 1+directions*(sum(muValues1Mean(:))-1);
mu2all_mean = 1+directions*(sum(muValues2Mean(:))-1);
mu1all_mean2 = 1+directions*(sum(muValues1Mean(:).^2)-1);

%Extract only reuse factors smaller or equal to 4
reuseIndices = find(reuseFactor>0 & reuseFactor<=4);
for j = 1:length(reuseIndices);
    if sum(reuseFactor(reuseIndices(j))==reuseFactor(reuseIndices(1:j-1)))>0
        reuseIndices(j)=1;
    end
end
reuseIndices = reuseIndices(reuseIndices>1);


%Placeholders for storing spectral efficiencies
SE_MR_mean = zeros(length(Mvalues),length(reuseIndices));
SE_ZF_mean = zeros(length(Mvalues),length(reuseIndices));
SE_PZF_mean = zeros(length(Mvalues),length(reuseIndices));


%Go through the different reuse factors
for j = 1:length(reuseIndices);
    
    %Extract the reuse factor
    currentReuseFactor = reuseFactor(reuseIndices(j));
    
    %Extract sum of mu-values for current reuse factor for mean interference
    mu1reuse_mean = directions*reuseMu1Mean(reuseIndices(j));
    mu2reuse_mean = directions*reuseMu2Mean(reuseIndices(j));
    variance_mean = directions*reuseMuMeanVariance(reuseIndices(j));
    
    %Number of neighbors that use each of the other sets of pilots
    neighborsPerOtherPilot = directions/(currentReuseFactor-1);
    
    %Go through different number of antennas
    for n = 1:length(Mvalues)
        
        
        %Compute length of pilot signal
        B1 = currentReuseFactor*K;
        
        if B1 < S
            
            %Compute achievable spectral efficiency using the formula in
            %Theorem 1 for MR combining
            SINR_MR_mean = Mvalues(n)/( (mu1all_mean*(K) + 1/SNR)*((mu1reuse_mean+1)+1/(B1*SNR)) + Mvalues(n)*mu2reuse_mean + variance_mean/B1);
            SE_MR_mean(n,j) = K*(1-B1/S)*log2(1+SINR_MR_mean);
            
            
            %Compute achievable spectral efficiency using the formula in
            %Theorem 1 for ZF combining
            if Mvalues(n)-K>0
                term_ZF_mean = (directions*reuseMu1Mean2(reuseIndices(j))+1^2)/(B1*(mu1reuse_mean+1)+1/SNR);
                SINR_ZF_mean = B1/( mu2reuse_mean*B1 + B1*variance_mean/(Mvalues(n)-K) +  (K*(mu1all_mean - B1*term_ZF_mean) + 1/SNR )*(B1*(mu1reuse_mean+1)+1/SNR)/(Mvalues(n)-K) );
                SE_ZF_mean(n,j) = K*(1-B1/S)*log2(1+SINR_ZF_mean);
            end
            
            %Compute achievable spectral efficiency using the formula in
            %Theorem 1 for PZF combining
            if Mvalues(n)-B1>0
                term_PZF_mean = (directions*reuseMu1Mean2(reuseIndices(j))+1^2)/(B1*(mu1reuse_mean+1)+1/SNR) + directions*reuseMu1Mean2Next(reuseIndices(j)) /(B1*(neighborsPerOtherPilot*reuseMu1MeanNext(reuseIndices(j)))+1/SNR);
                SINR_PZF_mean = B1/(mu2reuse_mean*B1 + B1*variance_mean/(Mvalues(n)-B1) +  (K*(mu1all_mean - B1*term_PZF_mean) + 1/SNR )*(B1*(mu1reuse_mean+1)+1/SNR)/(Mvalues(n)-B1) );
                SE_PZF_mean(n,j) = K*(1-B1/S)*log2(1+SINR_PZF_mean);
            end
            
        end
    end
end



%Plot Figure 7
figure(7); hold on; box on;

%Plot analytical results from Theorem 1
plot(Mvalues,max(SE_PZF_mean,[],2),'r-','LineWidth',1);
plot(Mvalues,max(SE_ZF_mean,[],2),'k--','LineWidth',1);
plot(Mvalues,max(SE_MR_mean,[],2),'b-.','LineWidth',1);

%Plot Monte-Carlo results for PZF
SEreuse1 = K*(1-K/S)*mean(averageSEsReuse1_PZF,1);
SEreuse3 = K*(1-3*K/S)*mean(averageSEsReuse3_PZF,1);
SEreuse4 = K*(1-4*K/S)*mean(averageSEsReuse4_PZF,1);
SEoptimzed = max([SEreuse1(:) SEreuse3(:) SEreuse4(:)],[],2);

plot(Mvalues,SEoptimzed,'ro');

%Plot Monte-Carlo results for ZF
SEreuse1 = K*(1-K/S)*mean(averageSEsReuse1_ZF,1);
SEreuse3 = K*(1-3*K/S)*mean(averageSEsReuse3_ZF,1);
SEreuse4 = K*(1-4*K/S)*mean(averageSEsReuse4_ZF,1);
SEoptimzed = max([SEreuse1(:) SEreuse3(:) SEreuse4(:)],[],2);

plot(Mvalues,SEoptimzed,'kd');

%Plot Monte-Carlo results for MR
SEreuse1 = K*(1-K/S)*mean(averageSEsReuse1_MR,1);
SEreuse3 = K*(1-3*K/S)*mean(averageSEsReuse3_MR,1);
SEreuse4 = K*(1-4*K/S)*mean(averageSEsReuse4_MR,1);
SEoptimzed = max([SEreuse1(:) SEreuse3(:) SEreuse4(:)],[],2);

plot(Mvalues,SEoptimzed,'b*');

axis([0 1000 0 80]);
xlabel('Number of BS Antennas (M)');
ylabel('Spectral Efficiency (SE) [bit/s/Hz/cell]');

legend('P-ZF','ZF','MR','Location','NorthWest');
