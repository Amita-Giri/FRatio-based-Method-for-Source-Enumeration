function [OptimumThresh,BestAccuracy,Accuracy] = AccuracyandThreshold(MaxnumSources,THRESH,MC_repetitions,Fratio_RAP,Fratio_AP)
for True_numS1 = 0:MaxnumSources
    %% Estimated number of sources with varying threshold value
    temp = 0;
    for Threshold = THRESH
        temp = temp+1;
        for True_nsources = True_numS1
            for idxMC= 1:MC_repetitions
                %% RAP
                for k= 0:length(Fratio_RAP{True_nsources+1,idxMC})-1
                    if Fratio_RAP{True_nsources+1,idxMC}(k+1,1) < Threshold
                        Est_no_sources_RAP{temp,1}(True_nsources+1,idxMC)=k;
                        break
                    else
                        Est_no_sources_RAP{temp,1}(True_nsources+1,idxMC)=length(Fratio_RAP{True_nsources+1,idxMC})-1;
                    end
                end
                %% AP
                for k= 0:length(Fratio_AP{True_nsources+1,idxMC})-1
                    if Fratio_AP{True_nsources+1,idxMC}(k+1,1) < Threshold
                        Est_no_sources_AP{temp,1}(True_nsources+1,idxMC)=k;
                        break
                    else
                        Est_no_sources_AP{temp,1}(True_nsources+1,idxMC)=length(Fratio_AP{True_nsources+1,idxMC})-1;
                    end
                end
            end
        end
    end

    %% Accuracy with varying threshold value
    temp = 0;
    for Threshold = THRESH
        temp = temp+1;
        for True_nsources = True_numS1
            Count_TRAP = 0;  Count_RAP = 0;  Count_AP = 0;
            for idxMC= 1:MC_repetitions
                %% RAP
                if Est_no_sources_RAP{temp,1}(True_nsources+1,idxMC)==True_nsources
                    Count_RAP = Count_RAP+1;
                end
                %% AP
                if Est_no_sources_AP{temp,1}(True_nsources+1,idxMC)==True_nsources
                    Count_AP = Count_AP+1;
                end
            end
            Accuracy{True_nsources+1,temp} = [Count_TRAP Count_RAP Count_AP];
        end
    end

    %% Optimum threshold and accuracy value
    temp = 0;
    for Threshold = THRESH
        temp = temp+1;
        sum_RAP = 0; sum_AP = 0;
        for True_nsources = True_numS1
            sum_RAP = sum_RAP + Accuracy{True_nsources+1, temp}(1,2);
            sum_AP = sum_AP + Accuracy{True_nsources+1, temp}(1,3);
        end
        Summ_RAP(temp,1) = sum_RAP;
        Summ_AP(temp,1) = sum_AP;
    end

    [A2,B2] = sort(Summ_RAP,'descend');
    [A3,B3] = sort(Summ_AP,'descend');

    OptimumThresh(True_nsources+1,:) = [THRESH(1,B2(1,1)) THRESH(1,B3(1,1))];
    BestAccuracy(True_nsources+1,:) = [Accuracy{True_nsources+1, B2(1,1)}(1,2) Accuracy{True_nsources+1, B3(1,1)}(1,3)];
end