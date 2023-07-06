function [Accuracy] = Thresholdfunction(MaxnumSources,THRESH,MC_repetitions,Fratio_AP)
for nsources=0:MaxnumSources
    temp = 0;
    for Threshold = THRESH
        temp = temp+1;
        true1=0; false1=0;
        tee = 6;
        for k = nsources:MaxnumSources
            if k == nsources
                for j = 0:nsources-1
                    if Fratio_AP{k+1,idxMC}(j+1,1)> Threshold
                        true1 = true1+1;
                    else
                        false1 = false1+1;
                    end
                end
                for idxMC= 1:MC_repetitions
                    if Fratio_AP{k+1,idxMC}(nsources+1,1) < Threshold
                        true1 = true1+1;
                    else
                        false1 = false1+1;
                    end
                end
            else
                for idxMC= 1:MC_repetitions
                    if Fratio_AP{k+1,idxMC}(nsources+1,1) > Threshold
                        true1 = true1+1;
                    else
                        false1 = false1+1;
                    end
                end
            end
        end
        Accuracy(nsources+1,temp) = true1./tee; false11(nsources+1,temp) = false1;
    end
end












% for True_nsources=0:MaxnumSources
%     temp = 0;
%     for Threshold = THRESH
%         temp = temp+1;
%         true1=0; false1=0;
%         tee = 0;
%         for k = True_nsources:MaxnumSources
%             tee = tee+1;
%             if k == True_nsources
%                 for idxMC= 1:MC_repetitions
%                     if Fratio_AP{k+1,idxMC}(True_nsources+1,1) < Threshold
%                         true1 = true1+1;
%                     else
%                         false1 = false1+1;
%                     end
%                 end
%             else
%                 for idxMC= 1:MC_repetitions
%                     if Fratio_AP{k+1,idxMC}(True_nsources+1,1) > Threshold
%                         true1 = true1+1;
%                     else
%                         false1 = false1+1;
%                     end
%                 end
%             end
%         end
%         Accuracy(True_nsources+1,temp) = true1./tee; false11(True_nsources+1,temp) = false1;
%     end
% end





% for True_nsources=0:MaxnumSources
%     temp = 0;
%     for Threshold = THRESH
%         temp = temp+1;
%         true1=0; false1=0;
%         tee = 0;
%         for k = True_nsources:MaxnumSources
%             tee = tee+1;
%             if k == True_nsources
%                 for idxMC= 1:MC_repetitions
%                     if Fratio_AP{k+1,idxMC}(True_nsources+1,1) < Threshold
%                         true1 = true1+1;
%                     else
%                         false1 = false1+1;
%                     end
%                 end
%             else
%                 for idxMC= 1:MC_repetitions
%                     if Fratio_AP{k+1,idxMC}(True_nsources+1,1) > Threshold
%                         true1 = true1+1;
%                     else
%                         false1 = false1+1;
%                     end
%                 end
%             end
%         end
%         Accuracy(True_nsources+1,temp) = true1./tee; false11(True_nsources+1,temp) = false1;
%     end
% end
%
%




