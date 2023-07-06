clear
close all
clc

load('Result_Real300_330OSVol.mat')
Time_Range = 101:500;
%% SEE THE PROJECT AWAY Residuals PLOTS
y(1,:) = sum((Multfactor*Data(MEG_Range,Time_Range)).^2);
for j=1:length(S_Full)
    gain_idxF = [];
    for i = 1:length(S_Full{j,1})
        gain_idxF= [gain_idxF (S_Full{j,1}(1,i)-1)*3+1:S_Full{j,1}(1,i)*3];
    end
    q_Full = qq_Full{j,1};
    A_Full = Gain(:,gain_idxF)*q_Full; S_Est_Fullw = (A_Full'*A_Full)\(A_Full'*Data(MEG_Range,Time_Range));
    ab{j,1} = (Data(MEG_Range,Time_Range)-(A_Full*S_Est_Fullw));
    y(j+1,:) = sum((Multfactor*((ab{j,1}))).^2);
end

figure;
x = linspace(-100,300,400);
plot(x,y(1,:),'k','LineWidth', 2)
hold on
plot(x,y(2,:),'b','LineWidth', 2)
hold on;
plot(x,y(3,:),'r','LineWidth', 2)
hold on
plot(x,y(4,:),'g','LineWidth', 2)
hold off
legend('Original','1 Source Model','2 Source Model','3 Source Model')
ylabel('Sum of Squares of Residuals')
xlabel('Time')




