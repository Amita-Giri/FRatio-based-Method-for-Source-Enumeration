clear
close all
clc

%%
load('Result_Real300_330OSVol.mat')
S_Reduced{6,1} = S_Full{5,1};
qq_Reduced{6,1} = qq_Full{5,1};
S_Est_Reduced{6,1} = S_Est_Full{5,1};
%% Plot dipoles
for run = 1:4
    temps = S_Reduced{run,1};
    tempq = qq_Reduced{run,1};
    sig = sign(S_Est_Reduced{run,1});
    ab = length(temps);
    for kkk=1:ab
        origin = GridLoc(temps(1,kkk),:);
        Orient = sig(kkk,1)*tempq(3*(kkk-1)+1:3*kkk,kkk);
        [azimuth,elevation,r] = cart2sph(Orient(1,1),Orient(2,1),Orient(3,1));
        azimuth = azimuth*180/pi
        elevation = elevation*180/pi
        alpha=azimuth; beta=0; gamma=elevation;

        %% sphere & Cylinder Parameters
        rs = 0.01;
        npoints = 20;
        rc = 0.003; hc = 0.03;

        [Xs,Ys,Zs] = sphere(npoints);
        Xs = Xs * rs;
        Ys = Ys * rs;
        Zs = Zs * rs;

        %% Cylinder
        [Xc,Yc,Zc] = cylinder(rc);
        Zc = Zc*hc;

        %% Together
        RX = rotx(alpha);
        RY = roty(beta);
        RZ = rotz(gamma);
        R = RZ*RY*RX;

        for i=1:npoints+1
            for j= 1:npoints+1
                temp = R*[Xs(i,j); Ys(i,j); Zs(i,j)];
                Xs1(i,j) = temp(1,1)+origin(1,1);
                Ys1(i,j) = temp(2,1)+origin(1,2);
                Zs1(i,j) = temp(3,1)+origin(1,3);
            end
        end
        AS{kkk,1} = Xs1; AS{kkk,2} = Ys1; AS{kkk,3} = Zs1;

        for i=1:2
            for j =1:npoints+1
                temp = R*[Xc(i,j); Yc(i,j); Zc(i,j)];
                Xc1(i,j) = temp(1,1)+origin(1,1);
                Yc1(i,j) = temp(2,1)+origin(1,2);
                Zc1(i,j) = temp(3,1)+origin(1,3);
            end
        end
        AC{kkk,1} = Xc1; AC{kkk,2} = Yc1; AC{kkk,3} = Zc1;
    end

    openfig("cortex.fig")
    for kkk=1:ab
        if kkk==1
            hold on
            surf(AS{kkk,1},AS{kkk,2},AS{kkk,3}, 'edgecolor', 'none','facecolor','b')
            hold on
            surf(AC{kkk,1},AC{kkk,2},AC{kkk,3}, 'edgecolor', 'none','facecolor','b')
            hold on
        elseif kkk==2
            hold on
            surf(AS{kkk,1},AS{kkk,2},AS{kkk,3}, 'edgecolor', 'none','facecolor','r')
            hold on
            surf(AC{kkk,1},AC{kkk,2},AC{kkk,3}, 'edgecolor', 'none','facecolor','r')
            hold on
        else
            hold on
            surf(AS{kkk,1},AS{kkk,2},AS{kkk,3}, 'edgecolor', 'none','facecolor','g')
            hold on
            surf(AC{kkk,1},AC{kkk,2},AC{kkk,3}, 'edgecolor', 'none','facecolor','g')
            hold on
        end
    end
    hold off
    view(180,90)
end









