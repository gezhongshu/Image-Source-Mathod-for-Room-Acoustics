if overlap_only == 1 || overlap_only == 0
    [ma,R]=max(x); % Find the direct soudn power, which is the lasrgest one
    R2 = find(IR~=0); % Find the first cell that has value
    di = R-R2(1); % find the position difference between the model and real IR
    figure;
    if overlap_only == 1
        % match the direct sound position
        if di > 0
            t = 1:TSample;
            x = x(di:TSample+di);
            plot((t)/Fs,x(t))
            hold on
            plot(TimePoints/Fs,IR,'LineWidth',2)
            grid;
            hold off
        else
            t = 1:TSample;
            x = [zeros(-di,1);x(1:TSample+di)];
            plot((t)/Fs,x(t))
            hold on
            plot(TimePoints/Fs,IR,'LineWidth',2)
            grid;
            hold off
        end
    % plot three graphs including overlapped echograph, real IR and model
    % IR
    else if overlap_only == 0
            % match the direct sound position
            if di > 0
                subplot(3,1,1)
                t = 1:TSample;
                x = x(di:TSample+di);
                plot((t)/Fs,x(t))
                hold on
                plot(TimePoints/Fs,IR,'LineWidth',2)
                grid;
                hold off
            else
                subplot(3,1,1)
                t = 1:TSample;
                x = [zeros(-di,1);x(1:TSample+di)];
                plot((t)/Fs,x(t))
                hold on
                plot(TimePoints/Fs,IR,'LineWidth',2)
                grid;
                hold off
            end
            subplot(3,1,2)
            plot((t)/Fs,x(t))
            subplot(3,1,3)
            plot(TimePoints/Fs,IR,'LineWidth',2)
        end
    end
end
