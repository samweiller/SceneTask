%blah test
%% sceneTaskTimingAnalysis
%  This compiles all fixation and block timings into two matricies for
%  averaging.


subs2analyze = 9:15;

k = 1;
blockTime = [];
fixTime = [];
for i = 1:max(size(subs2analyze))
    for run = 1:2
        load(sprintf('SCTASK_S%d_C%d_A%d.mat', subs2analyze(i), run, run));
        
        
        for b = 1:25
            switch mod(b,2)
                case 0
                    blockTime(k, ceil(b/2)) = timeLogger.block(b).blockLength;
                case 1
                    fixTime(k, ceil(b/2)) = timeLogger.block(b).blockLength;
            end;
        end;
        
        k = k+1;
    end;
end;