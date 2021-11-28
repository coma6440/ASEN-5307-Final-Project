function out = ang_mean(x)
%Compute anagular mean by mean of sin and cos components
out = atan2(mean(sin(x)),mean(cos(x)));
end