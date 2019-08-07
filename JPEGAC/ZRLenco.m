function out=ZRLenco(scan)
% zero count should not be larger or equal to 16
out=[];
zero=0;
for i=1:length(scan)
    if ~any(scan(i:end))% EOB
        out = [out [0;0]];
        return;
    elseif scan(i)==0 && zero~=15
        zero=zero+1;
    else
        out=[out [zero;scan(i)]];
        zero=0;
    end
end