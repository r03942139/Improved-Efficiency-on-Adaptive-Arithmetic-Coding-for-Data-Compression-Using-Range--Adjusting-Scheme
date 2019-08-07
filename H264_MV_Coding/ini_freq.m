function P_out = ini_freq(range)

P_out = ones(1,range);
for x = 1:range
	P_out(x) = max(60 * exp(-7*(x-1)),1);
end


end