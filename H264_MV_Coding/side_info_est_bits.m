function [ num_bits ] = side_info_est_bits( table )
    num_bits = 0;
    
    % Error handling
    assert(size(table,1)==1, "Error: Table is not 1D \n");
    
    for i = 1:length(table)
        bits_needed = length(dec2bin(table(i))) + 1;
        % Indicator of how many bits
        num_bits = num_bits + length(dec2bin(bits_needed));
        % Freq count recording: ranging from 
        num_bits = num_bits + bits_needed;
    end


end

