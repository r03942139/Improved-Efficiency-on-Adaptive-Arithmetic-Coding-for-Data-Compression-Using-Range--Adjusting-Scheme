function [top,bot]=divRange(top0,bot0,C,sym)
    % sym:[0,1,2...]

    bot=bot0+(top0-bot0)*C(sym);
    top=bot0+(top0-bot0)*C(sym+1);
end