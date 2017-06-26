function   [fillMod] = padFill(sidepad,omod)% function [fillMod] = padFill(sidepad,omod)%% Fills in the sides of a velocity or loss model with boundary% data. Useful when created a padded zone for a Cerjan et.al. style% sponge implementation. % % Input Arguments% -----------------------------------------------------% 1. sidepad        : width of sidepad in samples (same for%                     top and sides)%% 2. omod           : the original model without sidepads%%  Output Variables% -----------------------------------------------------% 1. fillMod        : padded model%%% By Jonathan Ajo-Franklin (Started 11/14/04)% % [IMPORTANT] % This function and all associated data files included% are free to use and covered under creative commons% attribution/share-alike  2.5 license. This means you can use them% for  anything you like including commercial apps providing that% you  attribute the creation (to me - see below) and you make% derivative  products available under the same license. If you use% them  in a scientific paper or class project please acknowledge% the source  as,    %% [Dr. Jonathan Ajo-Franklin, Lawrence Berkeley National Laboratory]%% For license details see %% http://creativecommons.org/licenses/by-sa/2.5/ xc = size(omod,2);zc = size(omod,1);truex = xc+(sidepad*2);truez = zc+(sidepad*2);fillMod = zeros(truez,truex);% filling in the center with the true modelfillMod(sidepad+1:(zc+sidepad),sidepad+1:(xc+sidepad)) = omod;% filling the vertical sidepads with velocity and loss datafor I=sidepad:(xc+sidepad),    for J=1:sidepad,          fillMod(J,I) = fillMod(sidepad+1,I);      end;    for J=(zc+sidepad):truez, fillMod(J,I) = fillMod(zc+sidepad-1,I);   end;end;% filling the horizontal sidepads with velocity and loss datafor J=1:(zc+(2*sidepad)),    for I=1:sidepad,          fillMod(J,I) = fillMod(J,sidepad+1);      end;    for I=(xc+sidepad):truex, fillMod(J,I) = fillMod(J,xc+sidepad-1);   end;end;