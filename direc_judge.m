function [spk_coef] = direc_judge(src,rcv,direc)
% Judge the energy scalar according to 
ref_drc = rcv-src;
ang = acos(dot(ref_drc,direc)/(norm(ref_drc)*norm(direc)));
spk_coef = 1-(1-0.1)*ang/pi;