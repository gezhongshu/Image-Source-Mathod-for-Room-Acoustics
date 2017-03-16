function [spk_coef] = direc_judge(ref_direc,direc)
% Judge the energy scalar according to poistions of the mic and the speaker
% ref_drc = rcv-src;
ang = acos(dot(ref_direc,direc)/(norm(ref_direc)*norm(direc)));
spk_coef = 1-(1-0.1)*ang/pi;