load DATA

st_time = cputime
fprintf('start at %d seconds\n', st_time);
W = SOBI(DATA, [0, 1:10, 12:2:20, 25:5:100, 125:25:350], 1000, 1e-5);
ed_time = cputime
fprintf('end at %d seconds\n', ed_time);
fprintf('duration = %d seconds\n', ed_time - st_time);
%core code of SOBI function comes from Jean-Francoi s Cardoso, with
%modifications by various individuals worked in Prof Akaysha Tang's labs
