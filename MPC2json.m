% Takes a loaded mpc instance and creates a file that can be read by the
% OPF_Tools library

function res = MPC2json(mpc, filename)


fileID = fopen(filename, 'w');

text = jsonencode(mpc);
fprintf(fileID, text);

fclose(fileID);

