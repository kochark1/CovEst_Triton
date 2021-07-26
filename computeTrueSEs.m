function computeTrueSEs(path, interm_folder,out_folder)
    interm_folder=strcat(path, '/', interm_folder);
    out_folder=strcat(path, '/', out_folder);
    load(strcat(interm_folder,'/RQWfile.mat'));
    
    trueAndTheoretical;
end