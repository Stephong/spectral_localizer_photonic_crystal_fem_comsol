function [] = init_comsol_server(options)

arguments
    % options.comsol_server_dir = 'C:\Program Files\COMSOL\COMSOL61\Multiphysics\bin\win64';
    % options.comsol_connection_dir = 'C:\Program Files\COMSOL\COMSOL61\Multiphysics\mli';
    options.comsol_server_dir = 'C:\Program Files\COMSOL\COMSOL62\Multiphysics\bin\win64';
    options.comsol_connection_dir = 'C:\Program Files\COMSOL\COMSOL62\Multiphysics\mli';
end

% from https://www.youtube.com/watch?v=JBYKSigJvUo

current_dir = pwd;

% launch comsol server
cd(options.comsol_server_dir);
system('comsolmphserver.exe &');
cd(current_dir);

% establish connection
cd(options.comsol_connection_dir);
mphstart(2036);
cd(current_dir);

end




