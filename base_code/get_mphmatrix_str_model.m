function [mphmatrix_str, mphxmeshinfo_str] = get_mphmatrix_str_model(model, options)

arguments
    model, 
    options.do_plot = 0
end

std1 = model.study('std1');
std1.run;

% get mphmatrix_str
mphmatrix_str = mphmatrix(model, 'sol1', ...
    out={'K', 'L', 'M', 'N', 'D', 'E', 'NF', ...
         'NP', 'MP', 'MLB', 'MUB', ...
         'Kc', 'Lc', 'Dc', 'Ec', 'Null', 'Nullf', 'ud', ...
         'uscale'});

% getmphxmeshinfo_str
mphxmeshinfo_str = mphxmeshinfo(model, soltag='sol1');

if options.do_plot
    % plot geometry
    figure
    mphgeom(model)

    % plot mesh
    figure
    mphmesh(model)
end

end