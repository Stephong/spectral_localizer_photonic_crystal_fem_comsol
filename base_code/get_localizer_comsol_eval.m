function localizer_eval = get_localizer_comsol_eval(tab_x_lambd, tab_y_lambd, tab_E_lambd, Iden, X, Y, H, kappa, options)
%

arguments
    tab_x_lambd
    tab_y_lambd
    tab_E_lambd
    Iden
    X
    Y
    H
    kappa
    options.N_spec = 10
    options.DoClean = false
    options.do_L_spec = 1
    options.do_C_L = 1
    options.do_print = 1
end

N_spec = options.N_spec;

assert(issparse(X));
assert(issparse(Y));
assert(issparse(H));


N_x_lambd = length(tab_x_lambd);
N_y_lambd = length(tab_y_lambd);
N_E_lambd = length(tab_E_lambd);

% Iden = speye(size(X,1));

tab_L_spec = zeros(N_E_lambd, N_x_lambd, N_y_lambd, N_spec);
tab_L_spec_flag = zeros(N_E_lambd, N_x_lambd, N_y_lambd);
tab_L_gap = zeros(N_E_lambd, N_x_lambd, N_y_lambd);
tab_C_L = zeros(N_E_lambd, N_x_lambd, N_y_lambd);

if ~options.DoClean
    tic
end

iter = 1;
for i_E=1:N_E_lambd
    for i_x=1:N_x_lambd
        if options.do_print 
            fprintf(repmat('\b', 1, 200));
            disp(iter/(N_E_lambd*N_x_lambd));
        end
        for i_y=1:N_y_lambd

            x = tab_x_lambd(i_x);
            y = tab_y_lambd(i_y);
            E = tab_E_lambd(i_E);
            L = [
                (H - E*Iden), kappa*(X-x*Iden)-1i*kappa*(Y-y*Iden);
                kappa*(X-x*Iden)+1i*kappa*(Y-y*Iden), -(H - E*Iden)
                ];

            if options.do_L_spec
                [~, tab_L_spec_i, tab_L_spec_flag_i] = eigs(L, N_spec, 'sm');
                tab_L_spec_i = diag(tab_L_spec_i);
    
                tab_L_spec(i_E, i_x, i_y, :) = tab_L_spec_i;
                tab_L_spec_flag(i_E, i_x, i_y) = tab_L_spec_flag_i;
                tab_L_gap(i_E, i_x, i_y) = min(abs(tab_L_spec_i));
            end

            if options.do_C_L
                tab_C_L(i_E, i_x, i_y) = 0.5 * signature(L);
            end

        end
        iter = iter + 1;
    end
end

if options.do_print
    fprintf(repmat('\b', 1, 200));
end

if ~options.DoClean
    toc
end

localizer_eval.tab_L_spec = tab_L_spec;
localizer_eval.tab_L_spec_flag = tab_L_spec_flag;
localizer_eval.tab_L_gap = tab_L_gap;
localizer_eval.tab_C_L = tab_C_L;

end