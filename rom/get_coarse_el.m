function [E] = get_coarse_el(n_f, n_c, e)
    %Takes element number of full order model, gives element number of
    %coarse model
    %n_f = FOM number of elements, n_c = coarse number of elements, e =
    %fine element number
    
    
    %fine elements per coarse mesh
    fine_per_coarse = n_f./n_c;
    %must be integer
    assert(~any(mod(n_f, n_c)), 'Error: no integer number of fine elements within a coarse element')
    
    row_fine = floor((e - 1)/n_f(1) + 1);
    col_fine = mod((e - 1), n_f(1)) + 1;
    
    row_coarse = floor((row_fine - 1)/fine_per_coarse(2) + 1);
    col_coarse = floor((col_fine - 1)/fine_per_coarse(1) + 1);
        
    E = (row_coarse - 1)*n_c(1) + col_coarse;
    
    %checkerboard plots; only possible if vector e contains all elements
    plt = false;
    if plt
        
       el = reshape(e, n_f(1), n_f(2))';
       rf = reshape(row_fine, n_f(1), n_f(2))';
       cf = reshape(col_fine, n_f(1), n_f(2))';
       rc = reshape(row_coarse, n_f(1), n_f(2))';
       cc = reshape(col_coarse, n_f(1), n_f(2))';
       figure;
       subplot(3,2,1)
       pcolor(el)
       subplot(3,2,3)
       pcolor(rf)
       subplot(3,2,4)
       pcolor(cf)
       subplot(3,2,5)
       pcolor(rc)
       subplot(3,2,6)
       pcolor(cc)
        
    end
    
    
end

