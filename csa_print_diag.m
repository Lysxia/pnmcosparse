function [] = csa_print_diag( D,name )

    colormap(gray(100));
    image(100*D);
    drawnow;
    
    if ~strcmp(name,'')
        print('-dpng',strcat('diagrams/',name,'.png'));
    end

end

