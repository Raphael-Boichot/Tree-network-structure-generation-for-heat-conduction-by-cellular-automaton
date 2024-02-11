%https://doi.org/10.1016/j.enconman.2008.09.003
%https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton
function [boudary_conditions] = cellular_automaton(boudary_conditions,temp_grad_ratio, kp_k0,k0,conductive_cells,variation_rate,number_of_epoch,m)
[height,width,~]=size(boudary_conditions);
new_pos_out=zeros(height,width);
new_pos_in=zeros(height,width);
new_pos2=zeros(height,width);
%************calcul des gradients de température et des gradients
%contigus au conducteur*****
for k=2:1:(height-1)
    for l=2:1:(width-1)
        new_pos_out(k, l) = 0;
        new_pos_in(k, l) = 0;
        % calcul des positions autorisées pour les nouvelles cellules
        %************à l'extérieur du conducteur*********
        anti_warning = 0;
        if boudary_conditions(k, l) == k0
            if (boudary_conditions(k - 1, l) == kp_k0) && (boudary_conditions(k - 1, l - 1) == kp_k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k - 1, l) == kp_k0) && (boudary_conditions(k - 1, l + 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (boudary_conditions(k + 1, l) == kp_k0) && (boudary_conditions(k + 1, l - 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (boudary_conditions(k + 1, l) == kp_k0) && (boudary_conditions(k + 1, l + 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l - 1) == kp_k0) && (boudary_conditions(k - 1, l - 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l - 1) == kp_k0) && (boudary_conditions(k + 1, l - 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l + 1) == kp_k0) && (boudary_conditions(k - 1, l + 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l + 1) == kp_k0) && (boudary_conditions(k + 1, l + 1) == kp_k0) ; anti_warning = anti_warning + 1; end
            if (anti_warning >= 1); new_pos_out(k, l) = temp_grad_ratio(k, l); end
            %*****************enlève les cases isolées en priorité********************************************
            if (anti_warning == 8); new_pos_out(k, l) = temp_grad_ratio(k, l) * 1000; end
        end
        %*************puis à l'intérieur du conducteur*******
        anti_warning = 0;
        if boudary_conditions(k, l) == kp_k0
            if (boudary_conditions(k - 1, l) == k0) && (boudary_conditions(k - 1, l - 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k - 1, l) == k0) && (boudary_conditions(k - 1, l + 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k + 1, l) == k0) && (boudary_conditions(k + 1, l - 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k + 1, l) == k0) && (boudary_conditions(k + 1, l + 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l - 1) == k0) && (boudary_conditions(k - 1, l - 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l - 1) == k0) && (boudary_conditions(k + 1, l - 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l + 1) == k0) && (boudary_conditions(k - 1, l + 1) == k0); anti_warning = anti_warning + 1; end
            if (boudary_conditions(k, l + 1) == k0) && (boudary_conditions(k + 1, l + 1) == k0); anti_warning = anti_warning + 1; end
            if (anti_warning >= 1); new_pos_in(k, l) = temp_grad_ratio(k, l); end
            %*****************remplit les cases isolées en priorité*******************************************
            if (anti_warning == 8); new_pos_in(k, l) = temp_grad_ratio(k, l) / 1000; end
        end
    end
end

%************calcul des n meilleurs et plus mauvais gradients************************************
for p = 1:1:max(uint32(conductive_cells*variation_rate*(number_of_epoch-m)/(number_of_epoch)),1)
    high = 0;
    low = 1E+15;
    for k=2:1:(height-1)
        for l=2:1:(width-1)
            %Recherche le plus forts et plus faible gradient dans et hors des drains thermiques
            if (new_pos_out(k, l) > high)
                high = new_pos_out(k, l);
                pos_high_x = k;
                pos_high_y = l;
            end
            if (new_pos_in(k, l) ~= 0) && (new_pos_in(k, l) < low)
                low = new_pos_in(k, l);
                pos_low_x = k;
                pos_low_y = l;
            end
        end
    end
    %****donne les nouvelles coordonées des deux prochains déplacements de cellules
    new_pos2(pos_high_x, pos_high_y) = 1;
    new_pos2(pos_low_x, pos_low_y) = -1;
    new_pos_out(pos_high_x, pos_high_y) = 0;
    new_pos_in(pos_low_x, pos_low_y) = 0;
end
%************Changement des conditions limites*****************************
for k=2:1:(height-1)
    for l=2:1:(width-1)
        if k>height;disp(k);end
        if l>width;disp(l);end
        if new_pos2(k, l) == -1; boudary_conditions(k, l) = k0; end
        if new_pos2(k, l) == 1; boudary_conditions(k, l) = kp_k0; end
    end
end
%******************************************Fin de l'automate cellulaire****
end

