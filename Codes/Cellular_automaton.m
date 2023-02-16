function [condi_limites] = Cellular_automaton(condi_limites,note, cond_haute,cond_basse,nombre_pixels_conducteurs,taux_variac,nombre_images,m)
        
[hauteur,largeur,profondeur]=size(condi_limites);
new_pos_out=zeros(hauteur,largeur);
new_pos_in=zeros(hauteur,largeur);
new_pos2=zeros(hauteur,largeur);
%************calcul des gradients de température et des gradients
        %contigus au conducteur*****
        for k=2:1:(hauteur-1);
        for l=2:1:(largeur-1);
        
                new_pos_out(k, l) = 0;
                new_pos_in(k, l) = 0;
                
                % calcul des positions autorisées pour les nouvelles cellules
                %************à l'extérieur du conducteur*********
                
                anti_warning = 0;
                if condi_limites(k, l) == cond_basse;
                    if (condi_limites(k - 1, l) == cond_haute) && (condi_limites(k - 1, l - 1) == cond_haute); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k - 1, l) == cond_haute) && (condi_limites(k - 1, l + 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                    if (condi_limites(k + 1, l) == cond_haute) && (condi_limites(k + 1, l - 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                    if (condi_limites(k + 1, l) == cond_haute) && (condi_limites(k + 1, l + 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l - 1) == cond_haute) && (condi_limites(k - 1, l - 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l - 1) == cond_haute) && (condi_limites(k + 1, l - 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l + 1) == cond_haute) && (condi_limites(k - 1, l + 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l + 1) == cond_haute) && (condi_limites(k + 1, l + 1) == cond_haute) ; anti_warning = anti_warning + 1; end;
                
                    if (anti_warning >= 1); new_pos_out(k, l) = note(k, l); end;
                    
%*****************enlève les cases isolées en priorité********************************************
                    if (anti_warning == 8); new_pos_out(k, l) = note(k, l) * 1000; end;
                                        
                end;
                                
                %*************puis à l'intérieur du conducteur*******
                
                anti_warning = 0;
                if condi_limites(k, l) == cond_haute;
                    if (condi_limites(k - 1, l) == cond_basse) && (condi_limites(k - 1, l - 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k - 1, l) == cond_basse) && (condi_limites(k - 1, l + 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k + 1, l) == cond_basse) && (condi_limites(k + 1, l - 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k + 1, l) == cond_basse) && (condi_limites(k + 1, l + 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l - 1) == cond_basse) && (condi_limites(k - 1, l - 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l - 1) == cond_basse) && (condi_limites(k + 1, l - 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l + 1) == cond_basse) && (condi_limites(k - 1, l + 1) == cond_basse); anti_warning = anti_warning + 1; end;
                    if (condi_limites(k, l + 1) == cond_basse) && (condi_limites(k + 1, l + 1) == cond_basse); anti_warning = anti_warning + 1; end;
                
                if (anti_warning >= 1); new_pos_in(k, l) = note(k, l); end;
                
%*****************remplit les cases isolées en priorité*******************************************
                if (anti_warning == 8); new_pos_in(k, l) = note(k, l) / 1000; end;
                
                end;
       
        end;
        end;
        
%************calcul des n meilleurs et plus mauvais gradients************************************
            for p = 1:1:max(uint32(nombre_pixels_conducteurs*taux_variac*(nombre_images-m)/(nombre_images)),1);

                high = 0;
                low = 1E+15;

        for k=2:1:(hauteur-1);
        for l=2:1:(largeur-1);
            
%Recherche le plus forts et plus faible gradient dans et hors des drains thermiques     
                    if (new_pos_out(k, l) > high);
                        high = new_pos_out(k, l);
                        pos_high_x = k;
                        pos_high_y = l;
                    end;
                
                    if (new_pos_in(k, l) ~= 0) && (new_pos_in(k, l) < low);
                        low = new_pos_in(k, l);
                        pos_low_x = k;
                        pos_low_y = l;
                    end;
            
        end;
        end;
        
%****donne les nouvelles coordonées des deux prochains déplacements de cellules
                new_pos2(pos_high_x, pos_high_y) = 1;
                new_pos2(pos_low_x, pos_low_y) = -1;
                new_pos_out(pos_high_x, pos_high_y) = 0;
                new_pos_in(pos_low_x, pos_low_y) = 0;
                
            end;
        
%************Changement des conditions limites*****************************  
        for k=2:1:(hauteur-1);
        for l=2:1:(largeur-1);
        
            if k>hauteur;disp(k);end;
            if l>largeur;disp(l);end;
            
            if new_pos2(k, l) == -1; condi_limites(k, l) = cond_basse; end;             
            if new_pos2(k, l) == 1; condi_limites(k, l) = cond_haute; end;
                        
        end;
        end;

%******************************************Fin de l'automate cellulaire****
end

