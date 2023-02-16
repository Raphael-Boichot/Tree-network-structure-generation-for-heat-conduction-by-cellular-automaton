
function t_max_final = automate_cell_direct(cond_haute,taux_remplissage,start_image)
rng('shuffle', 'twister')
%*****************Automate cellulaire*********************INPG/BOICHOT/2008
%****Conductivité des drains thermiques W/m.K******************************
% cond_haute = 10;
cond_basse = 1;
temp_puits = 300;
pas_x = 0.001;
p_vol=1e5;
taux_variac = 0.5;
% taux_remplissage=0.3;
%****Image à traiter comme configuration initiale**************************
condi_limites_1=imread(start_image);

disp('lecture de l''image terminée..........................................')

%****Récupération du format de l'image*************************************
[hauteur,largeur,profondeur]=size(condi_limites_1);

nombre_images = max([hauteur,largeur]);
condi_limites = zeros(hauteur,largeur);
condi_temp = zeros(hauteur,largeur);
f_temp=zeros(nombre_images);
f_iter=1:1:nombre_images;



%****Assignation des conditions aux limites en fonction de la couleur des
%pixels de l'image*********************************************************
pixels_blancs=0;
pixels_noirs=0;
for k = 1:1:hauteur;
    for l = 1:1:largeur; 
    rouge = condi_limites_1(k,l,1);
    vert = condi_limites_1(k,l,2);
    bleu = condi_limites_1(k,l,3);
    
        if (rouge == 255) && (vert == 255) && (bleu == 255); 
            choix = cond_basse;
            pixels_blancs=pixels_blancs+1;
        end;
        if (rouge == 127) && (vert == 127) && (bleu == 127); 
            choix = -2;
        end;
        if (rouge == 0) && (vert == 0) && (bleu == 255); 
            choix = -3;
        end;
        if (rouge == 0) && (vert == 0) && (bleu == 0); 
            choix = cond_haute;
            pixels_noirs=pixels_noirs+1;
        end;
                
        condi_limites(k, l) = choix;
       
     end;
end;

nombre_pixels_conducteurs=ceil(pixels_blancs*taux_remplissage);
condi_limites=init_image(condi_limites,nombre_pixels_conducteurs, cond_basse, cond_haute);

disp('conversion des conditions limites terminée...........................');

%****Pré-allocation de la taille des matrices utilisées dans les boucles***
temp=ones(hauteur,largeur).*temp_puits;
condu_tab=zeros(hauteur,largeur,4);
new_pos_in=zeros(hauteur,largeur);
new_pos_out=zeros(hauteur,largeur);
new_pos2=zeros(hauteur,largeur);
gradients=zeros(hauteur,largeur);
note=zeros(hauteur,largeur);
condi_limites_2=condi_limites_1;
affichage=zeros(1,4);

    
%disp('entrée des conditions initiales terminée.............................');

for m=1:1:nombre_images*1.2;
tic

disp('Echange de position des cellules.....................................');

%************************************************************Début de l'automate cellulaire


%boucle interne
%************Calcul de temp_max, temp_min, grad_max,grad_min*****
[somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max_sortie(m),temp,grad, variance_grad]=finite_temp_direct_sparse(cond_haute,cond_basse,temp_puits,pas_x,p_vol,condi_limites);
gradients=zeros(hauteur,largeur);
gradients(2:hauteur-1,2:largeur-1)=grad;
grad_max = max(max(gradients));
grad_min = min(min(gradients));
temp_max = max(max(temp));
temp_min = min(min(temp));
gradients=(gradients-grad_min)/(grad_max-grad_min);
temp=(temp-temp_min)/(temp_max-temp_min);
% 
% parfor objvar=1:1:21;
% [somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max(objvar),temp_temp,grad, variance_grad]=finite_temp_direct_sparse(cond_haute,cond_basse,temp_puits,pas_x,p_vol,Cellular_automaton(condi_limites,gradients*(1-(objvar-1)/20)+temp*((objvar-1)/20), cond_haute,cond_basse,nombre_pixels_conducteurs,taux_variac,nombre_images,m));
% end
% 
% [val,index]=min(t_max);
% index_de_m(m)=index;
% obj= (index-1)/20;  
% note = gradients*(1-obj)+temp*(obj);
note=gradients;
[condi_limites] = cellular_automaton(condi_limites,note, cond_haute,cond_basse,nombre_pixels_conducteurs,taux_variac,nombre_images,m);
[somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max,temp,grad, variance_grad]=finite_temp_direct_sparse(cond_haute,cond_basse,temp_puits,pas_x,p_vol,condi_limites);
%****créé une image de sortie compatible avec l'image d'entrée*************
for k = 1:1:hauteur;
   for l = 1:1:largeur; 
       
    choix = condi_limites(k, l) ;
     
        if choix == cond_basse;
            rouge = 255;
            vert = 255;
            bleu = 255; 
        end;
        
        if choix == -2;
            rouge = 127; 
            vert = 127 ;
            bleu = 127 ;
        end;
        if choix == -3;
            rouge = 0;
            vert = 0 ;
            bleu = 255; 
        end;
        if choix == cond_haute;
            rouge = 0 ;
            vert = 0 ;
            bleu = 0;
        end;

    condi_limites_2(k,l,1)=rouge;
    condi_limites_2(k,l,2)=vert;
    condi_limites_2(k,l,3)=bleu;
     
   end;
end;


condi_limites_2=uint8(condi_limites_2);

% nom_sortie = ['sortie',num2str(m),'.bmp'];
nom_sortie =['sortie_kp_ko_',num2str(cond_haute),'_phi_',num2str(taux_remplissage),'.png'];
miroir=fliplr(condi_limites_2(1:hauteur,1:largeur-1,:));
miroir2=fliplr(miroir);
imwrite([miroir2,miroir],nom_sortie);
figure(1)
subplot(2,4,1:2);
plot(t_max_sortie,'.r');
title('Max temperature');

subplot(2,4,3:4);
imagesc([miroir2,miroir]);
title('Géométrie');

subplot(2,4,5);
plot(variance,'.m');
imagesc(log10(entropie(2:end-1,2:end-1)));
colormap hot
title('Log10 Entropie');

subplot(2,4,6);
imagesc(temp);
colormap hot
title('Température');

subplot(2,4,7);
imagesc(gradients);
colormap hot
title('Gradients thermiques');

subplot(2,4,8);
imagesc(note);
title('Note pondérée');
colormap hot

t_max_final=t_max_sortie(m)
condi_limites_1=condi_limites_2;
pause(0.01);
%saveas(gcf,['Z_figure_', num2str(m)],'png');
saveas(gcf,['figure_kp_ko_',num2str(cond_haute),'_phi_',num2str(taux_remplissage),'.png']);
toc

% figure(2)
% plot(index_de_m,'bd');
end;
disp('Calcul terminé !............................................');
