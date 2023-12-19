
function t_max_final = automate_cell_direct(obj,cond_haute,taux_remplissage,start_image)
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

disp('Reading bitmap image...')

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

disp('Converting bitmap image to surface conditions...');

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

for m=1:1:nombre_images;
    tic
    disp(['-------------------------------------------------------'])
    disp('Applying the Cellular Automaton...');
    
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

    note=gradients*(1-obj)+temp*(obj);
    [condi_limites] = cellular_automaton(condi_limites,note, cond_haute,cond_basse,nombre_pixels_conducteurs,taux_variac,nombre_images,m);
    disp('Calculating the temperature map...');
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
    miroir=fliplr(condi_limites_2(1:hauteur,1:largeur-1,:));
    miroir2=fliplr(miroir);
    imwrite([miroir2,miroir],['Output_kp_ko_',num2str(cond_haute),'_phi_',num2str(taux_remplissage),'.png']);
    imwrite([miroir2,miroir],['Topology/sortie_kp_ko_',num2str(cond_haute),'_phi_',num2str(taux_remplissage),'_',num2str(m,'%06.f'),'.png']);
    figure(1)
    colormap jet
    subplot(2,4,1:2);
    
    if m>1;
        res(m-1)=(t_max_sortie(m-1)-t_max_sortie(m))/(t_max_sortie(1)-t_max_sortie(2));
        semilogy(abs(res),'.r');
        title('Residuals');
    end
    
    subplot(2,4,3:4);
    imagesc([miroir2,miroir]);
    title('Topology');
    
    subplot(2,4,5);
    plot(variance,'.m');
    imagesc(log10(entropie(2:end-1,2:end-1)));
    title('Log10 Entropy');
    
    subplot(2,4,6);
    imagesc(temp);
    title('Temperature');
    
    subplot(2,4,7);
    imagesc(log10(gradients));
    title('Log10 Thermal gradients');
    
    subplot(2,4,8);
    imagesc(log10(note));
    title('Log10 Objective function');
    
    disp(['Maximum temperature: ',num2str(t_max_sortie(m))])
    disp(['End of epoch: ',num2str(m)])
    disp(['-------------------------------------------------------'])
    
    saveas(gcf,['Figure_kp_ko_',num2str(cond_haute),'_phi_',num2str(taux_remplissage),'.png']);
    saveas(gcf,['Figure/figure_kp_ko_',num2str(cond_haute),'_phi_',num2str(taux_remplissage),'_',num2str(m,'%06.f'),'.png']);
    toc
    
    % figure(2)
    % plot(index_de_m,'bd');
end;
t_max_final=max(max(temp));
disp('End of convergence !');
