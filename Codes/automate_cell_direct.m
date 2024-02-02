
function t_max_final = automate_cell_direct(obj,kp_k0,taux_remplissage,start_image)
rng('shuffle', 'twister')
%*****************Automate cellulaire*********************INPG/BOICHOT/2008
%****Conductivité des drains thermiques W/m.K******************************
temp_puits = 300;
pas_x = 0.001;
p_vol=1e5;
taux_variac = 0.5;

%****Image à traiter comme configuration initiale**************************
condi_limites_1=imread(start_image);

disp('Reading bitmap image...')

%****Récupération du format de l'image*************************************
[hauteur,largeur,~]=size(condi_limites_1);

nombre_images = max([hauteur,largeur]);
condi_limites = zeros(hauteur,largeur);

%first pass: count white pixels (non conductive)
pixels_blancs=0;
for k = 1:1:hauteur
    for l = 1:1:largeur
        vec = condi_limites_1(k,l,:);
        if vec==[255 255 255]
            pixels_blancs=pixels_blancs+1;
        end
    end
end

%second pass: fill the image with conductive pixels
nombre_pixels_conducteurs=ceil(pixels_blancs*taux_remplissage);
k=0;
while k<nombre_pixels_conducteurs
    h=ceil(rand*hauteur);
    l=ceil(rand*largeur);
    if condi_limites_1(h,l,:)==[255 255 255]
    condi_limites_1(h,l,:)=[0 0 0]; 
    k=k+1;
    end
end

%save the intial configuration
miroir=fliplr(condi_limites_1(1:hauteur,1:largeur-1,:));
miroir2=fliplr(miroir);
imwrite([miroir2,miroir],['Topology/sortie_kp_ko_',num2str(kp_k0),'_phi_',num2str(taux_remplissage),'_',num2str(0,'%06.f'),'.png']);

%Third pass: convert image into boundary conditions
pixels_blancs=0;
pixels_noirs=0;
for k = 1:1:hauteur
    for l = 1:1:largeur
        rouge = condi_limites_1(k,l,1);
        vert = condi_limites_1(k,l,2);
        bleu = condi_limites_1(k,l,3);
        
        if (rouge == 255) && (vert == 255) && (bleu == 255)
            choix = 1;%by default k0=1, change value here if necessary
            pixels_blancs=pixels_blancs+1;
        end
        if (rouge == 127) && (vert == 127) && (bleu == 127)
            choix = -2;
        end
        if (rouge == 0) && (vert == 0) && (bleu == 255)
            choix = -3;
        end
        if (rouge == 0) && (vert == 0) && (bleu == 0)
            choix = kp_k0;%by default only the ratio is provided
            pixels_noirs=pixels_noirs+1;
        end
        condi_limites(k, l) = choix;
    end
end

disp('Converting bitmap image to surface conditions...');

%****Pré-allocation de la taille des matrices utilisées dans les boucles***
temp=ones(hauteur,largeur).*temp_puits;
condi_limites_2=condi_limites_1;
t_max_sortie=zeros(nombre_images);
res=zeros(nombre_images);
%disp('entrée des conditions initiales terminée.............................');

for m=1:1:nombre_images
    tic
    disp('-------------------------------------------------------')
    disp('Applying the Cellular Automaton...');
    
    %************************************************************Début de l'automate cellulaire
    
    %boucle interne
    %************Calcul de temp_max, temp_min, grad_max,grad_min*****
    [~, ~, ~,~, ~,t_max_sortie(m),temp,grad, ~]=finite_temp_direct_sparse(kp_k0,1,temp_puits,pas_x,p_vol,condi_limites);
    gradients=zeros(hauteur,largeur);
    gradients(2:hauteur-1,2:largeur-1)=grad;
    grad_max = max(max(gradients));
    grad_min = min(min(gradients));
    temp_max = max(max(temp));
    temp_min = min(min(temp));
    gradients=(gradients-grad_min)/(grad_max-grad_min);
    temp=(temp-temp_min)/(temp_max-temp_min);

    note=gradients*(1-obj)+temp*(obj);
    [condi_limites] = cellular_automaton(condi_limites,note, kp_k0,1,nombre_pixels_conducteurs,taux_variac,nombre_images,m);
    disp('Calculating the temperature map...');
    [~, entropie, ~,variance, ~,~,temp,~, ~]=finite_temp_direct_sparse(kp_k0,1,temp_puits,pas_x,p_vol,condi_limites);
    %****créé une image de sortie compatible avec l'image d'entrée*************
    for k = 1:1:hauteur
        for l = 1:1:largeur
            
            choix = condi_limites(k, l) ;
            
            if choix == 1
                rouge = 255;
                vert = 255;
                bleu = 255;
            end
            
            if choix == -2
                rouge = 127;
                vert = 127 ;
                bleu = 127 ;
            end
            if choix == -3
                rouge = 0;
                vert = 0 ;
                bleu = 255;
            end
            if choix == kp_k0
                rouge = 0 ;
                vert = 0 ;
                bleu = 0;
            end
            
            condi_limites_2(k,l,1)=rouge;
            condi_limites_2(k,l,2)=vert;
            condi_limites_2(k,l,3)=bleu;
            
        end
    end
    
    
    condi_limites_2=uint8(condi_limites_2);
    
    miroir=fliplr(condi_limites_2(1:hauteur,1:largeur-1,:));
    miroir2=fliplr(miroir);
    imwrite([miroir2,miroir],['Output_kp_ko_',num2str(kp_k0),'_phi_',num2str(taux_remplissage),'.png']);
    imwrite([miroir2,miroir],['Topology/sortie_kp_ko_',num2str(kp_k0),'_phi_',num2str(taux_remplissage),'_',num2str(m,'%06.f'),'.png']);
    figure(1)
    colormap jet
    subplot(2,4,1:2);
    
    if m>1
        res(m-1)=(t_max_sortie(m-1)-t_max_sortie(m))/(t_max_sortie(1)-t_max_sortie(2));
        semilogy(abs(res(1:m)),'.r');
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
    disp('-------------------------------------------------------')
    
    saveas(gcf,['Figure_kp_ko_',num2str(kp_k0),'_phi_',num2str(taux_remplissage),'.png']);
    saveas(gcf,['Figure/figure_kp_ko_',num2str(kp_k0),'_phi_',num2str(taux_remplissage),'_',num2str(m,'%06.f'),'.png']);
    toc
    
    % figure(2)
    % plot(index_de_m,'bd');
end
t_max_final=max(max(temp));
disp('End of convergence !');
