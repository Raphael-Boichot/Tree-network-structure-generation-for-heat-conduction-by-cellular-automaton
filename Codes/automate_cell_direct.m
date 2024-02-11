%https://doi.org/10.1016/j.enconman.2008.09.003
%https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton
function t_max_final = automate_cell_direct(temp_grad_ratio,kp_k0,filling_ratio,heat_sink_temperature,delta_x,p_vol,variation_rate,start_image,verbose)
rng('shuffle', 'twister')
%*****************Automate cellulaire*********************INPG/BOICHOT/2008

%get the image pixels as boudary conditions
initial_boundary_conditions=imread(start_image);
disp(['Launching case with temperature to gradient attraction = ',num2str(temp_grad_ratio)]);
[height,width,~]=size(initial_boundary_conditions);
number_of_epoch = max([height,width]);
boundary_conditions = zeros(height,width);

%first pass: count white pixels (non conductive)
non_conductive_cells=0;
for k = 1:1:height
    for l = 1:1:width
        vec = initial_boundary_conditions(k,l,:);
        if vec==[255 255 255]
            non_conductive_cells=non_conductive_cells+1;
        end
    end
end

%second pass: fill the image with conductive pixels
conductive_cells=ceil(non_conductive_cells*filling_ratio);
k=0;
while k<conductive_cells
    h=ceil(rand*height);
    l=ceil(rand*width);
    if initial_boundary_conditions(h,l,:)==[255 255 255]
        initial_boundary_conditions(h,l,:)=[0 0 0];
        k=k+1;
    end
end

%save the intial configuration
if verbose==1
    figure('Position',[100 100 800 800]);
    mirror=fliplr(initial_boundary_conditions(1:height,1:width-1,:));
    mirror2=fliplr(mirror);
    imwrite([mirror2,mirror],['Topology/sortie_kp_ko_',num2str(kp_k0),'_phi_',num2str(filling_ratio),'_',num2str(0,'%06.f'),'.png']);
end

%Third pass: congreen image into boundary conditions
non_conductive_cells=0;
black_pixels=0;
for k = 1:1:height
    for l = 1:1:width
        red = initial_boundary_conditions(k,l,1);
        green = initial_boundary_conditions(k,l,2);
        blue = initial_boundary_conditions(k,l,3);
        if (red == 255) && (green == 255) && (blue == 255)
            choice = 1;%by default k0=1, change value here if necessary
            non_conductive_cells=non_conductive_cells+1;
        end
        if (red == 127) && (green == 127) && (blue == 127)
            choice = -2;
        end
        if (red == 0) && (green == 0) && (blue == 255)
            choice = -3;
        end
        if (red == 0) && (green == 0) && (blue == 0)
            choice = kp_k0;%by default only the ratio is provided
            black_pixels=black_pixels+1;
        end
        boundary_conditions(k, l) = choice;
    end
end

temp=ones(height,width).*heat_sink_temperature;
boundary_conditions_to_pixels=initial_boundary_conditions;
t_max_sortie=zeros(number_of_epoch);
res=zeros(number_of_epoch);

for m=1:1:number_of_epoch
    tic
    if verbose==1
        disp('-------------------------------------------------------')
        disp('Applying the Cellular Automaton...');
    end
    [~,~,~,~,~,~,t_max_sortie(m),temp,grad,~]=finite_temp_direct_sparse(kp_k0,1,heat_sink_temperature,delta_x,p_vol,boundary_conditions);
    gradients=zeros(height,width);
    gradients(2:height-1,2:width-1)=grad;
    grad_max = max(max(gradients));
    grad_min = min(min(gradients));
    temp_max = max(max(temp));
    temp_min = min(min(temp));
    gradients=(gradients-grad_min)/(grad_max-grad_min);
    temp=(temp-temp_min)/(temp_max-temp_min);
    note=gradients*(1-temp_grad_ratio)+temp*(temp_grad_ratio);
    [boundary_conditions] = cellular_automaton(boundary_conditions,note, kp_k0,1,conductive_cells,variation_rate,number_of_epoch,m);
    % Variables output in this order :
    % 1. Distance of the hotest cell to the heat sink (scalar)
    % 2. Sum of cell entropy (scalar)
    % 3. Entropy map (matrix)
    % 4. Variance of temperatures accross the 1D adabatic borders (scalar)
    % 5. Variance of temperatures accross the 2D domain (scalar)
    % 6. Mean temperature (scalar)
    % 7. Maximal temperature accross the 2D domain (scalar)
    % 8. Map of temperatures (matrix)
    % 9. map of thermal gradients (matrix)
    % 10. Variance of gradients across the 2D domain (scalar)
    [~,~,entropy_map,~,variance,~,~,temp,~,~]=finite_temp_direct_sparse(kp_k0,1,heat_sink_temperature,delta_x,p_vol,boundary_conditions);
    if verbose==1
        for k = 1:1:height
            for l = 1:1:width
                choice = boundary_conditions(k, l) ;
                if choice == 1
                    red = 255;
                    green = 255;
                    blue = 255;
                end
                if choice == -2
                    red = 127;
                    green = 127 ;
                    blue = 127 ;
                end
                if choice == -3
                    red = 0;
                    green = 0 ;
                    blue = 255;
                end
                if choice == kp_k0
                    red = 0 ;
                    green = 0 ;
                    blue = 0;
                end
                boundary_conditions_to_pixels(k,l,1)=red;
                boundary_conditions_to_pixels(k,l,2)=green;
                boundary_conditions_to_pixels(k,l,3)=blue;
            end
        end
        
        boundary_conditions_to_pixels=uint8(boundary_conditions_to_pixels);
        mirror=fliplr(boundary_conditions_to_pixels(1:height,1:width-1,:));
        mirror2=fliplr(mirror);
        imwrite([mirror2,mirror],['Output_kp_ko_',num2str(kp_k0),'_phi_',num2str(filling_ratio),'.png']);
        imwrite([mirror2,mirror],['Topology/sortie_kp_ko_',num2str(kp_k0),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
        figure(1)
        colormap jet
        subplot(2,4,1:2);
        
        if m>1
            res(m-1)=(t_max_sortie(m-1)-t_max_sortie(m))/(t_max_sortie(1)-t_max_sortie(2));
            semilogy(abs(res(1:m)),'.r');
            title('Residuals');
        end
        
        subplot(2,4,3:4);
        imagesc([mirror2,mirror]);
        title('Topology');
        
        subplot(2,4,5);
        plot(variance,'.m');
        imagesc(log10(entropy_map(2:end-1,2:end-1)));
        title('Log10 Entropy');
        
        subplot(2,4,6);
        imagesc(temp);
        title('Temperature');
        
        subplot(2,4,7);
        imagesc(log10(gradients));
        title('Log10 thermal gradients');
        
        subplot(2,4,8);
        imagesc(log10(note));
        title('Log10 attraction function');
        
        disp(['Maximum temperature: ',num2str(t_max_sortie(m))])
        disp(['End of epoch: ',num2str(m)])
        disp('-------------------------------------------------------')
        
        saveas(gcf,['Figure_kp_ko_',num2str(kp_k0),'_phi_',num2str(filling_ratio),'.png']);
        saveas(gcf,['Figure/figure_kp_ko_',num2str(kp_k0),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
        toc
    end
end
t_max_final=max(max(temp));
disp(['Convergence of case with temperature to gradient attraction = ',num2str(temp_grad_ratio)]);
