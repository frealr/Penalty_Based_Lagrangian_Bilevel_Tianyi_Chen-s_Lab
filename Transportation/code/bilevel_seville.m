

%% Inicializacion red 9 nodos

clear all;
close all;
clc;


n = 9;

distance = 10000 * ones(n, n); % Distances between arcs

for i = 1:n
    distance(i, i) = 0;
end

distance(1, 2) = 0.75;
%distance(1,2) = 0.9;
distance(1, 3) = 0.7;
distance(1, 9) = 0.9;

distance(2, 3) = 0.6;
distance(2, 4) = 1.1;

distance(3, 4) = 1.1;
distance(3, 5) = 0.5;
distance(3, 9) = 0.7;

distance(4,5) = 0.8;
distance(4,6) = 0.7;
distance(4,8) = 0.8;
%distance(4,8) = 1.8;

distance(5,6) = 0.5;
distance(5,7) = 0.7;

distance(6,7) = 0.5;
distance(6,8) = 0.4;

for i = 1:n
    for j = i+1:n
        distance(j, i) = distance(i, j); % Distances are symmetric
    end
end

 d = distance;

 u = [0,1.6,0.8,2,1.6,2.5,3,2.5,0.8; ...
    2,0,0.9,1.2,1.5,2.5,2.7,2.4,1.8; ...
    1.5,1.4,0,1.3,0.9,2,1.6,2.3,0.9; ...
    1.9,2,1.9,0,1.8,2,1.9,1.2,2; ...
    3,1.5,2,2,0,1.5,1.1,1.8,1.7; ...
    2.1,2.7,2.2,1,1.5,0,0.9,0.9,2.9; ...
    2.8,2.3,1.5,1.8,0.9,0.8,0,1.3,2.1;...
    2.8,2.2,2,1.1,1.5,0.8,1.9,0,0.3; ...
    1,1.5,1.1,2.7,1.9,1.8,2.4,3,0];

 c = [0,1.7,2.7,0,0,0,0,0,2.9; ...
     1.7,0,2.1,3,0,0,0,0,0; ...
     2.7,2.1,0,2.6,1.7,0,0,0,2.5; ... 
     0,3,2.6,0,2.8,2.4,0,3.2,0; ...
     0,0,1.7,2.8,0,1.9,3,0,0; ...
     0,0,0,2.4,1.9,0,2.7,2.8,0; ...
     0,0,0,0,3,2.7,0,0,0; ...
     0,0,0,3.2,0,2.8,0,0,0; ...
     2.9,0,2.5,0,0,0,0,0,0];

c(c==0) = 1e1;

m = [0,2.5,2.5,5,5,7.5,7.5,7.5,2.5;...
    2.5,0,2.5,2.5,5,5,7.5,5,5;...
    2.5,2.5,0,2.5,2.5,5,5,7.5,2.5;...
    5,2.5,2.5,0,2.5,2.5,5,2.5,7.5;...
    5,5,2.5,2.5,0,2.5,2.5,5,2.5;...
    7.5,5,5,2.5,2.5,0,2.5,2.5,7.5;...
    7.5,7.5,5,5,2.5,2.5,0,5,7.5;...
    7.5,5,5,2.5,5,2.5,5,0,7.5;...
    2.5,5,2.5,5,5,7.5,7.5,7.5,0];

demand = 1e-1.*[0,9,26,19,13,12,13,8,11;
          11,0,14,26,7,18,3,6,12;
          30,19,0,30,24,8,15,12,5;
          21,9,11,0,22,16,25,21,23;
          14,14,8,9,0,20,16,22,21;
          26,1,22,24,13,0,16,14,12;
          8,6,9,23,6,13,0,11,11;
          9,2,14,20,18,16,11,0,4;
          8,7,11,22,27,17,8,12,0];



%% Inicializacion red Sevilla
  clear all; close all; clc;
% % 
 [n,link_cost,station_cost,link_capacity_slope,...
     station_capacity_slope,demand,prices,...
     op_link_cost,congestion_coef_stations,...
     congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
     a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();

     coordinates = readtable('../sevilla_net/coordinates.xlsx');
    coor_x = table2array(coordinates(1:24,3));
    coor_y = table2array(coordinates(1:24,7));




demand = demand./2e2;


adjacency_matrix = zeros(n, n);
k = 3;

for i = 1:n
    % Ordenar las distancias y obtener los índices correspondientes
    [~, sorted_indices] = sort(travel_time(i, :));
    % Conservar solo los k-vecinos más cercanos
    k_nearest_indices = sorted_indices(2:k+1); % El primer vecino es el propio nodo, por eso se omite
    % Establecer las conexiones en la matriz de adyacencia
    adjacency_matrix(i, k_nearest_indices) = 1;
    adjacency_matrix(k_nearest_indices, i) = 1; % La matriz de adyacencia es simétrica
end

adjacency_matrix = (adjacency_matrix + adjacency_matrix')./2;
adjacency_matrix(adjacency_matrix > 0.5) = 1;
adjacency_matrix(adjacency_matrix <= 0.5) = 0;
adjacency_matrix(travel_time > 7) = 0;

travel_time(adjacency_matrix == 0) = 1e2;

d = 0.25.*travel_time;
u = 0.25.*alt_time;
link_cost(adjacency_matrix == 0) = 1e8;
c = link_cost./2e2;

for i=1:n
    for j=1:n
        if adjacency_matrix(i,j) > 0
            adjacency_matrix(i,j) = d(i,j);
        end
    end
end

coor_y(9) = 37.3820073;
coor_x(9) = -5.9941967;
coor_x(21) = -5.9959814;
coor_x(17) = coor_x(17) + 0.004;
coor_y(7) = coor_y(7) - 0.002;

scaler = 0.7;
I = imread('../sevilla_net/sevilla.png', 'BackgroundColor', [1 1 1]);
%colormap gray;
I = imresize(I, scaler);
[filas, columnas, ~] = size(I);
figure('Position', [100, 100, columnas, filas]);  
g = graph(adjacency_matrix);
h = plot(g,'XData',scaler.*(coor_x-mean(coor_x)),'YData',scaler.*(coor_y-mean(coor_y)),'MarkerSize',5,'NodeFontSize',10,'Linewidth',2,'Interpreter','latex');
hold on
h = image(xlim+scaler.*0.008,-0.8*ylim-0.001.*scaler,I); 
uistack(h,'bottom');
set(h,'AlphaData',0.8);
xticks([]); yticks([]); 

% Obtener dimensiones de la imagen
%
% Configurar los ejes para que coincidan con las dimensiones de la imagen
%axis([0.5, columnas+0.5, 0.5, filas+0.5]);
axis off; % Desactivar los ejes
saveas(gcf, 'sevilla_bilevel.png');
dis_matrix = zeros(n);
for i=1:n
    for j=[1:(i-1),(i+1):n]
        [~, dis_matrix(i,j)] = shortestpath(g, i, j);
    end
end

rng(1);
b1 = rand(n);
b1(b1 > 0.5) = 1;
b1(b1 <= 0.5) = -1;
b2 = rand(n);
b2(b2 > 0.5) = 2;
b2(b2 <= 0.5) = -2;

m = dis_matrix + b1 + b2;
m = max(m,0.25);
%m = 0.75.*m;m



%% Check computational times for different initializations for inner problem

x = [15,8,5,1];
x = [10,15,20];
x = [4];


comp_time = zeros(1,1);
obj_vecs = zeros(1,1.1e6);
gamma = 3;
rng(1);

for i=1:1
    a_start = 15 + 25.*rand(n);
   %a_start = a;
  %  a_start (d > 1e2) = 0.01;
    tic;
    beta = 1.6e-4;
    %beta = 3.2e-5;
    %beta = 3.2e-5;
  % [f,fij,lam,mu] = A2_ods(n,d,u,a,zeros(n),zeros(n,n,n,n),zeros(n),zeros(n,1));
  %  [f,fij,lam,mu] = A2_f_ods(n,d,u,a,gamma,m,zeros(n),zeros(n,n,n,n),zeros(n),zeros(n,1));

    [a,fg,fijg,lamf,muf,a_vec,f_vec,fij_vec,lamf_vec,muf_vec,obj_vec,obj_vec_f] = A3_ods(n,d,u,gamma,m,c,a_start,beta,demand);

    obj_vecs(i,1:length(obj_vec)) = obj_vec;
    %obj_vecs(i,length(obj_vec)+1:end) = max(obj_vec).*ones(1,length(obj_vecs)-length(obj_vec));
    comp_time(i) = toc;
    disp(['i = ',num2str(i),', comp_time = ',num2str(comp_time(i))]);
end

% %% Results representation
% 
% 
% close all;
% 
% map_x = [2,6,6,11,11,13,14,14,1];
% map_y = [4,6,2,6,2,4,1,7,1];
% 
% 
% figure(1);
% for i=1:length(x)
%     subplot(2,2,i);
%    % imagesc(a_vecs(:,:,i));
%    a = a_vecs(:,:,i);
%    a(a < 0.01) = 0;
%    a = a  + a';
%    g = graph(a);
%    % h = plot(g);
% 
%     h = plot(g,'XData',map_x,'YData',map_y,'LineWidth',0.7.*g.Edges.Weight.^0.7,'NodeFontSize',5,...
%     'EdgeColor','#0072BD','EdgeAlpha',0.8,'interpreter','latex');%, ...
%          %   'MarkerSize',0.7*(s_h+s).^0.5 +1e-2);%,'LineWidth', ...
%             %0.1.*g.Edges.Weight,'NodeColor',colores,'EdgeColor',colores_edg,'EdgeAlpha',0.7,'NodeFontSize',8);
%         xticks([]); yticks([]); 
% 
%    tit = sprintf('$ \\gamma = %d $, $ A_0 =  %.2f$, obj = $ %.2f$, t = $ %.2f$',x(i),1,max(obj_vecs(i,2:end)),comp_time(i));
%    title(tit,'Interpreter','latex','FontSize',5);
% 
%    % title(tit,'FontSize',8);
% end
% 
% figure(2);
% for i=1:length(x)
%     subplot(2,2,i);
%     imagesc(f_vecs(:,:,i));
%     colorbar;
%     hold on;
%     for ii = 1:n
%         for jj = 1:n
%             text(jj, ii, num2str(f_vecs(ii, jj,i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white');
%             hold on;
%         end
%     end
% end





% 
% %% Lower level for different configurations of the network
% 
% close all;
% xx = 0:0.01:1; x_x = 1-xx;
% y1 = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,2) + d(2,4) + d(4,8)).*xx;
% xx = fg_vec(1); x_x = 1-xx;
% y1_op = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,2) + d(2,4) + d(4,8)).*xx;
% 
% 
% 
% xx = 0:0.01:1; x_x = 1-xx;
% y2 = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,4) + d(4,8)).*xx;
% xx = fg_vec(2); x_x = 1-xx;
% y2_op = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,4) + d(4,8)).*xx;
% 
% 
% xx = 0:0.01:1; x_x = 1-xx;
% y3 = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,5) + d(5,6) + d(6,8)).*xx;
% xx = fg_vec(3); x_x = 1-xx;
% y3_op = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,5) + d(5,6) + d(6,8)).*xx;
% 
% 
% xx = 0:0.01:1; x_x = 1-xx;
% subplot(311); 
% plot(xx,y1); hold on; scatter(fg_vec(1),y1_op,'red','filled'); xlabel('F'); ylabel('PAX(A,F)'); title('A_0 = 0.1')
% subplot(312); plot(xx,y2); hold on; scatter(fg_vec(2),y2_op,'red','filled'); xlabel('F'); ylabel('PAX(A,F)'); title('A_0 = 0.6')
% subplot(313); plot(xx,y3); hold on; scatter(fg_vec(3),y3_op,'red','filled'); xlabel('F');ylabel('PAX(A,F)'); title('A_0 = 1')



%% Functions - all ods



function [f,fij,lam,mu] = A2_ods(n,d,u,a,f_last,fij_last,lam_last,mu_last,demand)
    lam = lam_last;
    f = f_last;
    fij = fij_last;
    mu = mu_last;

    mu_o = zeros(n);
    mu_d = zeros(n);

    mu_i = zeros(n,n,n,n);
    mu_j = mu_i;

    lam_prev = -10*ones(n); %i,j
    f_prev = -10*ones(n); %od
    fij_prev = -10*ones(n,n,n,n); %i,j,o,d
    mu_prev = -10*ones(n,n,n); %i,o,d

    beta_1 = 1e-2;
    beta_2 = 4e-1;
    q = 0;
    max_dif = 1;

    while ( mean(mean(abs(f-f_prev)))/(beta_1) > 2e-2 ) ...
       || ( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1) > 1e-2 ) ...
       || (cons_f_a > 10) ...
       || (max_dif > 2e-1)
               % (mean(min((squeeze(sum(sum(fij,4),3)) < (a+3e-2)))) < 0.5) ... %sin demanda

        lam_prev = lam;
        f_prev = f;
        fij_prev  = fij;
        mu_prev = mu;
   
     %   disp(max_dif)

       

       mu = mu + beta_2.*(squeeze(sum( fij , 2 ))  - ...
          squeeze(permute(sum(fij , 1 ),[2 1 3 4])));

%         for o=1:n
%             for des=2
%                 mu(:,o,des) = mu(:,o,des) + beta*(sum( fij(:,:,o,des),2 ) - sum(fij(:,:,o,des) , 1)' );
%             end
%         end

        for o=1:n
            for des=1:n
                mu(o,o,des) = mu(o,o,des) + beta_2*(-f(o,des));
                mu_o(o,des) = mu(o,o,des);
                mu(des,o,des) = mu(des,o,des) + beta_2*(f(o,des));
                mu_d(o,des) = mu(des,o,des);
            end
        end

     %   mu = max(mu,0);

        for o=1:n
            for des=1:n
                mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
                mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
            end
        end

        for i=1:n
            for j=1:n %cambiar con los candidatos
                lam(i,j) = lam(i,j) + beta_2.*(   sum(sum(squeeze(permute(fij(i,j,:,:),[3,4,1,2]).*demand))) - a(i,j) );
            end
        end

%         for i=1:n
%             for j=1:n
%                 lam(i,j) = lam(i,j) + beta*(  sum(sum(squeeze(fij(i,j,:,:)))) - a(i,j));
%             end
%         end

        lam = max(0,lam);

        for i=1:n
            lam(i,i) = 100;
        end

        f = f - beta_1*(demand.*log(f./(1-f))  - demand.*u - mu_o + mu_d );

%         for o=1:n
%             for des=1:n
%                 f(o,des) = f(o,des) - beta*(log(f(o,des)./(1-f(o,des)))  - u(o,des) - mu(o,o,des) + mu(des,o,des) ); %puedo quitar algun bucle?
%             end
%         end

        
        f = min(0.99,f);
        f = max(0.01,f);

        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( demand(o,des).*d + lam + mu_i(:,:,o,des) - mu_j(:,:,o,des) );
            end
        end

%         for o=1:n
%             for des=1:n
%                 for i=1:n
%                     for j=1:n
%                         fij(i,j,o,des) = fij(i,j,o,des) - beta*(d(i,j) + lam(i,j) + mu(i,o,des) - mu(j,o,des)  );
%                     end
%                 end
%                % fij(:,:,o,des) = fij(:,:,o,des) - beta*(d + lam + mu_i - mu_j);
%             end
%         end

        fij = min(1,fij);
        fij = max(0,fij);

        for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
            mu(:,i,i) = 0;
        end

        cons_f_a = 0;
        for i=1:n
            for j=1:n %cambiar con los candidatos
                if squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) >= (a(i,j) + 3e-1)
                    cons_f_a = cons_f_a + 1;
                end
            end
        end



        max_dif = 0;
        
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd))  );
                end
                
            end
        end


        if q == 100
%             disp('fij 1,2 = ');
%             disp(fij(:,:,1,2));
%             disp('f 1,2 = ');
            %disp(f);
            q = 0;
        else
            q = q+1;
        end
        
        
    end

end

function [f,fij,lam,mu] = A2_f_ods(n,d,u,a,gamma,m,f_last,fij_last,lam_last,mu_last,demand)
    
    lam = lam_last;
    f = f_last; 
    fij = fij_last;
    mu = mu_last;

    lam_prev = -10*ones(n);
    f_prev = -10*ones(n);
    fij_prev = -10*ones(n,n,n,n);
    mu_prev = -10*ones(n,n,n);

    beta_1 = 1e-2/gamma; %gamma 10
    beta_2 = 4e-1/gamma; %gamma 10

    max_dif = 1;
    cons_f_a  = 20;
   

    while ( mean(mean(abs(f-f_prev)))/(beta_1*gamma) > 3e-2 ) ...
       || (cons_f_a  > 10) ...
       || (max_dif > 3e-1)  ...
          || ( mean(mean(max(max(abs(fij-fij_prev)))))/(beta_1*gamma) > 1e-1 )
        %|| (mean(mean((squeeze(sum(sum(fij,4),3)) < (a+1e-1)))) < 0.5) ...



        lam_prev = lam;
        f_prev =f;
        fij_prev = fij;
        mu_prev = mu;

       % disp(max_dif)
       %disp(['fij 1-2 = ',num2str(max(max(fij(:,:,1,2)))),', f 1-2 = ',num2str(f(1,2))]);
    

        mu = mu + beta_2.*(squeeze(sum( fij , 2 ))  - ...
          squeeze(permute(sum(fij , 1 ),[2 1 3 4])));

        for o=1:n
            for des=1:n
                mu(o,o,des) = mu(o,o,des) + beta_2.*(-f(o,des));
                mu_o(o,des) = mu(o,o,des);
                mu(des,o,des) = mu(des,o,des) + beta_2.*(f(o,des));
                mu_d(o,des) = mu(des,o,des);
            end
        end

      %  mu = max(mu,0);

        for o=1:n
            for des=1:n
                mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
                mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
            end
        end

        for i=1:n
            for j=1:n %cambiar a candidatos
                lam(i,j) = lam(i,j) + beta_2.*(   sum(sum(squeeze(permute(fij(i,j,:,:),[3,4,1,2]).*demand))) - a(i,j) );
            end
        end
        lam = max(lam,0);

        for i=1:n
            lam(i,i) = 100;
        end

        f = f - beta_1*(-m  + gamma.*demand.*(log(f./(1-f))  - u) - mu_o + mu_d );
        f = min(0.99,f);
        f = max(0.01,f);



        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( gamma*demand(o,des).*d + mu_i(:,:,o,des) - mu_j(:,:,o,des) + lam);
            end
        end

        fij = min(1,fij);
        fij = max(0,fij);

        cons_f_a = 0;
        for i=1:n
            for j=1:n %cambiar a candidatos
                if squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) >= (a(i,j) + 3e-1)
                    cons_f_a = cons_f_a + 1;
                end
            end
        end

        max_dif = 0;
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd))  );

                end
                
            end
        end
      %  disp(dif_ok./total_od);
       % disp(max_dif);

        for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
            mu(:,i,i) = 0;
        end

    end

end


function [f,fij] = A1_ods(n,d,u,a,f_last,fij_last,lam,mu)
    f = f_last;
    fij = fij_last; 

    mu_o = zeros(n);
    mu_d = zeros(n);

    mu_i = zeros(n,n,n,n);
    mu_j = mu_i;

    for o=1:n
        for des=1:n
            mu_o(o,des) = mu(o,o,des);
            mu_d(o,des) = mu(des,o,des);
        end
    end

    for o=1:n
        for des=1:n
            mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
            mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
        end
    end

    f_prev = -10*ones(n);
    fij_prev = -10*ones(n,n,n,n); %i,j,o,d

    beta_1 = 1e-2;

    q = 0;
    max_dif = 1;

    while ( max(max(abs(f-f_prev)))/(beta_1) > 4e-2 ) ...
   || ( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1) > 1e-2 ) ...
   || (max_dif > 7e-2)

        f_prev = f;
        fij_prev  = fij;

        % disp(max_dif);
        % disp( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1)   );

        f = f - beta_1*(log(f./(1-f))  - u - mu_o + mu_d );
        f = min(0.99,f);
        f = max(0.01,f);

        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( d + lam + mu_i(:,:,o,des) - mu_j(:,:,o,des) );
            end
        end
        fij = min(1,fij);
        fij = max(0,fij);

       for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
        end


        max_dif = 0;
        
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd)));
                end
                
            end
        end


    end
    disp('converge A1');


end


function [a_best,fg_best,fijg_best,lamf_best,muf_best,a_vec,f_vec,fij_vec,lamf_vec,muf_vec,obj_vec,obj_vec_f] = A3_ods(n,d,u,gamma,m,c,a_start,beta_base,demand)

    a = a_start;

    fg = 0.01*ones(n);
    fijg = zeros(n,n,n,n);
    lamg = zeros(n);
    mug = zeros(n,n,n);

    ff = fg;
    fijf = fijg;
    lamf = lamg;
    muf = mug;

    a_prev = -10*ones(n);
    f_prev = -10*ones(n);
    fij_prev = -10*ones(n,n,n,n);
    lam_prev = -10*ones(n);
    mu_prev = -10*ones(n,n,n);

    % beta = 1e-4; %gamma 10
    % beta = 1e-5;
    % beta = 1e-6;
    
    beta = beta_base;

    numdatos = 3e2;
    obj_vec_lasts_g = zeros(1,numdatos);
    obj_vec_lasts_f = zeros(1,numdatos);
    
    hplot = plot(1:numdatos,obj_vec_lasts_g,'-b',1:numdatos,obj_vec_lasts_f,'--b','LineWidth',1.5);
    xlabel('Last 300 results');
    ylabel('F(a,f(a))');
    ylim([-10,inf]);
    %legend('F(a,f)','Obj(a,\theta)','Location','northwest');



    a_vec = [a];
    f_vec = [fg];
    fij_vec = [fijg];
    lamf_vec = [lamf];
   % lamg_vec = [lamg];
    muf_vec = [muf];
   % mug_vec = [mug];
    obj_vec = [];
    obj_vec_f = [];
    lamg_prev = -10*ones(n);
    max_iters = 1.1e6;
    obj_best = -1e6;
    k = 0;
    q = 0;
    control = 0;

    % while ((  max(max(abs(fg - f_prev))) > 1e-2) || ( sum(sum((abs(a-a_prev)./a_prev)/beta))  > 1) || ...
    %       (max(max(max(max(abs(fijg-fij_prev))))) > 1e-2) || ...
    %       (max(max(abs(lamf-lam_prev)))/beta > 1e-2) ||...
    %       (max(max(max(abs(muf-mu_prev))))/beta > 1e-2)) && ((k < max_iters)) %  ||  (obj <= (max(obj_vec) - 0.5 )) ) && (k < (max_iters + 2e4))
    while (k < max_iters)

        %disp(['a gap = ',num2str( sum(sum((abs(a-a_prev)./a_prev)/beta)) )]);

        
     %   disp( (abs(a-a_prev)./a_prev) );

     
        a_prev = a;
        f_prev = fg;
        ff_prev = ff;
        fij_prev = fijg;
        fijf_prev = fijf;
        lam_prev = lamf;
        lamg_prev = lamg;
        mu_prev = muf;
        mug_prev = mug;

        %disp(lamg);

        if (max(max(lamg_prev - 1000*eye(n))) < 1) && (k > 1)
            beta = 1e-2;
          %  disp('acelero');
        else
            % beta = 1e-4; %gamma 10
            % beta = 1e-5;
            % beta = 1e-6;
            beta = beta_base;
     
        end
       

        [fg,fijg,lamg,mug] = A2_ods(n,d,u,a,f_prev,fij_prev,lamg_prev,mug_prev,demand);
       % disp('converge el lower');
        [ff,fijf,lamf,muf] = A2_f_ods(n,d,u,a,gamma,m,ff_prev,fijf_prev,lam_prev,mu_prev,demand);
        gt = c + gamma.*lamg - lamf;
        a = a - beta.*gt;


         %  disp(['dlamg/dx = '   ,num2str(max(max(abs(lamg-lamg_prev)./abs(a-a_prev))))]);

        a = max(0.001,a);
      %  a = min(0.999,a);
      %  a_vec = [a_vec a(1,3)];
      %  f_vec = [f_vec fg];

      aa_ob = a;
      aa_ob(aa_ob < 1e-2) = 0;
      ff_ob = fg;
      ff_ob_f = ff;
      ff_ob(ff_ob < 2e-2) = 0;
      ff_ob_f(ff_ob_f < 2e-2) = 0;

        obj = sum(sum(m.*demand.*ff_ob)) - sum(sum(aa_ob.*c));
        obj_f = sum(sum(m.*demand.*ff_ob_f)) - sum(sum(aa_ob.*c));

        obj_vec = [obj_vec obj];
        obj_vec_f = [obj_vec_f obj_f];


        if (max(obj_vec) == obj)
            a_best = a;
            fg_best = fg;
            fijg_best = fijg;
            lamf_best = lamf;
            muf_best = muf;
        end

      %  lamf_vec = [lamf_vec lamf];
      %  lamg_vec = [lamg_vec lamg];
      %  muf_vec = [muf_vec muf];
     %   mug_vec = [mug_vec mug];

        % if (beta == 1e-2)
        %     k = k + 100;
        %     q = q + 100;
        %     control = control + 100;
        % else
            k = k+1;
            q = q+1;
            control = control + 1;
        % end
    %    if q >= 20000 gamma 10
         if (q >= 20000)  %gamma 20
%             disp('fijg = ');
%             disp(fijg);
%             disp('fijf = ');
%             disp(fijf);
            disp('fg = ');
            disp(fg);
            disp('a = ');
            disp(a);
            q = 0;


            % if ( max( obj_vec ) <= (obj_best*1.005) )% && (k < max_iters))
            %      k = max_iters;
            % end
            obj_best = max(obj_vec);
            disp(['k = ',num2str(k),', best obj = ',num2str(obj_best)]);
        end

        if control >= 1e3
            disp(['k = ',num2str(k),', obj = ',num2str(obj),', beta = ',num2str(beta)]);
            ares = a;
            ares(ares < 0.1) = 0;
            fgres = fg;
            ffres = ff;
            fgres(fgres<0.02) = 0;
            ffres(ffres<0.02) = 0;
            res_g = sum(sum(m.*fgres.*demand))-sum(sum(c.*ares));
            res_f = sum(sum(m.*ffres.*demand))-sum(sum(c.*ares));
            obj_vec_lasts_g = [obj_vec_lasts_g(2:end),res_g];
            obj_vec_lasts_f = [obj_vec_lasts_f(2:end),res_f];

            set(hplot(1),'YData',obj_vec_lasts_g);
            set(hplot(2),'YData',obj_vec_lasts_f);
            drawnow;
           % disp(a);
           % disp(max(max(lamg-1000.*eye(n))));
            control = 0;
        end
        
       % disp(['k = ',num2str(k),', x = ',num2str(x),', y = ',num2str(yg), ', lam = ',num2str(lamf) , ', lam_op = ',num2str(lamg)]);
    end
   % disp(a);
   % disp(fg);
   % disp(fijg);
    %disp(['a = ',num2str(x),', y = ',num2str(yg), ', lam = ',num2str(lamf) , ', lam_op = ',num2str(lamg)]);
end





function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
    a_max,candidates,pond_coefs_mat] = parameters_sevilla_network()

    n = 24;
    
    %candidates to construct a link for each neighbor
    candidates = {
    [4,23,2,18,3,9,8,10,24,5,17];
    [1,4,23,18,3,9,8,5,24];
    [1,2,18,16,14,21,19,9,8,24];
    [1,2,23,16,15,20,11,22,6];
    [1,2,10,9,21,14,17,12,13,24];
    [4,22,16,14,19,7,20];
    [12,13,17,19,21,14,11,22,6,20];
    [1,2,3,9,10,24];
    [1,2,3,16,21,20,19,17,5,10,8,12];
    [1,8,9,20,19,12,17,5,24,21];
    [4,22,20,7,19,14,16,15,23];
    [13,5,17,10,9,19,14,7,20];
    [24,5,17,19,7,12];
    [21,3,18,16,11,22,6,20,7,12,19,17,5,15];
    [18,23,4,22,11,20,14,16];
    [18,23,4,15,11,6,14,21,9,3,22];
    [5,24,1,10,9,21,14,19,7,12,13];
    [2,23,15,16,14,21,19,3,1];
    [13,17,10,9,3,18,21,14,11,6,20,7,12];
    [6,22,11,4,15,14,21,9,10,19,7,12];
    [3,18,23,16,14,20,7,19,17,5,10,9];
    [4,6,20,7,14,16,11,15,23];
    [1,4,22,11,15,16,21,18,24,2];
    [1,2,23,3,8,10,17,5,13];
    };
    

    population_file = readtable('../sevilla_net/population.xlsx');
    population = table2array(population_file(1:24,2));
    coordinates = readtable('../sevilla_net/coordinates.xlsx');
    coor_x = table2array(coordinates(1:24,3));
    coor_y = table2array(coordinates(1:24,7));
    rng(1,"twister"); %seed
    distance = 1e6.*ones(n);
    for i=1:n
        distance(i,i) = 0;
        cand = candidates(i);
        cand = cand{1};
        cand = cand(cand > i);
        for j=i+1:n
            if sum(j == cand) > 0
                distance(i,j) = haversine(coor_y(i), coor_x(i), coor_y(j), coor_x(j));
                distance(j,i) = distance(i,j);
            end
            non_stop = rand < 0.4;

            alt_cost(i,j) = haversine(coor_y(i), coor_x(i), coor_y(j), coor_x(j));
            alt_cost(i,j) = alt_cost(i,j) + 0.2*non_stop*alt_cost(i,j);
            alt_cost(j,i) = alt_cost(i,j);
        end
    end
    demand = [0, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           272, 0, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           327, 327, 0, 327, 327, 664, 327, 327, 327, 664, 664, 327, 664, 327, 664, 664, 664, 327, 1125, 1125, 1125, 1125, 1125, 1125;
           185, 185, 185, 0, 185, 376, 185, 185, 185, 376, 376, 185, 376, 185, 376, 376, 376, 185, 637, 637, 637, 637, 637, 637;
           272, 272, 272, 272, 0, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           225, 225, 225, 225, 225, 0, 225, 225, 225, 188, 188, 225, 188, 225, 188, 188, 188, 225, 284, 284, 284, 284, 284, 284;
           283, 283, 283, 283, 283, 575, 0, 283, 283, 575, 575, 283, 575, 283, 575, 575, 575, 283, 975, 975, 975, 975, 975, 975;
           272, 272, 272, 272, 272, 553, 272, 0, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           272, 272, 272, 272, 272, 553, 272, 272, 0, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           511, 511, 511, 511, 511, 428, 511, 511, 511, 0, 428, 511, 428, 511, 428, 428, 428, 511, 645, 645, 645, 645, 645, 645;
           225, 225, 225, 225, 225, 188, 225, 225, 225, 188, 0, 225, 188, 225, 188, 188, 188, 225, 284, 284, 284, 284, 284, 284;
           272, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 0, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           306, 306, 306, 306, 306, 257, 306, 306, 306, 257, 306, 257, 0, 257, 306, 306, 306, 257, 387, 387, 387, 387, 387, 387;
           294, 294, 294, 294, 294, 597, 294, 294, 294, 597, 597, 294, 597, 0, 597, 597, 597, 297, 1012, 1012, 1012, 1012, 1012, 1012;
           409, 409, 409, 409, 409, 342, 409, 409, 409, 342, 342, 409, 342, 409, 0, 342, 342, 409, 516, 516, 516, 516, 516, 516;
           511, 511, 511, 511, 511, 428, 511, 511, 511, 428, 428, 511, 428, 511, 428, 0, 428, 511, 645, 645, 645, 645, 645, 645;
           429, 429, 429, 429, 429, 360, 429, 429, 429, 306, 360, 429, 360, 429, 360, 360, 0, 429, 542, 542, 542, 542, 542, 542;
           272, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 0, 937, 937, 937, 937, 937, 937;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 660, 660;
           879, 879, 879, 879, 879, 952, 879, 879, 879, 952, 952, 879, 952, 879, 952, 952, 952, 879, 860, 0, 860, 860, 860, 860;
           715, 715, 715, 715, 715, 775, 715, 715, 715, 775, 775, 715, 775, 715, 775, 775, 775, 715, 700, 700, 700, 700, 700, 700;
           511, 511, 511, 511, 511, 553, 511, 511, 511, 553, 553, 511, 553, 511, 553, 553, 553, 511, 500, 500, 500, 0, 500, 500;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 0, 660;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 660, 0];

    pond_coefs = [1.54535294, 1.54535294, 1.11218247, 1.72777031, 2.75490793, ...
        2.14294615, 1.41994992, 1, 1.11218247, 1.41994992,...
           2.14294615, 2.35701223, 2.75490793, 2.21242079, 3.        ,...
           1.59083705, 1.41994992, 1.01231062, 1.41994992, 2.35701223,...
           1.11218247, 2.14294615, 1.72777031, 1.45204153]';
    
    pond_coefs_tens(1,:,:) = pond_coefs.*ones(1,n);
    pond_coefs_tens(2,:,:) = permute(pond_coefs_tens(1,:,:),[1 3 2]);
    pond_coefs_mat = squeeze(permute(max(pond_coefs_tens(1,:,:),pond_coefs_tens(2,:,:)),[2,3,1]));

    crec_coefs = [1.6, 1.6, 1.1469802107427398, 0.9, 1.2081587290019313, 0.9661783579052806, 1.0802156586966714, ...
        0.9, 1.1469802107427398, 1.0802156586966714, 0.9661783579052806, 0.9989307690507944,...
        1.2081587290019313, 0.9, 1.0730507219686063, 0.9, 1.0802156586966714, 1.3414585127763425, ...
        1.0802156586966714, 0.9989307690507944, 1.1469802107427398, 0.9661783579052806, 0.9, 1.1224800389355853];

    % for o=1:n
    %     for d=1:n
    %         demand(o,d) = crec_coefs(o)*crec_coefs(d)*demand(o,d);
    %     end
    % end

   % demand = demand.*pond_coefs_mat;
    %demand = max(demand,demand');
    
    %fixed cost for constructing links
    link_cost = 1e6.*distance./(365.25*25);
    
    %fixed cost for constructing stations
    station_cost = 1e3.*population./(365.25*25);
    
    link_capacity_slope = 0.3.*link_cost; 
    station_capacity_slope = 0.2.*station_cost;
    
    
    % Op Link Cost
    op_link_cost = 4.*distance;
    
    % Congestion Coefficients
    congestion_coef_stations = 0.1 .* ones(1, n);
    congestion_coef_links = 0.1 .* ones(n);
    
    % Prices
    prices = 0.1.*(distance).^(0.7);
    %prices = zeros(n);
    
    % Travel Time
    travel_time = 60 .* distance ./ 30; % Time in minutes
    
    % Alt Time
    alt_time = 60 .* alt_cost ./ 30; % Time in minutes
    alt_price = 0.1.*(alt_cost).^(0.7); %price
    
    
    a_nom = 588;             
    
    tau = 0.57;
    sigma = 0.25;
    a_max = 1e9;
    eps = 1e-3;
end

function distancia = haversine(lat1, lon1, lat2, lon2)
    % Convierte las coordenadas de grados a radianes
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Diferencias en coordenadas
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;

    % Fórmula haversine
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));

    % Radio de la Tierra en kilómetros (aproximado)
    radio_tierra = 6371;

    % Calcula la distancia
    distancia = radio_tierra * c;
end
