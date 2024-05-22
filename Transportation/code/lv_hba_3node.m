

%% LV-HBA algorithm applied to transportation problem
%Fernando Real Rojas
%Universidad Rey Juan Carlos
%Contact: f.real.2018@alumnos.urjc.es


% Clear the workspace
clear all; close all; clc;

rng(1);
max_iters = 8e3;
obj_vecs = zeros(10,max_iters);
comp_time = zeros(10,1);

for i=1:10  
    tic;
    [a,f,fij,a_vec,f_vec,ovec_f,ovec_theta] = lv_hba();
    obj_vecs(i,:) = ovec_f;
    comp_time(i) = toc;
end
% plot(a_vec); hold on; plot(f_vec);


%% Computed gradients - Transportation problem

% F(a,f) = -m*f + c*a
% f(a,f,fij) = f*(log(f)-1) + (1-f)*(log(1-f)-1) + sum_{od}(d*fij^od) ...
%  + u*(1-f)
% g(a,fij) = a - fij <= 0 \forall (i,j) --> lam_{ij}
% fij - fji = {f,0,-f} if {i=o,i!=o!=d,i=d} \forall (i,o,d) --> (mu_i^{od})

% \nabla_a f
function val = dfa(c) 
    val = c;
end

% \nabla_f f
function val = dff(m)
    val = -m;
end

% \nabla_a g
function val = dga()
    val = 0;
end

% \nabla_f g toy example
function val = dgf(f,u)
    val = log(f./(1-f)) - u;
end

function val = dgfij(d,n)
    val = zeros(n,n,n,n);
    for oo=1:n
        for dd=1:n
            val(:,:,oo,dd) = d;
        end
    end
end


function val = glam(a,fij,n)
    val = zeros(n);
    for i=1:n
        for j=1:n
            val(i,j) = squeeze(sum(sum(fij(i,j,:,:),4),3)) - a(i,j);
        end
    end
end

function val = gmu(n,fij,f)
    val = zeros(n,n,n);

    for oo=1:n
        for dd=1:n
            if (oo ~= dd)
                for i=1:n
                    val(i,oo,dd) = sum(fij(i,:,oo,dd)) - sum(fij(:,i,oo,dd));
                    if i==oo
                        val(i,oo,dd) = val(i,oo,dd) - f(oo,dd);
                    end
                    if i==dd
                        val(i,oo,dd) = val(i,oo,dd) + f(oo,dd);
                    end
                end
            end
        end
    end

end

function a = box_a(a)
    a = max(a,0.01);
    a = min(a,10);
end

function f = box_f(f)
    f = max(f,0.1);
    f = min(f,0.9);
end

function z = box_z(z,z_max)
    z = min(z,z_max);
    z = max(z,0);
end


%% LV-HBA algorithm

function [a,f,fij,a_vec,f_vec,ovec_f,ovec_theta] = lv_hba()

    max_iters = 8e3;
    n = 3;

    d = 2.*(ones(n) - eye(n));

    d = [0,1,10;1,0,2;10,2,0];

    u = 3.*(ones(n) - eye(n));

    c = ones(n) + 10.*eye(n);
    c = [0,1,10;1,0,3;10,3,0];


    m = 2.*(ones(n) - eye(n));
    m = [0,2,6;2,0,1;6,1,0];


    a = 0.7.*rand(n,n) + 0.3;
    for i=1:n
        a(i,i) = 0;
    end

    f = 0.5.*(ones(n) - eye(n));
    fij = zeros(n,n,n,n);


    z = -ones(n);
    zmu = -ones(n,n,n);
    theta = f;
    theta_ij = fij;
    lam = -ones(n);
    mu = -ones(n,n,n);

    z_max = 1e3;

    gamma_1 = 5;
    gamma_2 = 5;

    eta = 1e-2;
    %eta = 1e-3;
    %eta = 1e-4;
    
    %eta = 5e-2;
    %eta = 3e-2;


    alfa = 8e-3;
    alfa = 1e-2;
    %alfa = 1e-3;
    %alfa = 1e-4;
    %alfa = 1e-2;
    
    beta = 2e-2;
    %beta = 2e-3;
    %beta = 2e-4;
    %beta = 5e-2;
    %beta = 2e-2;


    k = 0;
    
    %ck = 1e2;
    q = 0;

    a_prev = -10.*ones(n); f_prev = -10.*ones(n); fij_prev = -10.*ones(n,n,n,n);
    a_vec = [a];
    f_vec = [f];

    ovec_f = [];
    ovec_theta = [];




    Aeq = [];
    vec0 = zeros(1,2*n^2 + n^4);
    for i=1:n
        vec = vec0;
        vec((i-1)*n+i) = 1; % a_ii = 0;
        Aeq_p = [Aeq;vec];
        if rank(Aeq_p) > rank(Aeq)
            %Aeq = Aeq_p;
        end


        for oo=1:n
            for j=1:n
                vec = vec0;
                vec(2*n^2  + (i-1)*n^3 + (j-1)*n^2 + (oo-1)*n + j) = 1; % f^{oo}_{ij} = 0;
                Aeq_p = [Aeq;vec];
                if rank(Aeq_p) > rank(Aeq)
                   % Aeq = Aeq_p;
                end
            end

            for dd=1:n
                if (oo ~= dd)
                    vec = vec0;
                    vec(2*n^2 + (i-1)*n^3 + (i-1)*n^2 + (oo-1)*n + dd) = 1; % f^{od}_{ii} = 0;
                    Aeq_p = [Aeq;vec];
                    if rank(Aeq_p) > rank(Aeq)
                      %  Aeq = Aeq_p;
                    end
                end
            end
        end
    end

    % for oo=1:n
    %     vec = vec0;
    %     vec(n^2 + (oo-1)*n + oo) = 1; % f^{oo} = 0;
    %     Aeq_p = [Aeq;vec];
    %     if rank(Aeq_p) > rank(Aeq)
    %         Aeq = Aeq_p;
    %     end
    % end

    for oo=1:n
        for dd=[1:(oo-1),(oo+1):n]
            vec = vec0;
            vec(n^2 + (oo-1)*n + dd) = -1;
            for j=[1:(oo-1),(oo+1):n]
            %for j=1:n
                vec(2*n^2 + (oo-1)*n^3 + (j-1)*n^2 + (oo-1)*n + dd) = 1; % sum_j f_{oj}^{od} - sum_j f_{jo}^{od} - f^{od} = 0;
                vec(2*n^2 + (j-1)*n^3 + (oo-1)*n^2 + (oo-1)*n + dd) = -1;
            end
            Aeq_p = [Aeq;vec];
            if rank(Aeq_p) > rank(Aeq)
                Aeq = Aeq_p;
            end
            

            vec = vec0;
            vec(n^2 + (oo-1)*n + dd) = 1;
            for j=[1:(dd-1),(dd+1):n]
            %for j=1:n
                vec(2*n^2 + (dd-1)*n^3 + (j-1)*n^2 + (oo-1)*n + dd) = 1; % sum_j f_{dj}^{od} - sum_j f_{jd}^{od} + f^{od} = 0;
                vec(2*n^2 + (j-1)*n^3 + (dd-1)*n^2 + (oo-1)*n + dd) = -1;
            end
            Aeq_p = [Aeq;vec];
            if rank(Aeq_p) > rank(Aeq)
                Aeq = Aeq_p;
            end

            for i=1:n
                if (i~=oo) && (i~=dd)
                    vec = vec0;
                    for j=[1:(i-1),(i+1):n]
                        vec(2*n^2 + (i-1)*n^3 + (j-1)*n^2 + (oo-1)*n + dd) = 1; % sum_j f_{ij}^{od} - sum_j f_{ji}^{od} = 0;
                        vec(2*n^2 + (j-1)*n^3 + (i-1)*n^2 + (oo-1)*n + dd) = -1;
                    end
                    Aeq_p = [Aeq;vec];
                    if rank(Aeq_p) > rank(Aeq)
                        Aeq = Aeq_p;
                    end

                end
            end
        end
    end

    A = [];
    for i=1:n
        for j=1:n
            vec = vec0;
            vec((i-1)*n + j) = -1;

            for oo=1:n
                for dd=[1:(oo-1),(oo+1):n]
                    vec(2*n^2 + (i-1)*n^3 + (j-1)*n^2 + (oo-1)*n + dd) = 1;
                end
            end
            A = [A;vec];
        end
    end

    [rows_A,~] = size(A);
    [rows_Aeq,~] = size(Aeq);
    

    % while ((abs(a-a_prev)./alfa > 1e-6) || ...
    %         (abs(f-f_prev)./alfa > 1e-6)) &&  ...
    %         (k < max_iters)
    numdatos = 4e2;
    obj_vec_f = zeros(1,numdatos);
    obj_vec_theta = zeros(1,numdatos);
    
    hplot = plot(1:numdatos,obj_vec_f,'--b',1:numdatos,obj_vec_theta,'-r','LineWidth',1.5);
    xlabel('Iterations x 100');
    ylabel('Objective');
    legend('Obj(a,f)','Obj(a,\theta)','Location','northwest');
    tk = 0;
    max_iters = 8e3;
    
    %while (tk < 3.5e2) && (  (max(max(abs(a-a_prev)))/alfa > 1e-3) ||  (max(max(abs(theta-f))) > 1e-2)  )
    while (k < max_iters)
        ck = 0.1.*(k+1).^0.4;
        %ck = 0.025.*(k+1).^0.3;
        %ck = 0.025.*(k+1).^0.7;
        % gamma_1 = 1e6.*(k^(-1.1));

        % a_vec = [a_vec a]; f_vec = [f_vec f];

        
        a_prev = a;
        f_prev = f;
        dtheta = zeros(n);
        dtheta_ij = zeros(n,n,n,n);
        %disp(f)
        for oo=1:n
            for dd=1:n
                if (oo~=dd)
                    dtheta(oo,dd) = dgf(theta(oo,dd),u(oo,dd)) + ...
                    mu(dd,oo,dd) - mu(oo,oo,dd) + (1./gamma_1).*(theta(oo,dd) - f(oo,dd));
                end
            end
        end

        for oo=1:n
            for dd=1:n
                if (oo~=dd)
                    for i=1:n
                        for j=1:n
                            if i~=j
                                dtheta_ij(i,j,oo,dd) = d(i,j) + ...
                                mu(i,oo,dd) - mu(j,oo,dd) + lam(i,j) + ...
                                (1./gamma_1).*(theta_ij(i,j,oo,dd)-fij(i,j,oo,dd));
                            end
                        end
                    end
                end
            end
        end

       % disp(dtheta)
        dlam = -glam(a,theta_ij,n) + (1./gamma_2).*(lam - z); 
        dmu = -gmu(n,theta_ij,theta) + (1./gamma_2).*(mu - zmu);
        
    
        theta = theta - eta.*dtheta;
        theta = max(theta,0.0002);
        theta = min(theta,0.98);
        theta_ij = theta_ij - eta.*dtheta_ij;
        theta_ij = max(theta_ij,0);
        theta_ij = min(theta_ij,1);
        mu = mu - eta.*dmu;
        lam = lam - eta.*dlam;

        
        
        if min(min(lam)) < -z_max
            lam(lam < -z_max) = -z_max;
        else
            if max(max(lam)) > z_max
                lam(lam > z_max) = z_max;
            end
        end

        if min(min(min(mu))) < -z_max
            mu(mu < -z_max) = -z_max;
        else
            if max(max(max(mu))) > z_max
                mu(mu > z_max) = z_max;
            end
        end

        
        lam = max(lam,0);

        da = (1./ck).*(dfa(c)) + lam;  
        %disp(da)

        df = zeros(n);

        for oo=1:n
            for dd=1:n
                if (oo~=dd)
                    df(oo,dd) = dgf(f(oo,dd),u(oo,dd)) - ...
                    (1./gamma_1).*(f(oo,dd) - theta(oo,dd)) + ...
                    (1./ck).*(-m(oo,dd));
                end
            end
        end



        %df = (1./ck).*dff(m) + dgf(f,u) - (1./gamma_1).*(f - theta) ;
        dfij = dgfij(d,n) - (1./gamma_1).*(fij - theta_ij);


        dz = -(1./gamma_2).*(lam - z);
        dzmu = -(1./gamma_2).*(mu-zmu);
    
        a = a - alfa.*da;
        f = f - alfa.*df;
        fij = fij - alfa.*dfij;
       % a = 10.*ones(n);
       % f = ones(n);

        a0 = reshape(a,1,[]);
        f0 = reshape(f,1,[]);

        % fij = zeros(n,n,n,n);
        % fij(1,2,1,2) = 1;
        % fij(1,3,1,3) = 1;
        % fij(2,1,2,1) = 1;
        % fij(2,3,2,3) = 1;
        % fij(3,1,3,1) = 1;
        % fij(3,2,3,2) = 1;

        

        fij0mat = fij;
        fij0 = reshape(fij,1,[]);
        



        x0 = [a0,f0,fij0];

        func = @(x) norm(x - x0).^2;
        options = optimoptions('fmincon','Display','off','ConstraintTolerance',1e-4,'OptimalityTolerance',1e-4);




        sol = fmincon(func,x0,A,zeros(rows_A,1),Aeq,zeros(rows_Aeq,1),0.0002*ones(1,2*n^2+n^4),[1e3*ones(1,n^2),0.98*ones(1,n^2+n^4)],[],options);
        
        a_vec = sol(1:n^2);
        a = reshape(a_vec,n,n);

        f_vec = sol(n^2+1:2*n^2);
        f = reshape(f_vec,n,n);

        fij_vec = sol(2*n^2+1:end);
        fij = reshape(fij_vec,n,n,n,n);



        % a = max(0.0002,a);
        % a = min(0.98,a);
        % f = max(0.0002,f);
        % f = min(0.98,f);
        % a(1,[4:8]) = 0;
        % a(2,[5:9]) = 0;
        % a(3,[6:8]) = 0;
        % a(4,[1,9,7]) = 0;
        % a(5,[1,2,8,9]) = 0;
        % a(6,[1,2,3,9]) = 0;
        % a(7,[1:4,8:9]) = 0;
        % a(8,[1:3,9,5,7]) = 0;
        % a(9,[4:8]) = 0;
        % 
        % 
        % 
         for i=1:n
             fij(i,i,:,:) = 0;
             fij(:,:,i,i) = 0;
             f(i,i) = 0;
             a(i,i) = 0;
         end
        % 
        % 
        % fij = min(0.98,fij);
        % fij = max(0.0002,fij);

    
        tz = z - beta.*dz;
        if tz < -z_max
            z = -z_max;
        else
            if tz > z_max
                z = z_max;
            else
                z = tz;
            end
        end

        tzmu = zmu - beta.*dzmu;
        if tzmu < -z_max
            zmu = -z_max;
        else
            if tzmu > z_max
                zmu = z_max;
            else
                zmu = tzmu;
            end
        end

        k = k + 1;

        f_ob = f;
        f_ob(f_ob < 1e-2) = 0;
        theta_ob = theta;
        theta_ob(theta_ob < 1e-2) = 0;
        a_ob = a;
        a_ob(a_ob < 1e-2) = 0;


        ob_f = sum(sum(m.*f_ob)) - sum(sum(c.*a_ob));
        ob_theta = sum(sum(m.*theta_ob)) - sum(sum(c.*a_ob));

        ovec_f = [ovec_f ob_f];
        ovec_theta = [ovec_theta ob_theta];
        
        if q == 1e2
          %  disp(da);
          disp(['k = ',num2str(k)]);


            % disp(['k = ',num2str(k),', f 1-2 = ',num2str(f(1,2)),', a 1-2= ',num2str(a(1,2)), ...
            %     ', lam 1-2 = ',num2str(lam(1,2)),', z = ',num2str(z(1,2)), ...
            %     ', theta 1-2 = ',num2str(theta(1,2)), ...                
            %     ', f 2-1 = ', num2str(f(2,1)), ...
            %     ', theta 2-1 = ', num2str(theta(2,1)), ...
            %     ', a 2-1 = ', num2str(a(2,1))]);
            % disp('fij^{1-2} = ');
            % disp(fij(:,:,1,2));
            % 
            % 
            % disp('fij^{2-1} = ');
            % disp(fij(:,:,2,1));
            % 
            % disp('f21^{od}');
            % disp(permute(fij(2,1,:,:),[3,4,1,2]));

            disp('aij');
            disp(a);

            disp('f');
            disp(f);

            disp('theta');
            disp(theta);

            tk = toc;

            disp(tk);

            % disp('fij^{32}');
            % disp(fij(:,:,3,2));
            % 
            % disp('fij^{12}');
            % disp(fij(:,:,1,2));
            % 
            % disp('fij0^{12}');
            % disp(fij0mat(:,:,1,2));


           % disp('vv');
          %  disp(-f(1,2)+fij(1,2,1,2)-fij(1,2,2,1)+fij(1,2,1,3)-fij(1,2,3,1));

            %numdatos = numdatos+1;
            res_f = sum(sum(m.*f)) - sum(sum(c.*a));
            res_theta = sum(sum(m.*theta)) - sum(sum(c.*a));
            obj_vec_f = [obj_vec_f(2:end),res_f];
            obj_vec_theta = [obj_vec_theta(2:end),res_theta];

            set(hplot(1),'YData',obj_vec_f);
            set(hplot(2),'YData',obj_vec_theta);
            drawnow;
            
        %    disp(['df =',num2str(df)]);
            q = 0;
        end
        q = q+1;
        %disp(q);
    end

end








