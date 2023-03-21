clear;clc;close all;

% Jacobian - posteriors from DCM average
A = [1.23, -21.77, -4.29; -5.70, 0.45, -7.04; -6.52, 54.12, 4.38];

% coarse grain centers
range = -0.5:0.08:0.5;

% cycle through 3D grid
c_ii = 0;
Pz = 0.1232;
dt = 0.01;

count=0;
for ii = 1:numel(range)
    c_ii = c_ii + 1;
    x = range(ii);

    c_jj = 0;
    for jj = 1:numel(range)
        c_jj = c_jj + 1;
        y = range(jj);

        c_kk = 0;
        for kk = 1:numel(range)
            c_kk = c_kk + 1;
            z = range(kk);

            % the vector at time t
            v_t = [x; y; z];
            if sum(v_t.^2).^0.5<0.5
                count=count+1;
                % store the vector at time t
                %v_t_store(c_ii,c_jj,c_kk,:) = v_t;
                v_t_store(count,:) = v_t;
    
                % use first-order model to step ahead one point in time and
                % store the vector at time t + 1
                %v_tplus1(c_ii,c_jj,c_kk,:) = (A + 1)*v_t;
                %v_tplus1(count,:) = (A + 1)*v_t;
                %v_tplus1(count,:) = (A*Pz*dt + 1)*v_t;
                
                v_tplus1(count,:) = (A*Pz*dt+1)*v_t;
                v_delt(count,:) = A*Pz*v_t*dt;
            end
            %display(v_tplus1)

        end
    end
end

% normalize v_tplus1 to correlation space
%v_tplus1 = normalize(v_tplus1,'range',[-1 1]);
% 
% for ii = 1:3
%     v_tplus1(:,ii) = v_tplus1(:,ii)%normalize(v_tplus1(:,ii),'range',[-1 1]);
%     
% 
% end


%writematrix(v_t_store,'~/Desktop/vt.txt')8
%writematrix(v_tplus1,'~/Desktop/vt_plus1.txt')
%writematrix(v_t_store,'~/Desktop/vt.txt')
%writematrix(v_tplus1_alt,'~/Desktop/vt_plus1_alt.txt')
x=v_t_store(:,1);
y=v_t_store(:,2);
z=v_t_store(:,3);
s=1;
u=v_delt(:,1);
v=v_delt(:,2);
w=v_delt(:,2);

quiver3(x,y,z,u,v,w,1.3,'k','LineWidth',1);
alpha = 0.1;
view(124,-37)
xlabel('x')
ylabel('y')
zlabel('z')
% set(gcf,'Renderer','Painter')
% hgexport(gcf,'~/Desktop/k.eps');
% close all