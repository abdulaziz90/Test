
% clear
% close all
% clc
%%
addpath data/

disp '/*----------------------------------------------------------------*/'
disp ' example: Simulation of wide band radio interferometric data'
disp ' Physical model: power low spectra'
disp '/*----------------------------------------------------------------*/'

% First comment

load data64_new.mat

%/Input model image
im=fitsread('M31.fits');
im=((imresize(im,[256,256])));
im=(im+abs(im))./2;% Enforce positivity
im=im./max(im(:));
[Nix, Niy]=size(im);%Image dimensions

% rng(0);

[X,Y] = meshgrid(-Nix/2+1:Niy/2, -Nix/2+1:Niy/2);

Z = exp((-X.^2- Y.^2)./4);

Z1 = exp((-X.^2- Y.^2)./8);

alphamin=-0.01;
alphamax=0.01;
alpha_0= (alphamin + (alphamax-alphamin).*rand(Nix,Niy));

alpha_corr=conv2(alpha_0,Z,'same');
alpha_corr=alpha_corr./max(alpha_corr(:));

alpha_im=conv2(im,Z1,'same');%+alpha_correlated;
alpha_im=alpha_im./max(alpha_im(:));

alpha1=2*alpha_im+0.5*alpha_corr;

% figure(2),
% imagesc(alpha1), colorbar, axis image, colormap 'jet',

% pos = [(Nix/2)-50 (Niy/2)-50 100 100];
% rectangle('Position',pos,'LineWidth',2);

betamin=-0.01;
betamax=0.01;
beta_0= (betamin + (alphamax-betamax).*rand(Nix,Niy));

beta_corr=conv2(beta_0,Z,'same');
beta_corr=beta_corr./max(beta_corr(:));

beta_im=conv2(im,Z1,'same');%+alpha_correlated;
beta_im=beta_im./max(beta_im(:));

beta1=2*beta_im+0.5*beta_corr;


%%
xf_org = xf;
XF_org = XF;
[~,s_org,~]=svd(XF, 'econ');
s_org = diag(s_org)
r_org = rank(XF)

%%
% sim_mask_ar;
load indmap_final.mat;

%%
% band = zeros(count,2);
% band(j,:) = [band_low band_high];
% k = 64;
% for j = 1 : 8 :count
%     band(j:j+7,:) = repmat([k k-4],8,1);
%     k = k-4;
% end
% % for j = 1:count
% % band_low = randi(length(f),1);
% % band_high = band_low + randi(5,1);
% % if band_high > length(f)
% %     band_high = length(f);
% % end
% % end
%%
band = [];
band_low = 1;
for j = 1 : 50
    band(j,:) = [band_low band_low+8];
    band_low = band_low + 1;
end
%%
record_val = [];
record_band = [];
% index = img(:);
% index = d1(:);
index = indmap(:);
% val = [0 0.25 0.5 0.75 0.75];
val = 0;
% band = [5 10;
%         45 50;
%         40 45;
%         50 55;
%         10 15];
% val = repmat(val,1,round(count/numel(val))+1);
js = [];
for j = 1 : 5 : max(indmap(:))
    
    js = [js j];
    %     band_low = randi(length(f),1);
    %     band_high = band_low + 5; %randi(5,1);
    %     if band_high > length(f)
    %         band_high = length(f);
    %     end
    
    %     rand_val = val(randi(length(val),1));
    %     rand_val = j/10;
    rand_val = val;
    b = band(randi(size(band,1),1),:);
    record_band(j,:) = [b(1) b(2)];
%     record_val(j) = rand_val;
    
        XF(index == j,b(1):b(2)) = rand_val .* XF(index == j,b(1):b(2));
%     XF(index == j,band_low:band_high) = rand_val .* XF(index == j,band_low:band_high);
end
record_val = val
record_band 

for i = 1 : size(XF,2)
    xf(:,:,i) = reshape(XF(:,i),Nix,Niy);
end

[~,s_new,~]=svd(XF, 'econ');
s_new = diag(s_new)
r_new = rank(XF)

%%
% for k = 1 : max(d1(:))
%     A = XF(index == k,:);
%     for i = 1 :size(A,1)
%     figure(3),plot(A(i,:));
%     hold on;
%     %         waitforbuttonpress;
%     end
% end

%%
% load data_64r64final100s777.mat

ind1 = zeros(size(indmap));
for i = 1 : length(js)
    m = find(indmap == js(i));
    ind1(m) = i;
end
figure(1)
imagesc(ind1), colorbar, axis image, colormap 'jet',

%%
center_x = Nix/2;
center_y = Niy/2;
offset = 50;

rec = zeros(Nix,Niy);
rec (center_y-offset:center_y+offset,center_x-offset:center_x+offset) = 1;

% B = XF;
% B(~logical(rec)) = 0;
% B(all(B==0,2),:)=[];
% for i = 1 : size(B,1)
%     figure(10);
%     plot(B(i,:)),hold on;
% end


rec = logical(rec(:));
a = ind1;
% a(~rec) = 0;
% figure(2),
% imshow(label2rgb(a));

%%
load data_64r64final100s5.mat

XF = xsol2;
a = indmap;
c1 = 0;
for i = 1 : max(a(:))
    b = find(a == i);
    if ~isempty(b)
        c = XF(a == i,:);
        for j = 1 : size(c,1)
            figure(3),plot(c(j,:));hold on
            c1 = c1 + 1;
        end
    end
end
%%
% figure,imagesc(log10(XF+1e-4)),colormap jet
%%
% figure(5),
% subplot 121,imagesc(log10(xf(:,:,1)+1e-4)), axis square, colorbar, colormap 'jet',
% subplot 122,imagesc(log10(xf(:,:,end)+1e-4)), axis square, colorbar, colormap 'jet',

%%
% center_x = Nix/2;
% center_y = Niy/2;
% 
% figure(5),
% imagesc(indmap), axis square, colorbar, colormap 'jet',hold on
% box off;set(gcf, 'Position', [0 400 800 800]);
% % rectangle('Position',[center_x-64 center_x-64 128 128],'LineWidth',2,'EdgeColor','r');
% rectangle('Position',[center_x-10 center_x-10 20 20],'LineWidth',2,'EdgeColor','r');
% 
% axis([0 256 0 256]);
%%
N = Nix * Niy;
Nx = Nix;
Ny = Niy;
c = length(f);

% save data_64r64final100v777.mat xf_org XF_org s_org r_org xf XF s_new r_new Nx Ny N c f0 f indmap index alpha1 beta1 c1 record_val record_band js


