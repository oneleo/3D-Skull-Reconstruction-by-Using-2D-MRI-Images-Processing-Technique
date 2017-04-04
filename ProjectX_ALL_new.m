clear;
close;

[fname,fpath]=uigetfile('*.*','select file');
all_file_data=dir(fpath);
[all_file_num_add2 x]=size(all_file_data);
im_name=cell(all_file_num_add2-2,1);
for i=3:all_file_num_add2
    dot_pos=find(all_file_data(i).name=='.');
    file_num=str2num(all_file_data(i).name((dot_pos(4)+1):(dot_pos(5)-1)));
    im_name{file_num,1}=[fpath all_file_data(i).name];
    im_file(:,:,file_num)=dicomread(im_name{file_num,1});
end

%------------------------------
%���o�Ҧ��v���̤j�ƭȡA�]�w����x�}�j�p�קK���O
%------------------------------

im_max_value=max(max(max(im_file)));

%------------------------------
%���o�Ҧ��v����histogram�]�̪�ɶ��G�ܤ�90��^
%------------------------------

for j=0:im_max_value
   counts(j+1,1)=sum(sum(sum(im_file==j)));
end

%------------------------------
%��otsu�ȶiidx��
%------------------------------

%get otsu num
    p=counts/sum(counts);
    omega=cumsum(p);
    %im_max_value+1 type is uint16 not p's type double, so turn to double type.
    mu=cumsum(p.*double((1:im_max_value+1)'));
    mu_t=mu(end);
    sigma_b_squared=(mu_t*omega-mu).^2./(omega.*(1-omega));
    maxval=max(sigma_b_squared);
    idx=mean(find(sigma_b_squared == maxval));
%get end

%------------------------------
%�ϥ�otsu�ȨӨ��ohistogram��Ӥs����m�P�s����m�]�]otsu�j��s���p��ĤG�s���^
%------------------------------

%�s���G�Ĥ@smaller_position=9�P�ĤGgreater_position=161
%�s���Gmedium_position=86

counts_smaller_idx(1:idx)=counts(1:idx);
counts_greater_idx(1:im_max_value-idx)=counts(idx+1:im_max_value);

[smaller_value smaller_position]=max(counts_smaller_idx);
[greater_value greater_position_sub_idx]=max(counts_greater_idx);
greater_position=greater_position_sub_idx+idx;

smaller_to_idx(smaller_position-smaller_position+1:idx-smaller_position+1)=counts(smaller_position:idx);
[medium_value medium_position_sub_smaller_add_1]=min(smaller_to_idx);
medium_position=medium_position_sub_smaller_add_1+smaller_position-1;


%------------------------------
%���]idx���Ĥ@�s���P�i�������]�]�I�����Ӧh�n�L�դG�Ȥƪ��ȡ^
%------------------------------

idx_new=(smaller_position+medium_position)/2;

%------------------------------
%�v���B�z
%------------------------------

%���ưѼ�
%fspecial('average',xxx)���Pones(xxx,xxx)/sum(sum(ones(xxx,xxx)))
se_average=fspecial('average',3);

%�U�Q�Ѽ�
%se_unsharp = fspecial('unsharp',1);%1 is max unsharp

%open�Ѽ�
se_open = strel('disk',2);

%close�Ѽ�
se_close = strel('disk',5);

for k=1:all_file_num_add2-2

    %����
    blurred_im_file(:,:,k) = imfilter(im_file(:,:,k),se_average,'replicate');

    %�G��
    threshold_im_file(:,:,k)=(blurred_im_file(:,:,k)>idx_new);

    %open
    open_im_file(:,:,k)=imopen(threshold_im_file(:,:,k),se_open);
    
    %close
    close_im_file(:,:,k)=imclose(open_im_file(:,:,k),se_close);

    %�̪����
    boundaries_B=boundaries(close_im_file(:,:,k));
    boundaries_B_cell_1=boundaries_B{1};
    [close_im_file_row close_im_file_column]=size(close_im_file(:,:,k));
    external_edge_im_file(:,:,k)=bound2im(boundaries_B_cell_1,close_im_file_row,close_im_file_column,min(boundaries_B_cell_1(:,1)),min(boundaries_B_cell_1(:,2)));

    %��
    fill_im_file(:,:,k)=imfill(external_edge_im_file(:,:,k),'holes');

end

%�M�ܭ��
final_im_file=uint16(fill_im_file).*im_file;

%�N�v����3D���ƳB�z

final_im_file=smooth3(final_im_file);

%���ڭ̭n�ݪ���1��150�i

final_im_file_150=final_im_file(:,:,1:150);

%------------------------------
%�x�s�ܼ�final_im_file�î�����L�L���ܼƥH�קK�O���餣��
%------------------------------

save('project20090622_final_im_file_150.mat','final_im_file_150');
clear;
load project20090622_final_im_file_150;

%------------------------------
%3D�v������
%------------------------------

hiso = patch(isosurface(final_im_file_150,5),'FaceColor',[1,.75,.65],'EdgeColor','none');
lightangle(45,30); 
set(gcf,'Renderer','zbuffer');
lighting phong;
isonormals(final_im_file_150,hiso);
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);

view(45,30);

%�{������