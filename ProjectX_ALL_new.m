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
%取得所有影像最大數值，設定之後矩陣大小避免浪費
%------------------------------

im_max_value=max(max(max(im_file)));

%------------------------------
%取得所有影像的histogram（最花時間：至少90秒）
%------------------------------

for j=0:im_max_value
   counts(j+1,1)=sum(sum(sum(im_file==j)));
end

%------------------------------
%取otsu值進idx內
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
%使用otsu值來取得histogram兩個山頂位置與山谷位置（因otsu大於山谷小於第二山頂）
%------------------------------

%山頂：第一smaller_position=9與第二greater_position=161
%山谷：medium_position=86

counts_smaller_idx(1:idx)=counts(1:idx);
counts_greater_idx(1:im_max_value-idx)=counts(idx+1:im_max_value);

[smaller_value smaller_position]=max(counts_smaller_idx);
[greater_value greater_position_sub_idx]=max(counts_greater_idx);
greater_position=greater_position_sub_idx+idx;

smaller_to_idx(smaller_position-smaller_position+1:idx-smaller_position+1)=counts(smaller_position:idx);
[medium_value medium_position_sub_smaller_add_1]=min(smaller_to_idx);
medium_position=medium_position_sub_smaller_add_1+smaller_position-1;


%------------------------------
%重設idx為第一山頂與波谷之間（因背景佔太多要微調二值化的值）
%------------------------------

idx_new=(smaller_position+medium_position)/2;

%------------------------------
%影像處理
%------------------------------

%平滑參數
%fspecial('average',xxx)等同ones(xxx,xxx)/sum(sum(ones(xxx,xxx)))
se_average=fspecial('average',3);

%銳利參數
%se_unsharp = fspecial('unsharp',1);%1 is max unsharp

%open參數
se_open = strel('disk',2);

%close參數
se_close = strel('disk',5);

for k=1:all_file_num_add2-2

    %平滑
    blurred_im_file(:,:,k) = imfilter(im_file(:,:,k),se_average,'replicate');

    %二值
    threshold_im_file(:,:,k)=(blurred_im_file(:,:,k)>idx_new);

    %open
    open_im_file(:,:,k)=imopen(threshold_im_file(:,:,k),se_open);
    
    %close
    close_im_file(:,:,k)=imclose(open_im_file(:,:,k),se_close);

    %最長邊界
    boundaries_B=boundaries(close_im_file(:,:,k));
    boundaries_B_cell_1=boundaries_B{1};
    [close_im_file_row close_im_file_column]=size(close_im_file(:,:,k));
    external_edge_im_file(:,:,k)=bound2im(boundaries_B_cell_1,close_im_file_row,close_im_file_column,min(boundaries_B_cell_1(:,1)),min(boundaries_B_cell_1(:,2)));

    %填滿
    fill_im_file(:,:,k)=imfill(external_edge_im_file(:,:,k),'holes');

end

%套至原圖
final_im_file=uint16(fill_im_file).*im_file;

%將影像做3D平滑處理

final_im_file=smooth3(final_im_file);

%取我們要看的第1至150張

final_im_file_150=final_im_file(:,:,1:150);

%------------------------------
%儲存變數final_im_file並消除其他無用變數以避免記憶體不足
%------------------------------

save('project20090622_final_im_file_150.mat','final_im_file_150');
clear;
load project20090622_final_im_file_150;

%------------------------------
%3D影像重建
%------------------------------

hiso = patch(isosurface(final_im_file_150,5),'FaceColor',[1,.75,.65],'EdgeColor','none');
lightangle(45,30); 
set(gcf,'Renderer','zbuffer');
lighting phong;
isonormals(final_im_file_150,hiso);
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);

view(45,30);

%程式結束