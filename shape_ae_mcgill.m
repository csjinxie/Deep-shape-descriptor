clear
clc
str = 'D:\Data\McGillDBAll-255\';
str1='D:\Data\McGillDBAll-HKS-255\';
dir1 = dir(str); 
Allname = struct2cell(dir1);
Allname = Allname(:,3:length(Allname));
temp=1;
shape_feats=[];
flag=[];
for i=1:length(Allname)
    name=strcat(str, Allname{1,i});
    dir2=dir(name);
    Allname2 = struct2cell(dir2);
    for j=1:length(Allname2)-2
        load (strcat(str1, Allname{1,i}, '\', Allname2{1,j+2}(1:end-3), 'mat'));
        shape_feats{temp}=HKS;
        flag(temp)=i;
      types(temp).grp =Allname(1,i);
      types(temp).idx = i;
      temp=temp+1;
    end    
end


for i=1:length(Allname)
    [va,in]=find(flag==i);
    train_index(i)=in(1);
    class_num(i)=length(in);
end
train_num=10;
whole=[1:size(shape_feats,2)];
sample_index=[];
for i=1:length(train_index)
    sample_index=[sample_index  train_index(i): train_index(i)+train_num-1];
end
[index1, index2]=ismember(whole, sample_index);
test_index=whole(index1==0);



for i=1:size(shape_feats{1},1)
    temp=[];
    for j=1:size(shape_feats,2)
        temp=[temp;shape_feats{j}(i,:)];  
    end
   multi_histogram(:,:,i)=temp'; 
end


label=flag(sample_index);
lambda=0.001;
n=[128 1000 500 250 30];
num_iter=[30,30];
batch_size=[10,20];
err_type=2;

ff=1;
%for i=1:size(shape_feats{1},1)
for i=1:3:78
  train_data=multi_histogram(:,sample_index,i);
test_data=multi_histogram(:,test_index,i);  
[Wb, RBM_error, BP_error] = train_AE ( train_data, train_data, label,lambda, test_data, ...
    n, num_iter, batch_size,err_type );
%Y(:,:,i) = AE_forward (  multi_histogram(:,:,i), Wb, n );
Y(:,:,ff) = AE_forward (  multi_histogram(:,:,i), Wb, n );
ff=ff+1;
end


feature2=[];
for i=1:size(test_index,2)
    tm=[];
    for j=1:size(Y,3)
        tm=[tm;Y(:,test_index(i),j)];
    end
    feature2(i,:)=tm';   
end
d2=[];
for i=1:size(test_index,2)
        for j=i:size(test_index,2)
       d2(i,j) = norm(feature2(i,:)- feature2(j,:),'fro');
      end
end
cor1=d2+d2';
m = size(cor1,1);
yes = 0;
recall = [];
precision= [];
for i = 1:m
    r1 =[];
    p1 = [];
    %query = index2group(i,types);  %query grp name
    query=flag(test_index(i));
    score = cor1(i,:);
    [d idx] = sort(score);
    yes = 0;    
    for j = 1:m
            index = idx(j);
            %retrive = index2group(index,types);
            retrive=flag(test_index(index));
            if query==retrive
                yes = yes +1;
            end
            r1 = [r1 yes/length(find(flag(test_index)==query))];
            p1 = [p1 yes/j];
    end
    
   recall(i,:) = r1;
   precision(i,:) = p1;
end
 pre1 = sum(precision,1)/size(precision,1);
 rec1 = sum(recall,1)/size(recall,1);
  plot(rec1,pre1)
