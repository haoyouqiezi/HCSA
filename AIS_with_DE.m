function [] = AIS_with_DE()
clear all
close all
clc

global FES
global fun
global UB
global LB
global Ndim
global initial_flag
f_bias   = [-450,-450,-450,-450,-310,390,-180,-140,-330,-330,90,-460,-130,-300,120,120,120,10,10,10,360,360,360,260,260];
boundary = [-100,100; -100,100; -100,100; -100,100; -100,100; -100,100; 0,600; -32,32; -5,5; -5,5; -0.5,0.5; -pi,pi; -3,1; -100,100; -5,5; -5,5; -5,5; -5,5; -5,5; -5,5; -5,5; -5,5; -5,5; -5,5; -5,5];
Ndim = 10;                     %维数
mFES = 10000*Ndim;             %资源数
Np = Ndim*10;                  %个体数
Nr = 0.05;                     %基因重组率
m  = 0.3*Ndim;                 %基因重组的维数
% LB = -100*ones(1,Ndim);        %下边界
% UB = 100*ones(1,Ndim);         %上边界
H  = Np;                       %F记录的个体数
wf = 0.1315;                   %F的变异因子
p  = 0.2;                      %优秀个体的百分比
pp = 0.1;%3                    %确定变异的维数
stag_record = [];              %用于记录每个个体停滞的代数
stag_num    = 6;               %停滞超过stag_num的个体将被提出

   fun = 1;
for run=1:300
    LB = boundary(fun,1)*ones(1,Ndim);         %下边界
    UB = boundary(fun,2)*ones(1,Ndim);         %上边界
    
    ab = rand(Np,Ndim).*(ones(Np,1)*(UB-LB))+ones(Np,1)*LB;    %初始化种群
    initial_flag=0;
    for i=1:Np
        ab(i,Ndim+1)=benchmark_func(ab(i,1:Ndim),fun);
    end
    FES = Np;
    M_F= 0.4*ones(1,H);     %初始化M_F
    gen = 1;                %存储更新的代数
    k = 1;                  %用来存放该更新的F值的角标
    A = [];                 %A哟昂来存放呗更新过了的个体
    stag_record = zeros(1,Np); %记录个体停滞的代数
    %stag_pop    = [];          %记录停滞的种群
    trace = [];
    gradient = [];
    while(FES<mFES)         %资源还没有使用完时
        %%
        %基因重组部分
        %first:先进行基因重组
        X = [];
        time = 1;
        while(time<Np*Nr)   %Nr对应的是基因重组率
            aba = [ab;A];
            r = randperm(Np);
            r1= r(1);
            mm = size(aba,1);
            r = randperm(mm);
            r2= r(1); 
            X(1,1:Ndim+1) = aba(r1,1:Ndim+1);
            X(3,1:Ndim+1) = X(1,1:Ndim+1);
            X(2,1:Ndim+1) = aba(r2,1:Ndim+1);
            X(4,1:Ndim+1) = X(2,1:Ndim+1);
            MM = randperm(Ndim);
            MM = MM(1:m);
            for i=1:Ndim
                if(find(i==MM))
                    alpha = rand();
                    bata  = rand();
                    X(3,i) = alpha*X(1,i) - (1-alpha)*X(2,i);
                    X(4,i) = bata*X(2,i) - (1-bata)*X(1,i);
                end
            end
            %边界检测
            X(3,find(X(3,:)>UB(1)))=UB(1);
            X(3,find(X(3,:)<LB(1)))=LB(1);
            X(4,find(X(4,:)>UB(1)))=UB(1);
            X(4,find(X(4,:)<UB(1)))=UB(1);
            
            X(3,Ndim+1) = benchmark_func(X(3,1:Ndim),fun);  FES = FES+1;
            X(4,Ndim+1) = benchmark_func(X(4,1:Ndim),fun);  FES = FES+1;
            [A,best] = add_to_A(X(1,:),X(3,:),X(4,:),A);                    %在其子代中排序寻找个体加入A
            if(best==1)
                ab(r1,:) = X(best,:);
                stag_record(r1) = stag_record(r1)+1;           %停滞数加1
            else
                gradient(r1,1:Ndim) = X(best+1,1:Ndim) - ab(r1,1:Ndim);  %记录梯度信息
                ab(r1,:) = X(best+1,:);
                stag_record(r1) = 0;                           %停滞数复位
            end
            time = time+1;
        end
        %%
        %基因重组部分结束
        %DE变异公式的使用
        S_F = [];      %记录F的更新值
        dif = [];      %记录更新的差值
        for i=1:Np
            %使用DE的F进行变异
            %F的变异方式
            F = M_F + wf .* tan(pi * (rand(1, Np) - 0.5));
            pos = find(F <= 0);
            while ~ isempty(pos)
                F(pos) = M_F(pos) + wf .* tan(pi * (rand(1,length(pos)) - 0.5));
                pos = find(F <= 0);
            end
            F = min(F, 1);
            [val,index] = sort(ab(:,Ndim+1));
            ab_sort  = ab(index,:);
            percent  = randperm(floor(Np*p));         %p指的是优秀个体中的百分比
            ab_pbest = ab_sort(percent(1),:);         %在百分比的优秀个体中选择pbest
            r = randperm(Np);                         %在Np个种群个体中选择r1
            r1=r(1);
            x1=r(2);
            y1=r(3);
            z1=r(4);
            p1=r(5);
            q1=r(6);
            aba = [ab;A]; 
            mm  = size(aba,1);                        %在ab和A中整合的个体中选择r2
            r = randperm(mm);
            r2= r(1);
            x2= r(2);
            y2= r(3);
            z2= r(4);
            p2= r(5);
            q2= r(6);
            %v_ab(1,1:Ndim) = ab(i,1:Ndim)+F(i)*(ab_pbest(1:Ndim)-ab(i,1:Ndim))+F(i)*(ab(r1,1:Ndim)-aba(r2,1:Ndim));  %变异公式
            %v_ab(1,1:Ndim) = ab(i,1:Ndim)+rand*(ab_pbest(1:Ndim)-ab(i,1:Ndim))+rand*(ab(r1,1:Ndim)-aba(r2,1:Ndim));  %变异公式
            v_ab(1,1:Ndim) = ab(i,1:Ndim)+F(i)*(ab_pbest(1:Ndim)-ab(i,1:Ndim))+F(i)*(ab(r1,1:Ndim)-aba(r2,1:Ndim))+F(i)*(ab(x1,1:Ndim)-aba(x2,1:Ndim));
            %v_ab(1,1:Ndim) = ab(i,1:Ndim)+F(i)*(ab_pbest(1:Ndim)-ab(i,1:Ndim))+F(i)*(ab(r1,1:Ndim)-aba(r2,1:Ndim))+F(i)*(ab(x1,1:Ndim)-aba(x2,1:Ndim))+F(i)*(ab(y1,1:Ndim)-aba(y2,1:Ndim));%+F(i)*(ab(z1,1:Ndim)-aba(z2,1:Ndim));%+F(i)*(ab(p1,1:Ndim)-aba(p2,1:Ndim));
            %v_ab(1,1:Ndim) = ab(r1,1:Ndim) + (2*rand-1)*(ab(r1,1:Ndim)-aba(r2,1:Ndim));
            %v_ab(1,1:Ndim) = ab(i,1:Ndim)+F(i)*(ab_pbest(1:Ndim)-ab(i,1:Ndim));  %变异公式
            %利用原变异的方式确定变异的维数
            [index1,val] = find(ab(:,Ndim+1)==max(ab(:,Ndim+1)));
            [index2,val] = find(ab(:,Ndim+1)==min(ab(:,Ndim+1)));
            if(ab(index1(1),Ndim+1)==ab(index2(1),Ndim+1))
                f_n = 0.001;
            else
                f_n   =( ab(i,Ndim+1)-ab(index2(1),Ndim+1))/(ab(index1(1),Ndim+1)-ab(index2(1),Ndim+1));        %归一化操作
            end
             %fprintf('ab(i,Ndim+1)=%g max=%g min=%g\n',ab(i,Ndim+1),ab(index1(1),Ndim+1),ab(index2(1),Ndim+1));
%             if(f_n==0) 
%                 f_n=1; 
%             end
            alpha = exp(-1*pp*f_n);
            MM = floor(alpha*Ndim+1);                    %确定的变异维数
            if MM>=Ndim
                MM=Ndim;
            end
            m_rand = randperm(Ndim);
            M_Ndim = m_rand(1:MM);
            test_ab= ab(i,1:Ndim);
            %基因变异
            for j=1:Ndim
                if(find(j==M_Ndim))
                    test_ab(j) = v_ab(j);
                end
            end
            %边界检测
            test_ab(1,find(test_ab(1,:)>UB(1)))=UB(1);
            test_ab(1,find(test_ab(1,:)<LB(1)))=LB(1);
            
            test_ab(1,Ndim+1) = benchmark_func(test_ab(1,1:Ndim),fun);
            FES = FES+1;
            if(test_ab(Ndim+1)<ab(i,Ndim+1))
                dif=[dif;abs(test_ab(1,Ndim+1)-ab(i,Ndim+1))];
                A = [A;ab(i,:)];
                gradient(i,1:Ndim) = test_ab(1,1:Ndim) - ab(i,1:Ndim);    %记录梯度信息
                ab(i,:) = test_ab(1,:);
                S_F=[S_F;F(i)];
                stag_record(i) = 0;                         %复位停滞数
            else
                stag_record(i) = stag_record(i)+1;          %停滞数加1
            end
        end
        %运用DE公式的部分完，更新所使用的F的值
        if(~isempty(S_F))
            if(k>H)
                k=1;
            end
            sum_dif = sum(dif);
            dif_val = dif / sum_dif;
            meanw_SF=sum(dif_val.* S_F .*S_F)/ sum(dif_val.* S_F);
            M_F(k)=meanw_SF;
            k = k+1;
        end
        %%
        %DE的部分结束，开始自适应烟花搜索的部分
        %[ab] = fire_Search(ab,Np);
        [mmm,nnn] = size(A);
        if(mmm>Np)
            rrr = randperm(mmm);
            A = A(rrr(1:Np),:);
        end
        %[ab] = big_fire_Search(ab,Np);
        
        
        stag_index = [];
        stag_index_num = 0;
        for i_stag = 1:Np                                   %把停滞超过5带的个体
            if(stag_record(i_stag)>stag_num)
                      
%                     [ab(i_stag,:),F] = gradient_Search(ab(i_stag,:),gradient(i_stag,:));     %停滞的个体进行梯度搜素
                    %if(F==1) 
%%                        [ab(i_stag,:)] = Learn_to_best(ab,ab(i_stag,:),Np);
%此处对应的为基因konckout对应的对其进行隐藏
%%
%                         ab(i_stag,:) = fire_Search_singel(ab(i_stag,:));
                    %end
                    stag_record(i_stag) = 0;
% %                     ab(i_stag,:) = fire_Search_singel(ab(i_stag,:));
% %                   [ab_val,ab_index] = sort(ab(:,Ndim+1));
% %                   if(ab_index(1)==i_stag)
% %                       ab(i_stag,:) = fire_Search_singel(ab(i_stag,:));
% %                       stag_record(i_stag) = 1;
% %                   else
% %                       Y = ab(ab_index(1),:);
% %                       r = randperm(Np);
% %                       r1= r(1);
% %                       %r = randperm(size(A,1));
% %                       r2= r(2);
% %                       ab(i_stag,1:Ndim) =  Y(1,1:Ndim)+rand*(ab(r1,1:Ndim)-ab(r2,1:Ndim));
% %                       ab(i_stag,find(ab(i_stag,:)>UB(1,1))) = UB(1,1);
% %                       ab(i_stag,find(ab(i_stag,:)<LB(1,1))) = LB(1,1);
% %                       ab(i_stag,Ndim+1) = cec14_func(ab(i_stag,1:Ndim)',fun);
% %                       FES = FES+1;
% %                       stag_record(i_stag) = 0;
% %                   end
            end
        end
%         [ab] = fire_Search(ab,Np);
        [g_val g_index] = sort(ab(:,Ndim+1));
        for i=1:0.05*Np
%             [ab(g_index(i),:),F] = gradient_Search(ab(g_index(i),:),gradient(g_index(i),:));
%             if(F)
                ab(g_index(i),:) = fire_Search_singel(ab(g_index(i),:));
%             end
        end
        
        %%
        %记录当前值
        Z = [ab];
        [index,val] = find(Z(:,Ndim+1)==min(Z(:,Ndim+1)));
        trace(gen,1) = gen;
        trace(gen,2) = Z(index(1),Ndim+1)-f_bias(fun);
        trace(gen,2);
        gen = gen+1;
    end
    
%     figure
%     plot(trace(500:end,1),trace(500:end,2));
    [m_trace,n_trace] = size(trace);
    %trace(m_trace,2)    %-100*fun;
    
    [index_SOV,val_SOV] = find(ab(:,Ndim+1)==min(ab(:,Ndim+1)));
    if(mod(run,30))
        SOV(mod(run,30))=ab(index_SOV(1),Ndim+1)-f_bias(fun);   %-fun*100;
        if SOV(mod(run,30))<1e-8
            SOV(mod(run,30))=0;
        end
    else
        SOV(30)=ab(index_SOV(1),Ndim+1)-f_bias(fun);              %-fun*100;
        if SOV(30)<1e-8
            SOV(30)=0;
        end
    end
    
    Best=min(SOV);
    Worst=max(SOV);
    Median=median(SOV);
    Mean=mean(SOV);
    Std=std(SOV);
    if(mod(run,30)==0)
        fprintf('runtime=%g function=%g dimension=%g\n',run,fun,Ndim);
        fprintf('Best=%g Worst=%g Median=%g Mean=%g Std=%g\n',Best,Worst,Median,Mean,Std);
    end
    if(mod(run,30)==0)
        fun = fun+1;
        SOV = [];
    end
end
if(run>0)
    run
end