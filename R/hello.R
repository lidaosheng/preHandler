getData<-function(path,destdir){
  folders<-list.files(path)
  GBM_fpkm<-data.frame()
  #整合数据
  fd1<-folders[1]
  files_name<-list.files(paste(path,"/",fd1,sep = ""))
  files_name_gz<-files_name[grepl('.FPKM.txt',files_name)]
  if(length(files_name_gz)!=0){
    mydata<-read.table(gzfile(paste(path,'/',fd1,'/',files_name_gz,sep = "")))
    #列名
    names(mydata)<-c("ENSG_ID",fd1)
    GBM_fpkm<-mydata
  }


  for(fd in folders[2:length(folders)]){
    #解压文件
    files_name<-list.files(paste(path,"/",fd,sep = ""))
    files_name_gz<-files_name[grepl('.FPKM.txt',files_name)]
    if(length(files_name_gz)==0){
      next
    }
    mydata<-read.table(gzfile(paste(path,'/',fd,'/',files_name_gz,sep = "")))
    #列名
    names(mydata)<-c("ENSG_ID",fd)
    GBM_fpkm<-merge(GBM_fpkm,mydata,"ENSG_ID")
  }
  write.csv(GBM_fpkm,file = paste0(destdir,"/GBM_FPKM.csv"),row.names = FALSE)
}

#将csv文件读取并处理后返回结果,行样本，列基???
csvToEset<-function(fileDir){
  eset<-read.csv(fileDir,header = TRUE)
  tryCatch(
    {
      rownames(eset)<-eset[,1]
      eset<-eset[,2:ncol(eset)]
    },
    error=function(e){cat("fail to read goal csv file",conditionMessage(e),"\n\n")},
    finally={}
  )
  return(eset)
}

#turn esets that genes between which are different to same
comGenesEsets<-function(esetList){
  len=length(esetList)
  ifelse(
    len==0,
    stop("esetList is empty!"),
    {
      ifelse(len==1,
             return(esetList),
             {name<-lapply(esetList,function(x){colnames(x)})
             commGene<-name[[1]]
             for(i in 2:len){commGene<-intersect(commGene,name[[i]])} #get common genes between each eset
             esetlist2<-lapply(esetList,function(x){x[,commGene]}) #only keep commGene for each member of esetList
             rm(esetList,commGene,name,len)
             return(esetlist2)
             }
      )
    }
  )

}
reduceGeneNum<-function(eset){
  if(dim(eset)[2]<5000){return(eset)}
  eset <- eset[,order(apply(eset,2,mad), decreasing = T)[1:5000]]
  return(eset)
}

#remove the outliers and return the remain samples and labels
removeOutliers<-function(eset,label){
  sampleTree<-hclust(dist(eset),method="average")
  sizeGrWindow(12,9)
  par(cex=0.6)
  par(mar=c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
  #get user input
  again=TRUE
  while(again){
    tryCatch(
      { print("please input a height(numeric) to cut the tree")
        height<-scan("",what=numeric(0),nlines=1)
        abline(h=height,col="red")
        print("please input the min size to recognize a module")
        minSize<-scan("",what=integer(0),nlines=1)
        again=FALSE},
      error=function(e){print("there are some error，please input again...");again=TRUE},
      finally={}
    )
  }
  #cut tree
  clust = cutreeStatic(sampleTree, cutHeight = height, minSize = minSize)
  rm(sampleTree)
  again<-TRUE
  removeSample<-c(1:nrow(eset)) #init
  while(again){
    tryCatch(
      { print("please input the index of cluster to delete, press enter to complete:")
        removeClust<-scan("",what=integer(0),nlines = length(table(clust)))
        ifelse(
          #there is a bug to fix
          all(removeClust%in%c(0:length(table(clust))))==TRUE&length(removeClust)!=0,
          #表达式拼???
          {eps<-sapply(removeClust,function(x){paste("clust==",x,sep="")})
          eps<-parse(text=paste(eps,"",sep = "",collapse = "|"))
          removeSample<-eval(eps)
          rm(eps)
          again=FALSE},
          {print("some index is valid or empty input,please input again..")
            again=TRUE}
        )
      },
      error=function(e){print("输入有误，请重新输入...");again=TRUE},
      finally={}
    )

  }
  eset2<-eset[-which(removeSample),]
  label2<-label[-which(removeSample)]
  rm(eset,removeSample,again,clust)
  e_l<-list(eset=eset2,label=label2)
  return(e_l)
}

#a function to detect modules
moduleDetect<-function(eset,dissTOM){
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  minSize<-floor(dim(eset)[2]/100)
  if(minSize<4)
    minSize=4
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minSize)
  dynamicColors = labels2colors(dynamicMods)
  MEDissThres<-0.25
  # Call an automatic merging function
  merge = mergeCloseModules(eset, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Rename to moduleColors
  moduleColors = mergedColors
  return(moduleColors)
}

#返回训练好的moxing（实际是最后一折的）,以及十折交叉验证准确率
trainModel<-function(eset,label,model="NN"){
  if(model=="NN"){
    str = "nn <- nnet(label ~ .,data = trainset,size = 2,rang = 0.1,decay = 5e-4,maxit = 200,trace=F)"
    str2 = "acc <- getAcc(testset$label,predict)"
  }else if(model=="NB"){
    str = "nn <- NaiveBayes(label ~ .,data = trainset)"
    str2 = "acc <- getAcc(testset$label,predict$class)"
  }else{
    stop("参数model只能取值NN,NB..")
  }
  set.seed(100)
  eset<- as.data.frame(eset)
  label<-as.factor(label)
  maxs<-apply(eset,2,max)
  mins<-apply(eset,2,min)
  eset<-as.data.frame(scale(eset,center=mins,scale=maxs-mins))
  rm(mins,maxs)
  data<-cbind(eset,label=label)
  #十次
  #-----------------------------------------------------
  #----------------------------------------------------
  result<-sapply(1:10,function(x){
    #十折交叉
    folds<-createFolds(y=data$label,k=10)
    result1<-sapply(folds,function(x){
      #每次先选好训练集和测试集
      trainset<-data[-x,]
      testset<-data[x,]
      #训练网络
      # nn <- NaiveBayes(label ~ .,data = trainset)
      eval(parse(text=str))
      predict <- predict(nn,testset,type = "class")
      # acc <- getAcc(testset$label,predict$class)
      eval(parse(text=str2))
    })
    mean(result1)
  })
  print("十次十折交叉验证：")
  print(mean(result))
  # list1<-list(result=mean(result),nn=nn)
  return(mean(result))
}
#二分类，获取acc
getAcc<-function(label,predict){
  nt <- table(label,predict)
  acc<-0
  if(length(colnames(nt))==1){
    acc<-(table(label,predict)[colnames(table(label,predict)),colnames(table(label,predict))])/length(label)
    return(acc)
  }
  acc <- (nt[1,1]+nt[2,2])/(length(label))
  return(acc)
}

#验证用的方法，来看特征在模块中的位置
showCorPos<-function(eset,moduleColors,choose,label){
  colors<-unique(moduleColors)
  #获取choose所在模块颜色
  indexes<-match(choose,colnames(eset))
  choose_colors<-moduleColors[indexes]
  #获得各个模块表达谱子集,命名为color
  for(i in 1:length(colors)){
    str<-paste0(colors[i],'<-eset[,which(moduleColors=="',colors[i],'")]')
    eval(parse(text=str))
  }
  #获取各个模块和label的cor的降序基因名,命名为color_dec
  for(i in 1:length(colors)){
    str<-paste0('cor_',colors[i],'<-cor(',colors[i],',label)')
    eval(parse(text=str))
    str<-paste0('cor_',colors[i],'<-t(cor_',colors[i],')')
    eval(parse(text=str))
    str<-paste0('cor_',colors[i],'<-as.data.frame(cor_',colors[i],')')
    eval(parse(text=str))
    str<-paste0('cor_',colors[i],'<-sort(abs(cor_',colors[i],'),decreasing=T)')
    eval(parse(text=str))
    str<-paste0(colors[i],'_dec<-colnames(cor_',colors[i],')')
    eval(parse(text=str))
  }
  #获得各个模块的cor排名
  rank<-vector(length = length(choose),mode = "integer")
  ranks<-vector(length = length(choose),mode = "character")
  for(i in 1:length(choose)){
    str<-paste0('rank[i]<-match("',choose[i],'",',choose_colors[i],'_dec)')
    eval(parse(text=str))
    str<-paste0('len<-length(',choose_colors[i],'_dec)')
    eval(parse(text = str))
    str<-paste0(rank[i],'/',len)
    ranks[i]<-str
  }
  return(ranks)

}

#去冗余，每次去掉一个最差的特征
#终止条件：连续下降3次，或者单次下降2百分点
removeWF<-function(data,label,remainNum=2,model="NN"){
  #记录每次迭代次数，精度，去掉的特征
  len = ncol(data)-remainNum+1 #剩余迭代剩余次数+1
  iter = vector(mode = "integer",length = len)
  iter_acc = vector(mode = "numeric",length = len)
  iter_f = vector(mode = "character",length = len)
  #初始化数据，监控变量
  data1<-data
  count = 1
  isStop = 0
  #初次没有删除元素，但是也要记录
  acc<-trainModel(data1,as.factor(label),model) #记录初始精度
  iter[1]=0
  iter_f[1]="--"
  iter_acc[1]=acc
  #---------------------------------------------------------------
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores)
  clusterEvalQ(cl,library(caret))
  clusterEvalQ(cl,library(nnet))
  clusterEvalQ(cl,library(klaR))
  #-----------------------------------------------------------------
  while(len>1&&ncol(data1)>1){ #如果没有终止，一直迭代到剩下remainNum个特征
    accs<-parSapply(cl,1:ncol(data1),function(x){
      # accs<-sapply(1:ncol(data1),function(x){
      data2<-data1[,-x]
      acc_t<-trainModel(data2,as.factor(label),model)
    })
    #得到准确率提升最大的
    accs<-as.numeric(accs)
    remove_index=NULL
    remove_index<-which(accs==max(accs)) #那个去掉后，让整体性能提升最多的特征索引
    iter_f[count+1]<-colnames(data1)[remove_index[1]]
    iter_acc[count+1]<-accs[remove_index[1]]
    data1<-data1[,-remove_index[1]]
    len<-len-1
    iter[count+1]<-count
    count<-count+1
  }
  #--------------------------------------------------------
  stopCluster(cl)
  pos<-which(iter_acc==max(iter_acc))
  pos<-pos[length(pos)]#精度最高，长度最少
  rmGenes<-iter_f[1:pos]
  genelist<-colnames(data)
  x<-match(rmGenes,genelist)
  x<-x[-which(is.na(x))]
  if(length(x)>0)
    genelist<-genelist[-x]
  result<-list(genelist=genelist,acc=max(iter_acc),iter_f=iter_f,iter_acc=iter_acc)
  return(result)
}
#eset行样本，列特征
wgcnaPredict<-function(eset,label,stop_acc=1,model="NN",cor1=0.85){
  eset<-prepareData(eset,label,cor1)
  eset2<-scale(eset)#eset2是标准化的eset,仅用于聚类
  dissTOM<-1-cor(eset2)#相似矩阵化为相异矩阵，用于层次聚类
  moduleColors<-moduleDetect(eset2,dissTOM)#获取簇
  colors<-table(moduleColors)
  removeColors<-names(which(colors==1)) #移除只有一个元素的簇，否则会出错
  colors<-setdiff(names(colors),removeColors) #剩下模块的名字
  #和colors顺序一致
  colors<-lapply(colors,function(x){x1<-as.data.frame(eset[,which(moduleColors==x)])})
  #获取各个模块和label的cor的降序基因名,命名为color_dec
  #顺序与colors一致
  colors_dec<-lapply(colors,function(x){
    cor1<-cor(x,label)
    cor1<-t(cor1)
    cor1<-as.data.frame(cor1)
    cor1<-sort(abs(cor1),decreasing = T)
    cor1<-colnames(cor1)
  })
  rm(colors,eset2,dissTOM,removeColors)
  gc()
  #将各个模块cor第一名放进去，颜色顺序与colors同
  first<-sapply(colors_dec,function(x){x[1]})
  first<-unlist(first)

  #迭代替换基因
  print("Starting replace features in genelist ...")
  first<-replaceGene(first,colors_dec,eset,label,stop_acc,model)

  print("Starting remove least contribution feature ... ")
  result2<-removeWF(eset[,first],label,model=model)
  #用于测试目标基因在模块中的位置
  pos<-showCorPos(eset,moduleColors,result2$genelist,label)
  li <- list(result=result2,pos=pos)
  return(li)
}

#replace geneVector genes with gene in coresponding module,to reach the highest predict score
#first 由各个colors_dec第一个元素组成的基因列表
#colors_dec 某个模块基因降序排列（与label的cor）
#fast
replaceGene<-function(first,colors_dec,eset,label,end=1,model="NN"){
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores)
  clusterEvalQ(cl,library(caret))
  clusterEvalQ(cl,library(nnet))
  clusterEvalQ(cl,library(klaR))
  genelist<-as.character(first)
  acc<-trainModel(eset[,genelist],as.factor(label),model) #记录初始精度
  for(i in 1:length(genelist)){
    if(length(colors_dec[[i]])==1){next}
    # acc<-trainModel(eset[,genelist],as.factor(label),model) #记录初始精度
    #if the max acc bigger than end,stop the loop
    if(acc>=end)
      break
    accs<-parSapply(cl,colors_dec[[i]][-1],function(x){
      genelist[i]<-x
      acc2<-trainModel(eset[,genelist],as.factor(label),model)
    })
    #acc提升，替换
    if(max(accs)>acc){
      max_index<-1+which(accs==max(accs))[1]
      genelist[i]<-colors_dec[[i]][max_index]
      acc<-max(accs)
    }
  }
  stopCluster(cl)
  first<-genelist
  return(first)
}
#--------------------------------------------------------------------------
addGene<-function(geneVector,colors_dec,index,eset,label){
  for(i in 1:length(colors_dec)){
    geneVector1<-c(geneVector,co)
    sapply(colors_dec[i], function(x){
      geneVector1<-unique(c(x,geneVector))

    })


  }

}
#准备工作，数据清理，降维度
#eset行样本，列特征
#label为样本标签，integer格式
prepareData<-function(eset,label,cor1=0.85){
  #去掉缺失值--------------------------------------------------------------
  print("Removing feature with missing value...")
  if(length(is.na(eset))>0){
    missIndex<-apply(eset,2,function(x){length(which(is.na(x)))>0})
    missIndex<-as.integer(which(missIndex))
    if(length(missIndex)>0)
      eset<-eset[,-missIndex]
  }

  #去掉零方差------ --------------------------------------------------------
  print("Removing features with zero variance...")
  eset_p<-scale(eset)
  zerovar<-nearZeroVar(eset_p)
  if(length(zerovar)>0)
    eset<-eset[,-zerovar]

  if(dim(eset)[2]<4000){return(eset)}
  eset <- eset[,order(apply(eset,2,mad), decreasing = T)[1:4000]]

  #高相关过滤--------------------------------------------------------
  eset_p<-scale(eset)
  descrCorr<-cor(eset_p)
  highCorr<-findCorrelation(descrCorr,cor1)
  if(length(highCorr)>0)
    eset<-eset[,-highCorr]
  #t-test
  p.value.all.genes = apply(eset,2,function(x){t.test(x[which(label==1)],x[which(label==0)])$p.value})
  result.t_test<-subset(p.value.all.genes,p.value.all.genes<0.01)
  if(length(result.t_test)<500){
    result.t_test<-subset(p.value.all.genes,p.value.all.genes<0.5)
    if(length(result.t_test)<500){
      result.t_test<-p.value.all.genes
    }
  }
  eset<-eset[,names(result.t_test)]
  #离群样本
  #e_l<-removeOutliers(eset,label)
  return(eset)
}

