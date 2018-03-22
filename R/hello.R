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
#十折交叉验证
testbioPicker<-function(eset,label,model="NB",cor1=0.85,k=5,type="acc"){
  set.seed(100)
  data<-as.data.frame(cbind(eset,label=label))
  result<-sapply(1:3,function(x){
    #十折交叉
    folds<-createFolds(y=data$label,k=7)
    result1<-sapply(folds,function(x){
      #每次先选好训练集和测试集
      trainset<-data[-x,]
      testset<-data[x,]
      re<-wgcnaPredict(trainset[,-ncol(trainset)],trainset$label,cor1=cor1,model = "NB",k=k,type=type)
      result2<-trainModel(testset[,re$result$genelist],testset[,ncol(testset)],type = type)
      print(paste("nei bu :",re$result$acc," 外部交叉验证: ",result2,"--",dim(trainset)[1],"---",dim(testset)[1]))
      result2
    })
    mean(result1)
  })
  print(result)
  return(result)
}
trainModel4<-function(eset,label){
  # set.seed(100)
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
  result<-sapply(1:3,function(x){
    #十折交叉
    folds<-createFolds(y=data$label,k=4)
    result1<-sapply(folds,function(x){
      #每次先选好训练集和测试集
      trainset<-data[-x,]
      testset<-data[x,]
      #训练网络
      nn <- naiveBayes(label ~ .,data = trainset)
      predict <- predict(nn,testset[,-ncol(testset)],type = "class")
      acc <- getAcc(testset$label,predict)
      nn1 <- svm(label ~ ., data = trainset,cost=2)
      predict <- predict(nn1,testset[,-ncol(testset)],type = "class")
      acc1 <- getAcc(testset$label,predict)
      result2<-sapply(1:nrow(data),function(x){
        #每次先选好训练集和测试集
        trainset<-data[-x,]
        testset<-data[x,]
        predict <- knn(trainset[,-ncol(trainset)],testset[-length(testset)],trainset$label,k=7,prob=TRUE)
      })
      acc3 <- getAcc(label,result2,type = "bacc")
      acc<-(acc+acc1+acc3)/3
    })
    mean(result1)
  })
  print("十次十折交叉验证：")
  print(paste("dim(eset,label)----",dim(eset)[1],"---",mean(result)))
  # list1<-list(result=mean(result),nn=nn)
  return(mean(result))
}
trainModel3<-function(eset,label,k=5,type="acc"){
  eset<- as.data.frame(eset)
  label<-as.factor(label)
  maxs<-apply(eset,2,max)
  mins<-apply(eset,2,min)
  eset<-as.data.frame(scale(eset,center=mins,scale=maxs-mins))
  data<-as.data.frame(cbind(eset,label=label))
  #十次
  #-----------------------------------------------------
  #----------------------------------------------------
  #十折交叉
  predicts<-vector(length = nrow(data))
  result1<-sapply(1:nrow(data),function(x){
    #每次先选好训练集和测试集
    trainset<-data[-x,]
    testset<-data[x,]
    predict <- knn(trainset[,-ncol(trainset)],testset[-length(testset)],trainset$label,k=k,prob=TRUE)
  })
  acc <- getAcc(label,result1,type = type)
  return(acc)
}
trainModel2<-function(eset,label){
  set.seed(100)
  eset<- as.data.frame(eset)
  label<-as.factor(label)
  maxs<-apply(eset,2,max)
  mins<-apply(eset,2,min)
  eset<-as.data.frame(scale(eset,center=mins,scale=maxs-mins))
  data<-as.data.frame(cbind(eset,label=label))
  #十次
  result<-sapply(1:5,function(x){
    #十折交叉
    folds<-createFolds(y=data$label,k=10)
    result1<-sapply(folds,function(x){
      #每次先选好训练集和测试集
      trainset<-data[-x,]
      testset<-data[x,]
      #训练网络
      predict <- knn(trainset[,-ncol(trainset)],testset[,-ncol(testset)],trainset$label,k=7,prob=TRUE)
      acc <- getAcc(testset$label,predict,type = "bacc")
    })
    mean(result1)
  })
  return(mean(result))
}
#返回训练好的moxing（实际是最后一折的）,以及十折交叉验证准确率
trainModel<-function(eset,label,model="NN",type="acc"){
  set.seed(100)
  if(model=="NN"){
    str = "nn <- nnet(label ~ .,data = trainset,size = 2,rang = 0.1,decay = 15e-4,maxit = 350,trace=F)"
  }else if(model=="NB"){
    str = "nn <- naiveBayes(label ~ .,data = trainset)"
  }else if(model=="SVM"){
    str = "nn <- svm(label ~ ., data = trainset)"
  }else{
    stop("参数model只能取值NN,NB..")
  }
  # set.seed(100)
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
      predict <- predict(nn,testset[,-ncol(testset)],type = "class")
      acc <- getAcc(testset$label,predict,type)
    })
    mean(result1)
  })
  print("十次十折交叉验证：")
  print(paste("dim(eset,label)----",dim(eset)[1],"---",mean(result)))
  # list1<-list(result=mean(result),nn=nn)
  return(mean(result))
}
#二分类，获取acc
getAcc<-function(label,predict,type="acc"){
  nt <- table(label,predict)
  acc<-0
  if(type=="acc"){
    if(length(colnames(nt))==1){
      acc<-(table(label,predict)[colnames(table(label,predict)),colnames(table(label,predict))])/length(label)
      return(acc)
    }
    acc <- (nt[1,1]+nt[2,2])/(length(label))
    return(acc)
  }else if(type=="bacc"){
    bacc<-(nt[1,1]/(nt[1,1]+nt[1,2])+nt[2,2]/(nt[2,2]+nt[2,1]))/2
    return(bacc)
  }
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
removeWF<-function(data,label,remainNum=2,model="NN",k=5,type="acc"){
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
  # acc<-trainModel(data1,as.factor(label),model) #记录初始精度
  acc<-trainModel3(data1,label,k,type = type)
  # acc<-trainModel4(data1,label)
  iter[1]=0
  iter_f[1]="--"
  iter_acc[1]=acc
  #---------------------------------------------------------------
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores)
  clusterEvalQ(cl,library(caret))
  clusterEvalQ(cl,library(nnet))
  clusterEvalQ(cl,library(e1071))
  clusterEvalQ(cl,library(class))
  clusterEvalQ(cl,library(clv))
  #-----------------------------------------------------------------
  while(len>1&&ncol(data1)>1){ #如果没有终止，一直迭代到剩下remainNum个特征
    accs<-parSapply(cl,1:ncol(data1),function(x){
      # accs<-sapply(1:ncol(data1),function(x){
      data2<-data1[,-x]
      # acc_t<-trainModel(data2,as.factor(label),model)
      acc<-trainModel3(data2,label,k,type = type)
      # acc<-trainModel4(data2,label)
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
wgcnaPredict<-function(eset,label,stop_acc=1,model="NN",cor1=0.85,k=5,type="acc"){
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
    cor1<-sort(abs(cor1),decreasing = F)
    cor1<-colnames(cor1)
  })
  rm(colors,eset2,dissTOM,removeColors)
  gc()
  #将各个模块cor第一名放进去，颜色顺序与colors同
  first<-sapply(colors_dec,function(x){x[1]})
  first<-unlist(first)

  #迭代替换基因
  print("Starting replace features in genelist ...")
  first<-replaceGene(first,colors_dec,eset,label,stop_acc,model,k=k)

  print("Starting remove least contribution feature ... ")
  result2<-removeWF(eset[,first],label,model=model,k=k)
  #用于测试目标基因在模块中的位置
  #pos<-showCorPos(eset,moduleColors,result2$genelist,label)
  li <- list(result=result2)
  return(li)
}

#replace geneVector genes with gene in coresponding module,to reach the highest predict score
#first 由各个colors_dec第一个元素组成的基因列表
#colors_dec 某个模块基因降序排列（与label的cor）
#fast
replaceGene<-function(first,colors_dec,eset,label,end=1,model="NN",k=5,type="acc"){
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores)
  clusterEvalQ(cl,library(caret))
  clusterEvalQ(cl,library(nnet))
  clusterEvalQ(cl,library(e1071))
  clusterEvalQ(cl,library(class))
  clusterEvalQ(cl,library(clv))
  genelist<-as.character(first)
  # acc<-trainModel(eset[,genelist],as.factor(label),model) #记录初始精度
  acc<-trainModel3(eset[,genelist],label,k,type = type)
  # acc<-trainModel4(eset[,genelist],label)
  for(i in 1:length(genelist)){
    if(length(colors_dec[[i]])==1){next}
    #if the max acc bigger than end,stop the loop
    if(acc>=end)
      break
    accs<-parSapply(cl,colors_dec[[i]][-1],function(x){
      genelist[i]<-x
      # acc2<-trainModel(eset[,genelist],as.factor(label),model)
      acc2<-trainModel3(eset[,genelist],label,k,type=type)
      # acc2<-trainModel4(eset[,genelist],label)

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

  #与label的相关性过滤
  cor_label<-cor(eset,label)
  names(cor_label)<-colnames(eset)
  cor_label<-sort(abs(cor_label),decreasing = T)
  if(dim(eset)[2]<1000){return(eset)}
  eset<-eset[,names(cor_label)[1:1000]]

  # if(dim(eset)[2]<4000){return(eset)}
  # eset <- eset[,order(apply(eset,2,mad), decreasing = T)[1:4000]]

  #高相关过滤--------------------------------------------------------
  # eset_p<-scale(eset)
  # descrCorr<-cor(eset_p)
  # highCorr<-findCorrelation(descrCorr,cor1)
  # if(length(highCorr)>0)
  #   eset<-eset[,-highCorr]
  # #t-test
  # p.value.all.genes = apply(eset,2,function(x){t.test(x[which(label==1)],x[which(label==0)])$p.value})
  # q.value<-p.adjust(p.value.all.genes,method = "BH")
  # result.t_test<-subset(q.value,q.value<0.01)
  # if(length(result.t_test)<500){
  #   result.t_test<-subset(q.value,q.value<0.5)
  #   if(length(result.t_test)<500){
  #     result.t_test<-p.value.all.genes
  #   }
  # }
  # eset<-eset[,names(result.t_test)]
  #离群样本
  #e_l<-removeOutliers(eset,label)
  return(eset)
}
mcone<-function(eset,label,cor1=0.3,fm=1){
  corList<-abs(cor(eset,label))
  names(corList)<-colnames(eset)
  corList<-sort(corList,decreasing = T)
  #按相关性大小顺序排列
  corlist<-corList[which(corList>cor1)]
  genenames<-names(corlist)
  choose<-vector(length = length(genenames))
  choose[1]=genenames[1]
  len <- 1
  #筛选特征
  for(i in 2:length(genenames)){
    is_remove<-FALSE
    for(m in 1:len){
      if(cor(eset[,m],eset[,i])>corlist[i]/fm){
        is_remove<-TRUE
        break
      }
    }
    if(is_remove==FALSE){
      len<-len+1
      choose[len]<-genenames[i]
    }
  }
  choose<-choose[1:len]
  return(choose)
}
#=======================================================================
mcone3<-function(eset,label,r){
  micFC<-mine(eset,label)$MIC
  names(micFC)<-colnames(eset)
  choose<-micFC[which(micFC>=r)]
  Subset<-sort(choose,decreasing=T)
  Subset<-names(Subset)
  Subset<-match(Subset,colnames(eset))
  for(e in 1:length(Subset)){
    q=e+1
    while(q<=length(Subset)){
      if(mine(eset[,Subset[e]],eset[,Subset[q]])>micFC[Subset[q]])
        Subset<-Subset[-q]
      else
        q=q+1
    }
    e=e+1
  }
  return(colnames(eset)[Subset])
}
mcone2<-function(eset,label,cor1,fm){
  corList<-mine(eset,label)$MIC
  names(corList)<-colnames(eset)
  corList<-sort(corList,decreasing = T)
  #按相关性大小顺序排列
  corlist<-corList[which(corList>cor1)]
  genenames<-names(corlist)
  choose<-vector(length = length(genenames))
  choose[1]=genenames[1]
  len <- 1
  #筛选特征
  for(i in 2:length(genenames)){
    is_remove<-FALSE
    for(m in 1:len){
      if(mine(eset[,m],eset[,i])$MIC>corlist[i]/fm){
        is_remove<-TRUE
        break
      }
    }
    if(is_remove==FALSE){
      len<-len+1
      choose[len]<-genenames[i]
    }
  }
  choose<-choose[1:len]
  return(choose)
}
newone<-function(eset,label,cc){
  # cor_label<-cor(eset,label)
  cor_label<-mine(eset,label)$MIC
  names(cor_label)<-colnames(eset)
  cor_label2<-sort(abs(cor_label),decreasing = T)
  cor_label2<-names(cor_label2[cc])
}
#十折交叉验证
testNewOne<-function(eset,label,cc,cor1,fm){
  set.seed(100)
  data<-as.data.frame(cbind(eset,label=label))
  result<-sapply(1:10,function(x){
    #十折交叉
    folds<-createFolds(y=data$label,k=4)
    result1<-sapply(folds,function(x){
      #每次先选好训练集和测试集
      trainset<-data[-x,]
      testset<-data[x,]
      # choose<-mcone2(trainset[,-ncol(trainset)],trainset[,ncol(trainset)],cor1,fm)
      if(cc==0){cc=c(1:length(choose))}
      choose<-newone(trainset[,-ncol(trainset)],trainset[,ncol(trainset)],cc)
      re<-trainModel(testset[,choose],testset[,ncol(testset)])
    })
    mean(result1)
  })
  print(mean(result))
  return(result)
}
