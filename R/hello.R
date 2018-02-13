#读取TCGA文件
readTCGA<-function(path){
  options(stringsAsFactors=F)
  folders<-list.files(path)
  new_fpkm<-data.frame()
  #整合数据
  fd1<-folders[1]
  files_name<-list.files(paste(path,"/",fd1,sep=""))
  files_name_gz<-files_name[grepl('.FPKM.txt',files_name)]
  mydata<-read.table(gzfile(paste(path,"/",fd1,"/",files_name_gz,sep="")))
  #列名
  names(mydata)<-c("ENSG_ID",fd1)
  new_fpkm<-mydata

  for(fd in folders[2:length(folders)]){
    #解压文件
    files_name<-list.files((paste(path,"/",fd,sep = "")))
    files_name_gz<-files_name[grepl('.FPKM.txt',files_name)]
    mydata<-read.table(gzfile(paste(path,"/",files_name_gz,sep = "")))
    #给出列名
    names(mydata)<-c("ENSG_ID",fd)
    #整合到一个数据框
    new_fpkm<-merge(new_fpkm,mydata,by="ENSG_ID")
  }
  #输出
  write.csv(new_fpkm,"f:/data/output/output.csv",row.names = F)
}

#path为GEO的数据文件.gz的本地路径，输出行基因，列样本矩阵
getData<-function(path){
  library(GEOquery)
  options(stringsAsFactors = F)
  gset<-getGEO(filename = path)
  eset<-exprs(gset)
  feadata<-featureData(gset)
  feadata<-feadata@data
  if(length(which(feadata$ID!=rownames(eset)))==0){
    rownames(eset)=feadata$`Gene Symbol`
    eset<-eset[-which(rownames(eset)==""),] #去掉空的基因（即探针对应不上基因的条目）
    result<-table(rownames(eset))
    #获取重复基因名
    dup<-names(which(result>1))
    #基因的表达值合并后
    for(x in dup){
      #获取重复值为x的所有行
      set<-eset[which(rownames(eset)==x),]
      #将行叠加到一起
      temp <- 1.0*colSums(set)/nrow(set)
      #将eset中x的行都删掉，然后将set第一行补上去
      eset<-eset[-which(rownames(eset)==x),]
      eset<-rbind(eset,temp)
      rownames(eset)[nrow(eset)]<-x
    }

  }else{
    stop("eset的探针名和注释文件不对应")
  }
  result<-list(feadata=feadata,eset=eset)
  return(result)

}
#富集分析评估
enrichmentAssess<-function(eset,moduleColors){
  background<-colnames(eset)
  modules<-unique(moduleColors)
  all.genes<-apply(as.data.frame(modules),1,function(x){
    module<-colnames(eset)[which(moduleColors==x)]
    return(module)
  })

}
celToExprs<-function(fileDir){
  # filters <- matrix(c("CEL file", ".[Cc][Ee][Ll]", "All", ".*"), ncol = 2, byrow = T)
  # cel.files <-tk_choose.files(fileDir,caption = "Select CELs", multi = TRUE,filters = filters, index = 1)
  # data.raw <- ReadAffy(filenames = cel.files)
  # rm(filters,cel.files)
  setwd(fileDir)
  data.raw<-ReadAffy()
  #RMA标准???
  rma<- rma(data.raw)
  #获取表达矩阵
  eset<-exprs(rma)
  rm(rma)
  #获取至少在一个样本中表达的矩???
  data.mas5calls<-mas5calls(data.raw)
  rm(data.raw)
  eset.mas5calls<-exprs(data.mas5calls)
  AP <- apply(eset.mas5calls, 1, function(x)any(x=="P"))
  rm(data.mas5calls)
  rm(eset.mas5calls)
  present.probes<-names(AP[AP])
  rm(AP)
  eset<-eset[present.probes,]
  if(length(which(is.na(eset)))!=0) stop("有缺失???")
  rm(present.probes)
  #处理方差???0的基因，缺失???
  zerovar<-nearZeroVar(t(eset))
  eset2<-t(eset)
  rm(eset)
  if(length(zerovar)!=0){
    print("有方差为0")
    eset2<-eset2[,-zerovar]
  }
  rm(zerovar)
  write.csv(eset2,file=paste(fileDir,"/eset2.csv",sep = ""))
  return(eset2)
}
#重载celToExps方法，批量处???

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
  print(table(clust))
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
#build a network,eset has its labels
buildNetwork<-function(eset){
  again=TRUE
  powers<-NULL
  while(again){
    tryCatch(
      {print("please input a set of softpowers,press enter to complete???")
       powers<-scan("",what=integer(0))
       #powers should be all positive,and with length>0
       ifelse(length(powers)>0&all(powers>0)==TRUE,
              {sft_p = pickSoftThreshold(eset, powerVector = powers, verbose = 5)
              #draw a plot
              sizeGrWindow(9, 5)
              par(mfrow = c(1,2));
              cex1 = 0.9;
              plot(sft_p$fitIndices[,1], -sign(sft_p$fitIndices[,3])*sft_p$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
              text(sft_p$fitIndices[,1], -sign(sft_p$fitIndices[,3])*sft_p$fitIndices[,2],labels=powers,cex=cex1,col="red");
              abline(h=0.90,col="red")
              plot(sft_p$fitIndices[,1], sft_p$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
              text(sft_p$fitIndices[,1], sft_p$fitIndices[,5], labels=powers, cex=cex1,col="red")
              again<-FALSE},
              {print("powers should be all positive,and with length>0..")
               again<-TRUE}
          )
       },
       error=function(e){print("there are power which are valid,please input again");again=TRUE}
    )
  }
  #type a best power
  again=TRUE
  power=NULL
  while(again){
    tryCatch(
      {print("choose a best power???")
       power<-scan("",what=integer(0),nlines = 1)
       ifelse(power%in%powers&length(power)==1,again<-FALSE,{print("a vaild power,try again..");again<-TRUE})
      },
      error=function(e){print("a valid input,try again..");again=TRUE}
    )
  }
  #build adjacency matrix
  adjacency_p<-adjacency(eset,power=power)
  TOM_P<-TOMsimilarity(adjacency_p)
  dissTOM_P<-1-TOM_P
  return(dissTOM_P)
}
#a function to detect modules
moduleDetect<-function(eset,dissTOM){
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
  ## I'll add a check for the input valid later
  print("input a minSize(integer) for module detect:")
  minSize<-scan("",what = integer(0),nlines = 1)
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minSize)
  dynamicColors = labels2colors(dynamicMods)
  print(class(dynamicColors))

  MElist<-moduleEigengenes(eset,colors=dynamicColors)
  MEs <- MElist$eigengenes
  if(length(is.na(MEs$MEgrey))>0){
    MEs$MEgrey<-NULL
  }
  MEDiss<-1-cor(MEs)
  METree<-hclust(as.dist(MEDiss),method="average")
  #Graphical the result
  sizeGrWindow(7,6)
  plot(METree,main="Clustering of module eigengenes(Normal)")
  MEDissThres = 0.25
  print("input a number to cut the tree(recommend 0.25):")
  MEDissThres<-scan("",what = numeric(0),nlines = 1)
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(eset, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  #绘制融合???(Dynamic Tree Cut)和融合后(Merged dynamic)的聚类图
  #sizeGrWindow(12, 9)
  #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  #dev.off()
  # 只是绘制融合后聚类图
  plotDendroAndColors(geneTree,mergedColors,"Merged dynamic",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  #5.结果保存
  # Rename to moduleColors
  moduleColors = mergedColors
  return(moduleColors)
}
#模块去冗余
moduleER<-function(data,moduleColors,modules){
  lapply(moduleColors,function(x){genes<-colnames(data)[which(moduleColors==x)]})

}

#calculate the module preservation between two data set
preservation<-function(eset1,eset2,dynamicColors_p1){
  setLabels=c("ref","test")
  multiExpr=list(ref=list(data=eset1),test=list(data=eset2))
  rm(eset1,eset2)
  multiColor=list(ref=dynamicColors_p1)
  rm(dynamicColors_p1)
  nSets=2
  system.time( {
    mp = modulePreservation(multiExpr, multiColor,
                            referenceNetworks = 1,
                            nPermutations = 200,
                            randomSeed = 1,
                            verbose = 3)}
  )
  save(mp,file="./mp.rda")
  ref = 1
  test = 2
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
  write.csv(signif(statsZ[, "Zsummary.pres"], 2),file = "zscore.csv")
  print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

}

moduleToCyto<-function(eset,moduleColors,module,hasTOM,TOM,power){
  if(hasTOM!=TRUE){
    TOM=TOMsimilarityFromExpr(eset, power=power)
  }
  inModule = is.finite(match(moduleColors, module))
  modTOM<-TOM[inModule,inModule]
  probes<-colnames(eset)
  modProbes = probes[inModule];
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("edges-", module, ".txt", sep=""),
                                 nodeFile = paste("nodes-", module, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = NULL,
                                 nodeAttr = moduleColors[inModule]);
}

allModuleToCyto<-function(eset,moduleColors,hasTOM,TOM,power){
  if(hasTOM!=TRUE){
    TOM=TOMsimilarityFromExpr(eset, power=power)
  }
  modules = as.character(names(table(moduleColors)))
  len = length(modules)
  sapply(modules,function(x){moduleToCyto(eset,moduleColors,x,T,TOM,NULL)})
}

trainData<-function(eset,label){
  set.seed(500)
  maxs<-apply(eset,2,max)
  mins<-apply(eset,2,min)
  data<-as.data.frame(scale(eset,center=mins,scale=maxs-mins))
  rm(eset)
  gc()
  sample_data <- sample(2,nrow(data),replace=TRUE,prob=c(0.7,0.3))
  label<-as.factor(label)
  tra_label<-label[sample_data==1]
  test_label<-label[sample_data==2]
  tra_data<-data[sample_data==1,]
  test_data<-data[sample_data==2,]
  #训练???
  data2<-cbind(tra_data,label8=tra_label)
  #测试???
  data3<-cbind(test_data,label8=test_label)
  #参数选择
  min1<-1000
  index<-1
  for (i in 1:(ncol(data3)-1)){
    test_model <- randomForest(label8~.,data=data2,mtry=i)
    err <- mean(test_model$err)
    #print(err)
    if(err<min1){index=i;min1=err}
  }
  print(paste("the smallest err index:",index," min-err:",min1))
  #mtry<-scan("",what = integer(0),nlines = 1)
  mtry<-index
  tran_model <- randomForest(label8~.,data=data2,mtry=mtry,ntree=500)
  plot(tran_model)
  print("choose a number of tree:")
  ntree<-scan("",what = integer(0),nlines = 1)
  abline(v=ntree,col="red")
  tran_model <- randomForest(label8~.,data=data2,mtry=mtry,ntree=ntree)
  print(table(actual=data3$label8,predicted=predict(tran_model,newdata = data3,type = "class")))
  return(tran_model)
}

#返回训练好的神经网络（实际是最后一折的）,以及十折交叉验证准确率
trainModelNN<-function(eset,label){
  eset<- as.data.frame(eset)
  label<-as.factor(label)
  maxs<-apply(eset,2,max)
  mins<-apply(eset,2,min)
  eset<-as.data.frame(scale(eset,center=mins,scale=maxs-mins))
  rm(mins,maxs)
  data<-cbind(eset,label=label)
  #十次
  result = vector(length=10)
  for(m in 1:10){
    #十折交叉
    folds<-createFolds(y=data$label,k=10)
    result1 = vector(length=10)
    for(i in 1:10){
      #每次先选好训练集和测试集
      trainset<-data[-folds[[i]],]
      testset<-data[folds[[i]],]
      #训练网络
      nn <- nnet(label ~ .,data = trainset,size = 2,rang = 0.1,decay = 5e-4,maxit = 200,trace=F)
      predict <- predict(nn,testset,type = "class")
      acc <- getAcc(testset$label,predict)
      result1[i] <- acc
    }
    result[m]=mean(result1)
  }
  print("十次十折交叉验证：")
  print(mean(result))
  list1<-list(result=mean(result),nn=nn)
  return(list1)
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
relateMT<-function(eset,moduleColors,label){
  # Define numbers of genes and samples
  nGenes = ncol(eset);
  nSamples = nrow(eset);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(eset, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, label, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(label),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  return(moduleTraitCor)

}
#选出感兴趣的模块
chooseModuleByCor<-function(moduleTraitCor,threshold=0.7){
  index<-which(abs(moduleTraitCor)>threshold)#选出相关系数大于threshold的模块
  chooseModule<-rownames(moduleTraitCor)[index] #MEpink
  return(chooseModule)
}
#在感兴趣的模块中挑选出首批基因集合用于下一步的特征选择
getFirstGeneSet<-function(moduleList){



}
#去冗余，去掉模块多重共线基因
removeCL<-function(data,moduleColors,module=NULL){
  if(ncol(data)!=length(moduleColors)){
    stop("data特征个数与moduleColor长度不匹配~")
  }
  if(!is.null(module))
    data = data[,which(moduleColors==module)]
  comboInfo<-findLinearCombos(data)
  if(length(comboInfo$remove)!=0){
    data<-data[,-comboInfo$remove]
  }
  return(data)
}
#去冗余，每次去掉一个最差的特征
#终止条件：连续下降3次，或者单次下降2百分点
removeWF<-function(data,label,remainNum=2){
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
  list1<-trainModelNN(data1,as.factor(label)) #记录初始精度
  acc<-list1$result #去掉特征前的准确率
  iter[1]=0
  iter_f[1]="--"
  iter_acc[1]=acc
  #-----------------------------------------------------------------
  while(len>1){ #如果没有终止，一直迭代到剩下remainNum个特征
    index<-as.data.frame(c(1:ncol(data1))) #1-特征总数，将向量化为数据框
    accs<-apply(index,1,function(x){
      data2<-data1[,-x]
      acc_t<-trainModelNN(data2,as.factor(label))$result
    })
    #得到准确率提升最大的
    accs<-as.numeric(accs)
    #diff<-accs-acc #去掉每个特征后，性能提升情况
    remove_index=NULL
    #-----------------------------------------------------------------------------------
    # if(all(diff<0)){ #性能全部下降
    #   iter_f[count+1]<-'empty'
    #   iter_acc[count+1]<-acc #精度是之前的精度，因为没有移除
    # }else{#如果有一个是提升的，重置isStop
    remove_index<-which(accs==max(accs)) #那个去掉后，让整体性能提升最多的特征索引
    iter_f[count+1]<-colnames(data1)[remove_index]
    iter_acc[count+1]<-accs[remove_index]
    data1<-data1[,-remove_index]
    #acc<-max(accs)
    len<-len-1
    iter[count+1]<-count
    count<-count+1
    #-----------------------------------------------------------------------------------
    #删除特征
    #remove_index<-which(accs==max(accs)) #那个去掉后，让整体性能提升最多的特征索引
    #记录
    # iter[count+1]<-count
    # iter_f[count+1]<-colnames(data1)[remove_index]
    # #iter_acc[count+1]<-accs[remove_index]
    # #为下次迭代准备
    # len<-len-1
    # count<-count+1
    # data1<-data1[,-remove_index]
    # acc<-max(accs)
  }
  result<-list(iter=iter,iter_f=iter_f,iter_acc=iter_acc)
  return(result)
}

#-----------------------------------------------------------------------------------------
#eset是行基因，列样本
biomarkerPick<-function(eset,label){
  #去掉低方差，填补缺失值
  zerovar<-nearZeroVar(t(eset))
  eset<-t(eset)[,-zerovar]
  imp<-preProcess(eset,method="knnImpute",k=5)
  eset<-predict(imp,eset)

  eset<-reduceGeneNum(eset)      #保留变异系数前5000的基因
  el<-removeOutliers(eset,label) #
  eset<-el$eset
  label<-el$label

  dissTOM<-buildNetwork(eset)
  moduleColors<-moduleDetect(eset,dissTOM)
  moduleTraitCor<-relateMT(eset,moduleColors,label)
  chooseModule<-chooseModuleByCor(moduleTraitCor)
  genes<-getFirstGeneSet(chooseModule) #筛选出第一批基因作为候选基因

}
#将每个模块中互信息最高的挑选出来合并
getFirstSet<-function(data,moduleColors,r=0.2){
  moduleNames<-unique(moduleColors)
  dataList<-sapply(moduleNames,function(x){
    data1<-data[,which(moduleColors==x)]
  })
}
#--聚类效果评估--
assessCluster<-function(clust1,clust2){
  if(length(clust1)!=length(clust2))
    stop("clust1的长度和cluster2不等..")
  if(is.character(clust1)|is.character(clust2)){
    clust1<-as.factor(clust1)
    clust2<-as.factor(clust2)
    clust1<-as.integer(clust1)
    clust2<-as.integer(clust2)
  }
  external.ind<-std.ext(clust1, clust2)
  result<-list()
  result$Rand<-clv.Rand(external.ind)
  result$Jaccard<-clv.Jaccard(external.ind)
  result$Folkes_Mallows<-clv.Folkes.Mallows(external.ind)
  result$Phi<-clv.Phi(external.ind)
  result$Russel_Rao<-clv.Russel.Rao(external.ind)
  return(result)
}
mcone<-function(eset,label,r){
  #将所有特征与label的相关性存在micFC中，将相关性高的(>r)的索引记下，存于subset中
  micFC<-apply(eset,2,function(x){m=mine(x,label);m$MIC})
  subset<-micFC[micFC>r] #将与label相关性高于r的选出
  subset<-sort(subset,decreasing = T) #降序排列
  names1<-names(subset)
  print("到了这里1")
  subset<-sapply(names1,function(x){which(colnames(eset)==x)}) #降序后，基因在原eset中的位置索引
  subset<-as.integer(subset)
  numSubset = length(subset)
  e=1
  while(e<length(subset)){
    q=e+1
    while(q<=length(subset)){
      if(mine(eset[,subset[e]],eset[,subset[q]])$MIC>=micFC[subset[q]])
        subset=subset[-q]
      else
        q=q+1
    }
    e=e+1
  }
  #return(eset[,subset])
  subset_c<-as.character(colnames(eset)[subset])
  return(subset_c)
}

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#m是字符向量

#路径以目标文件夹加/结尾
test2<-function(fileDir,eset,label){
  files<-list.files(fileDir)
  #获取所有文件的路径
  str1<-paste(fileDir,files,sep="")
  result1<-sapply(str1,function(x){m<-csvToEset(x);m<-as.character(m$name);list1<-trainModelNN(eset[,m],label,NULL);list1$result})
  return(result1)
}


createConfusionMatrix <- function(act, pred) {
  # You've mentioned that neither actual nor predicted may give a complete
  # picture of the available classes, hence:
  numClasses <- max(act, pred)
  # Sort predicted and actual as it simplifies what's next. You can make this
  # faster by storing `order(act)` in a temporary variable.
  pred <- pred[order(act)]
  act  <- act[order(act)]
  sapply(split(pred, act), tabulate, nbins=numClasses)
}
#验证用的方法，来看特征在模块中的位置
showCorPos<-function(eset,moduleColors,choose,label){
  colors<-unique(moduleColors)
  #获取choose所在模块颜色
  indexes<-match(choose,colnames(eset))
  choose_colors<-moduleColors[indexes]
  #获得各个模块表达谱子集,命名为color
  for(i in 1:length(colors)){
    str<-paste0(colors[i],'<-',substitute(eset),'[,which(moduleColors=="',colors[i],'")]')
    eval(parse(text=str))
  }
  print(colors)
  print('------------美丽的分割线--------------------')
  print(ls())
  #获取各个模块和label的cor的降序基因名,命名为color_dec
  for(i in 1:length(colors)){
    print(paste(colors[i],"标记1"))
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
    print(ls(pattern = "cor_*"))
  }
  print('-------------美丽分割线----------------')
  print(ls(pattern = "cor_*"))
  #获得各个模块的cor排名
  rank<-vector(length = length(choose),mode = "integer")
  for(i in 1:length(choose)){
    str<-paste0('rank[i]<-match("',choose[i],'",',choose_colors[i],'_dec)')
    eval(parse(text=str))
  }
  return(rank)

}
#验证用的方法，来看特征在模块中的位置,choose:features are known
wgcnaPredict<-function(eset,moduleColors,label){
  colors<-table(moduleColors)
  removeColors<-names(which(colors==1))
  colors<-setdiff(names(colors),removeColors)
  for(i in 1:length(colors)){
    str<-paste0(colors[i],'<-as.data.frame(eset[,which(moduleColors=="',colors[i],'")])')
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
  #将各个模块cor第一名放进去，颜色顺序与colors同
  first<-vector(length = length(colors),mode = "character")
  for(i in 1:length(colors)){
    #应该使用color_dec，以缩小范围，这里暂时全部遍历，不影响
    str<-paste0('first[i]<-colnames(cor_',colors[i],')[1]')
    eval(parse(text = str))
  }
  #迭代替换基因
  genelist<-first
  for(i in 1:length(colors)){
    str<-paste0('genelist<-replaceGene(genelist,',colors[i],'_dec,i,eset,label)')
    eval(parse(text = str))
  }
  return(genelist)
}

#replace geneVector genes with gene in coresponding module,to reach the highest predict score
replaceGene<-function(geneVector,color_dec,index,eset,label){
  if(length(color_dec)==1){
    return(geneVector)
  }
  acc<-trainModelNN(eset[,geneVector],as.factor(label))$result #记录初始精度
  acc_max<-acc
  max_index<-1
  geneVector2<-geneVector
  #replace gene with other gene which belong to same module,and pick the highest score gene and replace
  for(i in 2:length(color_dec)){
    geneVector2[index]<-color_dec[i] #replace it with genes in the same module
    acc2<-trainModelNN(eset[,geneVector2],as.factor(label))$result #记录初始精度
    if(acc2>acc_max){
      max_index<-i
      acc_max<-acc2
    }
  }
  geneVector[index]=color_dec[max_index]
  acc<-trainModelNN(eset[,geneVector],as.factor(label))$result #记录初始精度
  return(geneVector)
}
