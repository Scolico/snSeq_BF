targetSets<-agingSets

###calculate gsva score

matrix<-data@assays$RNA@counts
matrix<-NormalizeData(matrix,normalization.method='RC',scale.factor=1e6)
matrix<-Matrix(log2(matrix+1),sparse=T)

aa<-tapply(rownames(mmeta),mmeta$subtype,list)  
gsva.score.aging.terms<-tapply(rownames(mmeta),mmeta$subtype,function(a){
  m<-matrix[,colnames(matrix) %in% a]
  gsva(m,targetSets,min.sz=5,kcdf='Gaussian') ###用counts kcdf='Poisson',else'Gaussian'
})

gsva.score.aging.mean<-lapply(gsva.score.aging.terms, function(a){
  a.mmeta<-mmeta[rownames(mmeta) %in% colnames(a),]
  a<-t(a)
  by(a,a.mmeta$time.sp,colMeans)
})

gsva.score.aging.mean.matrix<-lapply(gsva.score.aging.mean,function(a){
  states<-names(a)
  enrichTerms<-names(a[[1]])
  x<-matrix(unlist(a),nrow=lengths(a),byrow=F)
  x<-data.frame(x)
  rownames(x)<-enrichTerms
  x$diif<-apply(x,1,function(c){
    max(c)-min(c)
  })
  colnames(x)<-c(paste0('states_',states),'diff')
  return(x)
})