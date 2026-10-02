library(arrow)
PanCancer_Shared_Matrix<-read_parquet("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/PanCancer_Shared/PanCancer_Shared_Features.parquet")
quantile(PanCancer_Shared_Matrix$PPI_Agg_Dep,probs=99/100,na.rm=T)
  99% 
0.999 
quantile(PanCancer_Shared_Matrix$PCom_Agg_Dep,probs=99/100,na.rm=T)
  99% 
0.998 
quantile(PanCancer_Shared_Matrix$GO_Sim,probs=99/100,na.rm=T)
> quantile(PanCancer_Shared_Matrix$SubLoc,probs=99/100,na.rm=T)
99% 
  1 
> quantile(PanCancer_Shared_Matrix$Domain,probs=99/100,na.rm=T)
99% 
  0 
> quantile(PanCancer_Shared_Matrix$Conservation_score,probs=99/100,na.rm=T)
99% 
199 
> quantile(PanCancer_Shared_Matrix$Gene_Age,probs=99/100,na.rm=T)
 99% 
2555 
> quantile(PanCancer_Shared_Matrix$PPI_Union,probs=99/100,na.rm=T)
99% 
842 

