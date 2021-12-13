library(tidyverse)
library(data.table)
library(dplyr)
data = fread('C:\\Users\\User\\Desktop\\데이터 안심구역 데이터 활용 아이디어&시각화 경진대회\\QGIS\\클터.csv')
str(data)

is.integer64 <- function(x){
  class(x)=="integer64"
}

data = data %>% mutate_if(is.integer64,as.numeric)

data_scale = data %>% select(-c(시군구)) %>% scale()
data_scale = data.frame(data_scale)


model = data_scale %>% lm(온실가스배출량~.,data=.)
summary(model)
a = data_scale %>% select(-온실가스배출량) %>% cor()

library(corrplot)
corrplot(a,method="number")

library(car)
vif(model)

library(gvlma)
gvlma(model)

library(lmtest)
dwtest(model)

par(mfrow=c(2,2))
plot(model) 

coef(model)
cor(data_scale)


library(psych)
KMO(r=cor(data_scale))

z = data_scale %>% select(-온실가스배출량)
KMO(z)
f = prcomp(z)
#2개
plot(1:7,f$sdev,type = 'b')

secu_factanal <- factanal(z, factors = 2,rotation = "varimax", # "varimax", "promax", "none" 
                          scores="regression")
print(secu_factanal)
print(secu_factanal$loadings, cutoff=0)
secu_factanal$loadings

weight = c(0.740,0.706,0.761,0.751,0.484,-.039,0.829)

eco_demand = data.frame(t(z) * weight) %>% apply(2,sum)


sun = fread('C:\\Users\\User\\Desktop\\데이터 안심구역 데이터 활용 아이디어&시각화 경진대회\\QGIS\\경기도_태양광발전소설치현황.csv')
str(sun)
y = sun %>% group_by(시군명) %>% summarise(n=n())
j = c(121,142,268,591,53,156,58,458,698,458,82,140,674,39,277,36,431,461,
      310,43,33,553,827,54,93,608,878,96,36,844,676)
q = cbind(y,j,eco_demand)
q = q %>% mutate(면적당태양광 = n/j) %>% select(-c(n,j))
row.names(q) = q$시군명
q = q %>% select(-시군명)

library(caret)
library(corrplot)
library(cluster)
library(factoextra)
library(gridExtra)

plot1 = fviz_nbclust(q,kmeans,method='wss')
plot2 = fviz_nbclust(q,kmeans,method='silhouette')
grid.arrange(plot1,plot2,ncol=2)

plot(q$eco_demand,q$면적당태양광)

##
km = kmeans(q,centers = 3, nstart = 1, iter.max =30)
fviz_cluster(km,data=q,geom="point") + theme_bw()
q = q %>% mutate(클러스터 = km$cluster)
q %>% filter(클러스터 == 1)
gap_stat <- clusGap(q, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)



##k-medoids
plot3 = fviz_nbclust(q, pam, method = "wss")
plot4 = fviz_nbclust(q, pam, method = 'silhouette')
grid.arrange(plot3,plot4,ncol=2)


gap_stat <- clusGap(q,
                    FUN = pam,
                    K.max = 10, #max clusters to consider
                    B = 50) #total bootstrapped iterations

#plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
kmed <- pam(q, k = 3)
fviz_cluster(kmed, data = q)+ theme_bw()


#install.packages("fpc")
library(fpc)
pamk.result <- pamk(q)
pamk.result$nc
plot(pamk.result$pamobject)

##
ds <- dbscan(q, eps=0.8, MinPts=5)
plotcluster(q, ds$cluster)
fviz_cluster(ds, data = q)+ theme_bw()
plot1 = fviz_nbclust(q,dbscan,method='wss')
plot2 = fviz_nbclust(q,dbscan,method='silhouette')
grid.arrange(plot1,plot2,ncol=2)

##Hierarchical 클러스터링
hc <- hclust(dist(q), method="ave")
plot(hc, hang=-1)
