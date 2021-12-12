library(tidyverse)
library(data.table)
library(dplyr)
data = fread('행정구역별,지목별_국토이용현황_시군구_20211014193259.csv')
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
corrplot(a,method="number",type = 'lower')

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
plot(1:7,f$sdev,type = 'b',main = 'Scree Plot')

secu_factanal <- factanal(z, factors = 2,rotation = "varimax", # "varimax", "promax", "none" 
                          scores="regression")
print(secu_factanal)
print(secu_factanal$loadings, cutoff=0)
secu_factanal$loadings

weight = c(0.740,0.706,0.761,0.751,0.484,-.039,0.829)

eco_demand = data.frame(t(z) * weight) %>% apply(2,sum)


sun = fread('경기도_태양광발전소설치현황.csv')
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

##
km = kmeans(q,centers = 3, nstart = 1, iter.max =30)
fviz_cluster(km,data=q,geom="point") + theme_bw() 
q = q %>% mutate(클러스터 = km$cluster)
q %>% filter(클러스터 == 3)



###AHP 위한 전처리
setwd('.\\공모전\\옥상햇빛')

library(raster)
library(rgeos)
library(rgdal)
library(maptools)

#경기도 shp데이터 
map = readOGR('LARD_ADM_SECT_SGG_경기/LARD_ADM_SECT_SGG_41.shp')
map = spTransform(map, CRSobj = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
korea = fortify(map, region='SGG_NM')

#shp 데이터 전처리
change = function(x){
  for (i in 1:length(x)){
    if(x[i] %in% c("성남시분당구", "성남시수정구", "성남시중원구")){
      x[i] = '성남시'
    }else if(x[i] %in% c("수원시권선구", "수원시영통구", "수원시장안구", "수원시팔달구")){
      x[i] = '수원시'
      
    }else if(x[i] %in% c("안산시단원구", "안산시상록구")){
      x[i] = '안산시'
      
    }else if(x[i] %in% c("기흥구", "수지구","처인구")){
      x[i] = '용인시'
      
    }else if(x[i] %in% c("덕양구", "일산동구","일산서구")){
      x[i] = '고양시'
      
    }else if(x[i] %in% c("안양시동안구", "안양시만안구")){
      x[i] = '안양시'
      
    }else{
      x[i] = x[i]
    }
    
  }
  return (x)
  
}


ads = korea %>% mutate(id = change(id)) %>% dplyr::mutate(위도 = round(lat,3),경도 = round(long,3)) %>% dplyr::select(c(경도,위도,id))
rm(map, korea)

#일사량 데이터 
dat = fread('태양광기상자원지도(년)_2010.csv',skip=4)


#일사량 데이터에 좌표기반으로 시군구 붙이기
fin = inner_join(ads,dat)
colnames(fin)[3] = '시군구명'
fin = fin[,c(1,2,3,7,8,9)]
ahp = fin %>% group_by(시군구명) %>% dplyr::summarise(전천일사량 =mean(누적전천일사량_지형))#,직달일사량 =mean(누적직달일사량_지형), 산란일사량 = mean(누적산란일사량_지형) )
rm(fin,dat,ads)

#미세먼지, 태양광발전잠재량,기온 데이터 붙이기
aa = fread('AHP분석용 데이터.csv',encoding='UTF-8')
ahp = ahp %>% filter(시군구명 %in% aa$시군구명)
ahp = inner_join(ahp,aa)
rm(aa)

#AHP분석을 위한 최종 데이터
ahp_final = ahp %>% dplyr::select(-시군구명) %>% scale()

#변수별 순위
ahp_rank = ahp %>% dplyr::select(-시군구명) 
for (i in 1:length(ahp_rank)){
  if (i==3){#미세먼지 낮을수록 청정도 높음
    ahp_rank[i] = rank(ahp_rank[i])
  }else{
    ahp_rank[i] = rank(-ahp_rank[i])
  }
  
}
rownames(ahp_rank) = as.vector(ahp$시군구명)

#2,3. 쌍대비교 및 가중치 산출
weight_list = list(전천일사량=matrix(0,12,12),잠재량=matrix(0,12,12), 미세먼지=matrix(0,12,12), 평균기온=matrix(0,12,12), 강수량=matrix(0,12,12), 평균풍속=matrix(0,12,12))

for (k in 1:length(ahp_rank)){
  for (i in 1:nrow(ahp_rank)){
    for (j in 1:nrow(ahp_rank)){
      diff = ahp_rank[[i,k]] - ahp_rank[[j,k]] #순위 뺐을때 마이너스면 앞에꺼가 더 우위->음수일때 큰 가중치
      weight = case_when(
        diff < -9 ~ 9,
        diff < -8 ~ 8,
        diff < -6 ~ 7,
        diff < -5.5 ~ 6,
        diff < -4.5 ~ 5,
        diff < -4 ~ 4,
        diff < -2 ~ 3,
        diff <= -1.5 ~ 2,
        abs(diff) < 1.5 ~ 1,
        diff <= 2 ~ 1/2,
        diff <= 4 ~ 1/3,
        diff <= 4.5 ~ 1/4,
        diff <= 5.5 ~ 1/5,
        diff <= 6 ~ 1/6,
        diff <= 8 ~ 1/7,
        diff <= 9 ~ 1/8,
        diff <= 11 ~ 1/9,
      )
      weight_list[[k]][i,j] = weight
      
    }
  }

}

#쌍대비교 행렬 이름변경
for(i in 1:length(weight_list)){
  row.names(weight_list[[i]]) = row.names(ahp_rank)
  colnames(weight_list[[i]]) = row.names(ahp_rank)
  
}

weight_list

#정규화 행렬
weight_list_normal = list(전천일사량=matrix(0,12,12),잠재량=matrix(0,12,12), 미세먼지=matrix(0,12,12), 평균기온=matrix(0,12,12), 강수량=matrix(0,12,12), 평균풍속=matrix(0,12,12))

for(i in 1:length(weight_list)){
  weight_list_normal[[i]] = t(t(weight_list[[i]])/colSums(weight_list[[i]]))
}
weight_list_normal

#항목별 가중치 산출
weight_cate=list()

for(i in 1:length(weight_list)){
  weight_cate[[i]]= apply(weight_list_normal[[i]], 1, mean)
}
weight_cate

sum(weight_cate[[1]])

#3.일관성 평가

#일관성 수치 계산
lambda_max=list()
for(i in 1:length(weight_list)){
  lambda_max[[i]] = mean((weight_list[[i]] %*% weight_cate[[i]])/weight_cate[[i]])
}
lambda_max

#일관성 지수(CI)
CI=list()
for(i in 1:length(weight_list)){
  CI[[i]] = (lambda_max[[i]] - 12)/11
}
CI

#RI (평균 무작위 수준(Saaty, 1982)을 기준으로 함, n=12일때)
RI = 1.48

#일관성 비율(CR)
CR=list()
for(i in 1:length(weight_list)){
  CR[[i]] = CI[[i]] / RI
}
CR #CR<0.1이므로 일관성이 존재한다고 판단.

names(CR) = names(weight_list_normal)
CR

#4.최종중요도 도출
weight_total = matrix(c(weight_cate[[1]],weight_cate[[2]],weight_cate[[3]], weight_cate[[4]], weight_cate[[5]],weight_cate[[6]]),nrow=12)
weight_total

rownames(weight_total) = rownames(ahp_rank)

#변수끼리 쌍대비교 (근거는 상관계수보고 내 임의로)
#변수별 상관관계 
library(corrplot)
corrplot(cor(ahp_final),method = 'number',type = 'lower',tl.srt=360)

#쌍대행렬
v = matrix(c(1,1/3,1,1,2,1,3,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1/2,1,1/2,1,1,1,1,1,1,1,1,1),nrow=6) #누가누구에게 영향주는지 + 상관계수
v
v_n = t(t(v)/colSums(v))
v_weight = apply(v_n,1,mean) #변수끼리 가중치
ri = 1.24 #n=6
lam = mean((v %*% v_weight)/v_weight)
ci = (lam - 6)/5
ci
ci/ri #0.1보다 작음
v_weight

#최종 중요도 -> 화성, 성남, 수원
t(weight_total %*% v_weight)

