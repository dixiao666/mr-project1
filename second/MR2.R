install.packages("pacman")
library("pacman")
install.packages("MendelianRandomization")
library("MendelianRandomization")
library(TwoSampleMR)
library("ieugwasr")
install.packages("purrr")
library("purrr")
library(readr) 
library(ggplot2)   
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMDI0MDA1MTA2NjhAc3htdS5lZHUuY24iLCJpYXQiOjE3NDM0MDg2NzgsImV4cCI6MTc0NDYxODI3OH0.h1IYMe9C4Z02uN5tapNEE9RDM8TYbVCNVZemiFl2AT-WvwboOYVjIc7tWETE2-g1njDNkIcCN-9EPJBQJpk9N9niTh3kx0c_6t80jf_nCjh8sRFoFmK2a8v_U_BD7hKUUt9F8fh0wMhH3A3rf0vcTyd88diKkN0-Egy01HJbsEQLjdxJ9HTrXzgfylkKP-hiBsJ02a6qlJPCBEBINI4lDhiNcV5JmpdT4VxUdT8KT305jFetkf2Ux3Tms49LD8hVI5FgUOIt20YukMlHEPyF_5vU2_z9Al02DdL3cu8Nrd6865YrErs2_DqgZPQla1eGbxM3nvUJ4UopGxR0WvMalA")
# 使用token
normalizePath("~/.Renviron")
ieugwasr::get_opengwas_jwt()
api_status()
user()
ieugwasr::user()
#暴露数据
expofile="ieu-a-300"
#染色体位置
chr_pos <- 1
#开始位置
pos_start <- 55505221
#结束位置
pos_end <- 55530525
#结局数据
outcfile <- "ieu-a-7"

#在线读取gwas数据
#提取检测的snp，存到biof_exposure_dat中
biof_exposure_dat <-extract_instruments(outcomes = expofile,
                                        clump = FALSE)

#进行连锁不平衡剔除
biof_exposure_dat <- clump_data(biof_exposure_dat,
                                clump_kb = 100,#定义连锁不平衡窗口大小
                                clump_r2 = 0.3,#定义连锁不平衡的R平方阈值
                                clump_p1 = 5e-08,#保留p值最小的snp
                                clump_p2 = 5e-08) #剔除p值大于阈值的snp
#获取去除连锁不平衡后的行数和列数
dim(biof_exposure_dat)
#subset函数从biof_exposure_dat筛选SNP
#chr.exposure == chr_pos: chromosome与药物靶点chromosome相等
#pos.exposure > pos_start - 100000:snp位置大于靶点start位置左侧100kb
#pos.exposure < pos_end + 100000:SNP位置小于靶点end位置右侧100kb
#这样提取出药物靶点周围±100kb的snp
Drug_Target_SNP <- subset(biof_exposure_dat,
                          chr.exposure == chr_pos &
                            pos.exposure > pos_start - 100000 &
                            pos.exposure < pos_end +100000)

# 将结果写入csv文件
# Drug_Target_SNP:提取的药物靶点周围snp
# Drug_Target_SNP.csv：输出的csv文件名
write.csv(Drug_Target_SNP, file = "Drug_Target_SNP.csv")
#从outcome数据中提取与药物靶点SNP相关的表型数据
biof_Outcome_dat <- extract_outcome_data(
  snps =  Drug_Target_SNP$SNP,
  outcomes = outcfile
)
#harmonise and merge 数据
#确保SNP对暴露和结果的效应基于同一等位基因
harmonise_dat <- harmonise_data(
  exposure_dat = Drug_Target_SNP ,
  outcome_dat = biof_Outcome_dat ,
  action = 1)
#查看harmonise后的数据
View(harmonise_dat)

#将结果写入文件
write.table(harmonise_dat,
            "expo_outc_harmonise.csv",
            row.names = F,
            sep="\t",
            quote = F)

#去除混杂因素
p_load(MendelianRandomization,purrr,reader)
#设置变量grp_size,用于设置每个分组的SNP数量
grp_size <- 100
#将data中的SNP分成多个大小为grp_size的分组，保存到变量grps中
grps <- split(harmonise_dat$SNP,ceiling(seq_along(harmonise_dat$SNP)/grp_size))
#通过phenoscanner API对每个分组进行关联性分析
results <- map_dfr(grps, ~phenoscanner(snpquery = .x,
                                       catalogue = "GWAS",
                                       pvalue = 1e-05,
                                       proxies = "None",
                                       r2 = 0.8,
                                       build = 37) $results)




#将关联结果写入文件confounder.csv
write_csv(results,"confounder.csv")
View(rusults)
#设置变量remove_snps,用于保存要移除的混杂因素snp
remove_snps <- c("rs764406","rs2456973","rs9739070")
#移除含有混杂SNP的行
harmonise_dat <- harmonise_dat[!harmonise_dat$SNP %in% remove_snps,]
#将移除混杂因素后的结果写入clean_confounder07.csv
write_csv(harmonise_dat, "clean_confounder.csv")

res=mr(harmonise_dat)
res
write.csv(res, file="MRres.csv")
#mendelian randomization(MR) analysis
MR_result <- mr(harmonise_dat)
View(MR_result)

#view(MR_result)
result_OR=generate_odds_ratios(MR_result)
result_OR$or1=1/result_OR$or
result_OR$or_lci951=result_OR$or1-1.96*result_OR$se
result_OR$or_uci951=result_OR$or1+1.96*result_OR$se
write.table(result_OR[,5:ncol(result_OR)],"MR_OR.xls",row.names = F,sep="\t",quote = F)

mr_heterogeneity(harmonise_dat, method_list = c("mr_egger_regression","mr_ivw"))#异质性检验
#outlier test
run_mr_presso(harmonise_dat, NbDistribution = 1000)#偏倚检验
pleio <- mr_pleiotropy_test(harmonise_dat)#多效性——MR egger
#view (pleio)
write.csv(pleio, file="pleio.csv")#多效性结果
single <- mr_leaveoneout(harmonise_dat)
mr_leaveoneout_plot(single)#留一法检验敏感性
write.csv(single, file="single.csv")

p1 <-mr_scatter_plot(res,harmonise_dat)#散点图
p1
ggsave(p1[[1]],file="scatter.pdf",width=8,height=8)

#单snp分析
res_single=mr_singlesnp(harmonise_dat)
#绘制森林图，森林图可以表示每个snp与结果的关联效应
p2 <- mr_forest_plot(res_single)
p2[1]
ggsave(p2[[1]],file="forest.pdf",width=8,height=8)

#留一法敏感性分析
p4<- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonise_dat))
p4[[1]]
ggsave(p4[1], file="leaveoneout.pdf",width=8,height=8)







harmonise_dat$maf.exposure <- ifelse(
  harmonise_dat$eaf.exposure > 0.5, 
  1 - harmonise_dat$eaf.exposure, 
  harmonise_dat$eaf.exposure
  )

# 清理无效值
harmonise_dat_clean <- harmonise_dat %>%
  filter(
    !is.na(maf.exposure),
    maf.exposure > 0, 
    maf.exposure <= 0.5
    )
result <- coloc.abf(dataset1=list(pvalues=harmonise_dat_clean$pval.outcome, snp=harmonise_dat_clean$SNP, type="cc", s=0.33, N=50000), 
                    dataset2=list(pvalues=harmonise_dat_clean$pval.exposure, snp=harmonise_dat_clean$SNP, type="quant", N=10000), MAF=harmonise_dat_clean$maf.exposure )
library(dplyr)
need_result=result$results %>% filter(SNP.PP.H4 > 0.95)




